#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <float.h>

const double ln0 = -10000;

typedef struct {
  double ln;
  int sgn;
} lnum_t;

__inline static void lnum_set(lnum_t *x, double num)
{
  if ( num > 0 ) {
    x->ln = log(num);
    x->sgn = 1;
  } else if ( num < 0 ) {
    x->ln = log(-num);
    x->sgn = -1;
  } else {
    x->ln = ln0;
    x->sgn = 1;
  }
}

__inline static double lnum_get(lnum_t *x)
{
  return x->sgn * exp(x->ln);
}

__inline static lnum_t *lnum_copy(lnum_t *y, lnum_t *x)
{
  y->ln = x->ln;
  y->sgn = x->sgn;
  return y;
}

__inline static lnum_t *lnum_add(lnum_t *z, lnum_t *x, lnum_t *y)
{
  int sgn = x->sgn * y->sgn;
  double lnx = x->ln, lny = y->ln, del = lny - lnx;

  if ( fabs(del) < 1e-10 && sgn < 0 ) {
    z->sgn = 1;
    z->ln = ln0;
  } else if ( del < 0 ) {
    z->sgn = x->sgn;
    z->ln = lnx + log(1 + sgn * exp(del));
  } else {
    z->sgn = y->sgn;
    z->ln = lny + log(1 + sgn * exp(-del));
  }
  return z;
}

__inline static lnum_t *lnum_mul(lnum_t *z, lnum_t *x, lnum_t *y)
{
  z->sgn = x->sgn * y->sgn;
  z->ln = x->ln + y->ln;
  return z;
}



const int blksz = 64;

typedef struct {
  int n, ncap;
  lnum_t *a; /* coefficients */
} lpoly_t;

static lpoly_t *lpoly_open(void)
{
  lpoly_t *p;
  int i;

  p = calloc(1, sizeof(*p));
  p->n = 0;
  p->ncap = blksz;
  p->a = calloc(p->ncap, sizeof(lnum_t));
  for ( i = 0; i < p->ncap; i++ ) {
    lnum_set(p->a + i, 0);
  }
  return p;
}

static void lpoly_close(lpoly_t *p)
{
  free(p->a);
  free(p);
}

static void lpoly_print(lpoly_t *p)
{
  int i, n = p->n;
  double x;

  for ( i = 0; i < n; i++ ) {
    if ( p->a[i].ln < -10 ) continue;
    x = lnum_get(p->a + i);
    if ( i == 0 ) printf("%g", x);
    else if ( i == 1 ) printf("%+gx", x);
    else printf("%+gx^%d", x, i);
  }
  printf("\n");
}

/* resize and clear content */
static lpoly_t *lpoly_resize(lpoly_t *p, int n)
{
  int i;

  if ( n > p->ncap ) {
    p->ncap = ((n / blksz) + 1) * blksz;
    p->a = realloc(p->a, p->ncap * sizeof(lnum_t));
  }
  p->n = n;
  for ( i = 0; i < p->ncap; i++ )
    lnum_set(p->a + i, 0);
  return p;
}

/* q = p */
static lpoly_t *lpoly_copy(lpoly_t *q, const lpoly_t *p)
{
  int i, n = p->n;

  lpoly_resize(q, n);
  for ( i = 0; i < n; i++ )
    lnum_copy(q->a + i, p->a + i);
  return q;
}

/* p = a[0] + a[1]*x + a[2]*x^2 + ... */
static lpoly_t *lpoly_set(lpoly_t *p, int n, ...)
{
  int i;
  va_list vl;
  double x;

  lpoly_resize(p, n);
  va_start(vl, n);
  for ( i = 0; i < n; i++ ) {
    x = va_arg(vl, double);
    lnum_set(p->a + i, x);
  }
  va_end(vl);
  return p;
}

/* p = q * s */
static void lpoly_mullnum(lpoly_t *p, lpoly_t *q, lnum_t *ls)
{
  int i, n = q->n;

  lpoly_resize(p, n);
  for ( i = 0; i < n; i++ ) {
    p->a[i].sgn = q->a[i].sgn * ls->sgn;
    p->a[i].ln = q->a[i].ln + ls->ln;
  }
}

/* p *= s */
static void lpoly_imulnum(lpoly_t *p, double s)
{
  int i;
  lnum_t ls;

  lnum_set(&ls, s);
  for ( i = 0; i < p->n; i++ ) {
    p->a[i].sgn *= ls.sgn;
    p->a[i].ln += ls.ln;
  }
}

/* r = p + q */
static void lpoly_add(lpoly_t *r, lpoly_t *p, lpoly_t *q)
{
  int i, n = p->n, m = q->n, nmax = ( n > m ) ? n : m;

  lpoly_resize(r, nmax);
  for ( i = 0; i < nmax; i++ ) {
    if ( i < m ) {
      if ( i < n ) {
        lnum_add(r->a + i, p->a + i, q->a + i);
      } else {
        lnum_copy(r->a + i, q->a + i);
      }
    } else if ( i < n ) {
      lnum_copy(r->a + i, p->a + i);
    }
  }
}

/* r = p + s * q */
static void lpoly_sadd(lpoly_t *r, lpoly_t *p, lpoly_t *q, double s)
{
  int i, n = p->n, m = q->n, nmax = ( n > m ) ? n : m;
  lnum_t ls[1], t[1];

  lpoly_resize(r, nmax);
  lnum_set(ls, s);
  for ( i = 0; i < nmax; i++ ) {
    if ( i < m ) {
      lnum_mul(t, q->a + i, ls);
      if ( i < n ) {
        lnum_add(r->a + i, p->a + i, t);
      } else {
        lnum_copy(r->a + i, t);
      }
    } else if ( i < n ) {
      lnum_copy(r->a + i, p->a + i);
    }
  }
}

/* z = x * y */
static lpoly_t *lpoly_mul(lpoly_t *z, const lpoly_t *x, const lpoly_t *y)
{
  int i, j, n = x->n, m = y->n;
  lnum_t t[1], s[1];

  lpoly_resize(z, n + m - 1);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < m; j++ ) {
      lnum_mul(t, x->a + i, y->a + j);
      lnum_add(s, z->a + i + j, t);
      lnum_copy(z->a + i + j, s);
    }
  }
  return z;
}

/* y *= x */
static lpoly_t *lpoly_imul(lpoly_t *y, const lpoly_t *x, lpoly_t *t)
{
  lpoly_mul(t, y, x);
  lpoly_copy(y, t);
  return y;
}

/* y = x^n */
static lpoly_t *lpoly_pow(lpoly_t *y, const lpoly_t *x, int n, lpoly_t *t)
{
  int i;
  lpoly_copy(y, x);
  for ( i = 1; i < n; i++ )
    lpoly_imul(y, x, t);
  return y;
}

static lpoly_t *is2dos(int n, int m)
{
  lpoly_t *beta, *bm, *p1, *p2, *p3, *t, *ak, *x1p, *x1m;
  lpoly_t *c0, *s0, *cn, *sn, **csqr, **ssqr;
  lpoly_t **abpow, **akpow, *Z1, *Z2, *Z3, *Z4, *Z;
  lnum_t mfac[1], pow2m[1], tmp[1];
  int i, j, k = 1;

  beta = lpoly_open();
  bm = lpoly_open();
  ak = lpoly_open();
  p1 = lpoly_open();
  p2 = lpoly_open();
  p3 = lpoly_open();
  t = lpoly_open();
  x1p = lpoly_open();
  x1m = lpoly_open();
  c0 = lpoly_open();
  s0 = lpoly_open();
  cn = lpoly_open();
  sn = lpoly_open();
  Z1 = lpoly_open();
  Z2 = lpoly_open();
  Z3 = lpoly_open();
  Z4 = lpoly_open();
  Z = lpoly_open();
  csqr = calloc(n + 1, sizeof(*csqr));
  ssqr = calloc(n + 1, sizeof(*ssqr));
  abpow = calloc(m + 1, sizeof(*abpow));
  akpow = calloc(m + 1, sizeof(*akpow));
  for ( i = 0; i <= n; i++ ) {
    csqr[i] = lpoly_open();
    ssqr[i] = lpoly_open();
  }
  for ( i = 0; i <= m; i++ ) {
    abpow[i] = lpoly_open();
    akpow[i] = lpoly_open();
  }
  lpoly_set(beta, 4, 0.0, 2.0, 0.0, -2.0); /* 2x(1-x^2) */
  lpoly_pow(p1, beta, m, t); /* beta^m */
  lnum_set(pow2m, 1.0);
  pow2m->ln = (m-1)*log(0.5); /* 0.5^(m-1) */
  lpoly_mullnum(bm, p1, pow2m);

  lnum_set(mfac, 1.0);
  for ( i = 2; i <= m; i++ )
    mfac->ln += log(i*0.5);

  lpoly_set(p1, 2, 1.0, 1.0);
  lpoly_pow(x1p, p1, m, t); /* x1p = (1+x)^m */
  lpoly_set(p1, 2, 1.0, -1.0);
  lpoly_pow(x1m, p1, m, t); /* x1m = (1-x)^m */
  lpoly_resize(p1, m + 1);
  lnum_set(p1->a + m, 1); /* p1 = x^m */
  lpoly_mul(p2, p1, x1p); /* p2 = x^m (1 + x)^m */
  lpoly_mul(p3, p1, x1m); /* p3 = x^m (1 - x)^m */
  lpoly_add(c0, x1m, p2);
  lpoly_sadd(s0, x1m, p2, -1); /* s0 = (1-x)^m - x^m (1+x)^m */
  lpoly_add(cn, x1p, p3);
  lpoly_sadd(sn, x1p, p3, -1); /* sn = (1+x)^m - x^m (1-x)^m */

  for ( k = 1; k < n; k++ ) {
    lnum_set(tmp, -cos(M_PI*k/n));
    lpoly_mullnum(p1, beta, tmp); /* p1 = beta*cos(pi*k/n) */
    lpoly_set(p2, 5, 1.0, 0.0, 2.0, 0.0, 1.0); /* p2 = (1+x^2)^2 */
    lpoly_add(ak, p2, p1); /* ak = (1+x^2)^2 - beta*cos(pi*k/n) */

    /* abpow[j] = (ak^2 - beta^2)^j/(2j)! */
    lpoly_pow(p1, ak, 2, t);
    lpoly_pow(p2, beta, 2, t);
    lpoly_sadd(p3, p1, p2, -1); /* p3 = ak^2 - beta^2 */
    lpoly_set(abpow[0], 1, 1.0);
    for ( j = 1; j <= m/2; j++ ) {
      lpoly_mul(abpow[j], abpow[j-1], p3);
      lpoly_imulnum(abpow[j], 1.0/(2.*j*(2*j-1)));
    }
    /* akpow[j] = ak^(m-2j)m!/(m-2j)!/2^(m-1) */
    lpoly_set(akpow[0], 1, 1.0);
    lnum_copy(akpow[0]->a, mfac);
    for ( j = 1; j <= m; j++ ) {
      lpoly_mul(akpow[j], akpow[j-1], ak);
      lpoly_imulnum(akpow[j], 1./j);
    }

    /* ck^2 = bm + Sum_{j=0 to m/2} abpow[j]*akpow[m-2j] */
    lpoly_copy(p2, bm);
    for ( j = 0; j <= m/2; j++ ) {
      lpoly_mul(p1, abpow[j], akpow[m-2*j]);
      lpoly_add(p3, p2, p1);
      lpoly_copy(p2, p3);
    }
    lpoly_copy(csqr[k], p2);
    /* sk^2 = ck^2 - 2*(beta^m/2^(m-1)) */
    lpoly_sadd(ssqr[k], csqr[k], bm, -2.0);
  }

  if ( n % 2 == 0 ) { /* n is even */
    lpoly_set(Z1, 1, 1.0);
    lpoly_set(Z2, 1, 1.0);
    lpoly_mul(Z3, c0, cn);
    lpoly_mul(Z4, s0, sn);
    for ( k = 0; k < n/2; k++ ) {
      lpoly_imul(Z1, csqr[2*k+1], p1);
      lpoly_imul(Z2, ssqr[2*k+1], p1);
      if ( k > 0 ) {
        lpoly_imul(Z3, csqr[2*k], p1);
        lpoly_imul(Z4, ssqr[2*k], p1);
      }
    }
  } else { /* n is odd */
    lpoly_copy(Z1, cn);
    lpoly_copy(Z2, sn);
    lpoly_copy(Z3, c0);
    lpoly_copy(Z4, s0);
    for ( k = 0; k < (n-1)/2; k++ ) {
      lpoly_imul(Z1, csqr[2*k+1], p1);
      lpoly_imul(Z2, ssqr[2*k+1], p1);
      lpoly_imul(Z3, csqr[2*k+2], p1);
      lpoly_imul(Z4, ssqr[2*k+2], p1);
    }
  }
  lpoly_add(p1, Z1, Z2);
  lpoly_add(p2, p1, Z3);
  lpoly_add(Z, p2, Z4);
  lpoly_imulnum(Z, 0.5);

  lpoly_close(beta);
  lpoly_close(bm);
  lpoly_close(p1);
  lpoly_close(p2);
  lpoly_close(p3);
  lpoly_close(t);
  lpoly_close(ak);
  lpoly_close(x1p);
  lpoly_close(x1m);
  lpoly_close(c0);
  lpoly_close(s0);
  lpoly_close(cn);
  lpoly_close(sn);
  lpoly_close(Z1);
  lpoly_close(Z2);
  lpoly_close(Z3);
  lpoly_close(Z4);
  for ( i = 0; i <= n; i++ ) {
    lpoly_close(csqr[i]);
    lpoly_close(ssqr[i]);
  }
  for ( i = 0; i <= m; i++ ) {
    lpoly_close(abpow[i]);
    lpoly_close(akpow[i]);
  }
  free(csqr);
  free(ssqr);
  free(abpow);
  free(akpow);
  return Z;
}

static int is2dos_save(lpoly_t *p, int n, int m)
{
  char fn[32];
  FILE *fp;
  int i, digs;
  double x, lnimax;

  sprintf(fn, "is2lndos%dx%d.dat", n, m);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  lnimax = log((double) ULONG_MAX);
  digs = (int)(-log10(p->a[p->n/2].ln*DBL_EPSILON));
  if ( digs < 0 ) digs = 0;
  printf("digs %d\n", digs);
  for ( i = 0; i < p->n; i+= 2 ) {
    x = p->a[i].ln;
    if ( x < 0 ) x = ln0;
    else if ( x < lnimax ) /* round to nearest integer */
      x = log((unsigned long) (exp(x)+.5));
    fprintf(fp, "%d %.*f\n", -2*n*m+2*i, digs, x);
  }
  fclose(fp);
  return 0;
}

int main(int argc, char **argv)
{
  int n = 4, m = 4;
  lpoly_t *p;

  if ( argc == 2 ) {
    n = m = atoi(argv[1]);
  } else if ( argc == 3 ) {
    n = atoi(argv[1]);
    m = atoi(argv[2]);
  }
  p = is2dos(n, m);
  if ( n*m <= 32 ) lpoly_print(p);
  is2dos_save(p, n, m);
  lpoly_close(p);
  return 0;
}

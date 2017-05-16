#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <float.h>

const double ln0 = -10000;

/* log(exp(a) + exp(b)) */
__inline static double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = b - a) < ln0) ? a : a + log(1 + exp(c));
}

/* log(exp(a) - exp(b)), only works for a > b */
__inline static double lndif(double a, double b)
{
  double c;
  return ((c = b - a) < ln0) ? a : a + log(1 - exp(c));
}

/* log(exp(a)+b) */
__inline static double lnaddn(double a, double b)
{
  return (a > -ln0) ? a : a + log(1 + b*exp(-a));
}

/* the exact logarithm of the partition function, Z
 * the average energy and heat capacity of the Ising model */
__inline static double is2exact(int n, int m, double beta, double *eav, double *cv)
{
  double mh, nm, ex, f, th, sech, bet2, bsqr, log2, x;
  double lnz, lnz1, lnz2, lnz3, lnz4, dz, ddz;
  double z21, z31, z41, za1;
  double dr1, dr2, dr3, dr4, ddr1, ddr2, ddr3, ddr4;
  double g, g0, dg, ddg, dg0;
  double xn2b, sh2b, coth2b;
  double lnch2b, lncc2b, lncl, lnsl, cd, cdsqr, lnddcl;
  int k, sgn4 = 1;

  mh = .5*m;
  nm = m * n;
  log2 = log(2.0);
  bet2 = 2.*beta;
  bsqr = beta*beta;
  xn2b = exp(-bet2);

  lnz1 = lnz2 = lnz3 = lnz4 = 0;
  dr1 = dr2 = dr3 = dr4 = 0;
  ddr1 = ddr2 = ddr3 = ddr4 = 0;
  lnch2b = lnadd(bet2, -bet2) - log2;
  coth2b = 2./(1. - xn2b*xn2b) - 1.;
  lncc2b = lnch2b + log(coth2b); /* ln[ cosh(2b) * coth(2b) ] */
  g0 = bet2 + log(2./(1. + xn2b) - 1.);
  sgn4 = (g0 >= 0) ? 1 : -1;

  sh2b = 0.5*(1./xn2b - xn2b);
  dg0 = 2. + 2./sh2b;
  x = sh2b*sh2b;
  cd = 2. - 2./x; /* cl' = cd * cosh(2b) */
  cdsqr = cd*cd;
  lnddcl = lnaddn(lncc2b, 2.0/(x * sh2b)) + 2.*log2; /* log(cl'') */

  for ( k = 0; k < n; k++ ) { /* for odd number */
    /* cosh(gamma_k) = cosh^2(2K)/sinh(2K) - cos(pi*k/n) */
    lncl = lnaddn(lncc2b, -cos((2*k + 1)*M_PI/n));
    lnsl = lncl + 0.5*log(1. - exp(-2.*lncl));
    g = lnadd(lncl, lnsl); /* gamma_{2k+1} */
    f = mh * g;
    lnz1 += lnadd(f, -f);
    lnz2 += lndif(f, -f);

    dg = exp(lnch2b - lnsl)*cd; /* g' = cl'/sl; */
    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    x = mh*dg;
    dr1 += x*th;
    dr2 += x/th;
    printf("dg %g, th %g\n", dg, th);

    /* g''=cl''/sl - cl' ^2 *cl/sl^3; */
    ddg = exp(lnddcl - lnsl);
    ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    sech = 2.0*dg/(ex + 1.0/ex); /* g' * sech(0.5*m*g) */
    ddr1 += mh*(ddg*th + mh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex); /* g' * csch(0.5*m*g) */
    ddr2 += mh*(ddg/th - mh*(sech*sech));

    if ( k == 0 ) {
      g = g0;
    } else {
      lncl = lnaddn(lncc2b, -cos(2.0*M_PI*k/n));
      lnsl = lncl+0.5*log(1-exp(-2*lncl));
      g = lnadd(lncl, lnsl);
    }
    f = mh*g;
    lnz3 += lnadd(f, -f); /* log [2 cosh(f)] */
    lnz4 += (f < 0) ? lndif(-f, f) : lndif(f, -f); /* avoid neg. g0 */

    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    dg = (k == 0) ? dg0 : exp(lnch2b - lnsl)*cd;
    dr3 += mh*dg*th;
    dr4 += mh*dg/th;

    if ( k == 0 ) {
      ddg = -4*coth2b*coth2b*exp(-lnch2b);
    } else {
      ddg = exp(lnddcl - lnsl);
      ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    }
    sech = 2.0*dg/(ex + 1.0/ex);
    ddr3 += mh*(ddg*th + mh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex);
    ddr4 += mh*(ddg/th - mh*(sech*sech));
  }

  z21 = exp(lnz2 - lnz1);
  z31 = exp(lnz3 - lnz1);
  z41 = sgn4*exp(lnz4 - lnz1);
  za1 = 1.0 + z21 + z31 + z41;
  lnz = lnz1 + log(za1);
  lnz += .5*nm*log(2.*sh2b) - log2;
  dz = (dr1 + z21*dr2 + z31*dr3 + z41*dr4)/za1;
  printf("%g %g %g; %g %g %g %g\n",
      z21, z31, z41, dr1, dr2, dr3, dr4);
  if ( eav != NULL ) {
    *eav = - nm*coth2b - dz;
  }
  ddr1 += dr1 * dr1;
  ddr2 += dr2 * dr2;
  ddr3 += dr3 * dr3;
  ddr4 += dr4 * dr4;
  ddz = (ddr1 + z21*ddr2 + z31*ddr3 + z41*ddr4)/za1;
  if ( cv != NULL ) {
    *cv = bsqr * (-2.*nm/(sh2b*sh2b) + ddz - dz*dz);
  }
  return lnz;
}


typedef struct {
  double ln;
  int sgn;
} lnum_t;

__inline static lnum_t *lnum_set(lnum_t *x, double num)
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
  return x;
}

__inline static lnum_t *lnum_setln(lnum_t *x, double ln)
{
  x->ln = ln;
  x->sgn = 1;
  return x;
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

__inline static lnum_t *lnum_sub(lnum_t *z, lnum_t *x, lnum_t *y)
{
  lnum_t ny = {y->ln, -y->sgn};
  return lnum_add(z, x, &ny);
}

__inline static lnum_t *lnum_mul(lnum_t *z, lnum_t *x, lnum_t *y)
{
  z->sgn = x->sgn * y->sgn;
  z->ln = x->ln + y->ln;
  return z;
}

__inline static void lnum_imul(lnum_t *y, lnum_t *x)
{
  lnum_t z[1];
  lnum_mul(z, y, x);
  lnum_copy(y, z);
}

__inline static lnum_t *lnum_div(lnum_t *z, lnum_t *x, lnum_t *y)
{
  z->sgn = x->sgn * y->sgn;
  z->ln = x->ln - y->ln;
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
static void lpoly_resize(lpoly_t *p, int n)
{
  int i;

  if ( n > p->ncap ) {
    p->ncap = ((n / blksz) + 1) * blksz;
    p->a = realloc(p->a, p->ncap * sizeof(lnum_t));
  }
  p->n = n;
  for ( i = 0; i < p->ncap; i++ )
    lnum_set(p->a + i, 0);
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
static void lpoly_set(lpoly_t *p, int n, ...)
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
  for ( i = 0; i < p->n; i+= 2 ) {
    x = p->a[i].ln;
    if ( x < 0 ) x = ln0;
    else if ( x < lnimax ) /* round to nearest integer */
      x = log((double)((unsigned long) (exp(x)+.5)));
    fprintf(fp, "%d %.*f\n", -2*n*m+2*i, digs, x);
  }
  fclose(fp);
  return 0;
}


__inline static double is2_exact(int n, int m, double beta, double *eav, double *cv)
{
  lnum_t t1[1], t2[1], t3[1], t4[1], sh[1], ccs[1];
  lnum_t dg[1], ddg[1], *cg, *sg;
  lnum_t Y[1], Y1[1], Y2[1], Y3[1], Y4[1];
  double xp, gamk, *dgam, *ddgam, *tg, log2, ish2, r2, r3, r4;
  double dlnY1, dlnY2, dlnY3, dlnY4;
  double ddlnY1, ddlnY2, ddlnY3, ddlnY4;
  int k;

  log2 = log(2);
  /* sh = sinh(2*K) = exp(2*beta)(1-exp(-4*beta)/2); */
  xp = exp(-2*beta);
  lnum_set(t1, (1 - xp*xp)/2);
  lnum_imul(lnum_setln(sh, 2*beta), t1);
  lnum_copy(t3, sh);
  t3->ln = -t3->ln; /* t3 = 1/sinh(2*K); */
  lnum_add(ccs, sh, t3); /* ccs = sinh(2*K) + 1/sinh(2*K) */
  ish2 = exp(2*t3->ln); /* 1/sinh^2(2*K) */
  /* dg = ccs' = 2*cosh(2*K)*(1 - 1/sinh^2(2*K)) */
  lnum_set(t4, (1 + xp*xp)*(1 - ish2));
  lnum_imul(lnum_setln(dg, 2*beta), t4);
  /* ddg = ccs'' */
  lnum_setln(ddg, sh->ln + log(4*(1 + ish2 + ish2*ish2)));
  cg = calloc(2*n, sizeof(*cg));
  sg = calloc(2*n, sizeof(*sg));
  tg = calloc(2*n, sizeof(*tg));
  dgam = calloc(2*n, sizeof(*dgam));
  ddgam = calloc(2*n, sizeof(*ddgam));
  for ( k = 0; k < n * 2; k++ ) {
    if ( k == 0 ) {
      gamk = 2*beta + log((1 - xp)/(1 + xp));
      dgam[k] = 2 + 4*xp/(1 - xp*xp);
    } else {
      lnum_add(t2, ccs, lnum_set(t4, -cos(M_PI*k/n))); /* t2 = cosh(gamk) */
      /* gamk = arccosh(t2) = ln(t2 + sqrt(t2^2-1)) = ln(t2) + ln(1+sqrt(1-t2^(-2))) */
      xp = sqrt(1 - exp(-2*t2->ln));
      gamk = t2->ln + log(1 + xp);
      /* ln sinh(gamk) = ln(t2) + ln sqrt(1 - t2^(-2)); */
      t3->sgn = t2->sgn;
      t3->ln = t2->ln + log(xp);
      dgam[k] = lnum_get( lnum_div(t1, dg, t3) ); /* dg/sinh(gamk) */
      lnum_imul(lnum_set(t3, -dgam[k]*dgam[k]), t2);
      lnum_add(t4, ddg, t3); /* t4 = ddg - cosh(gamk)*dgam[k]^2 */
      ddgam[k] = t4->sgn * exp(t4->ln - t2->ln)/xp;
      printf("k %d, dgam %g, ddgam %g\n", k, dgam[k], ddgam[k]);
      /* (exp(ddg)-cosh(gammak)*exp(2*dgam[k]))/sinh(gammak) */
    }
    t1->sgn = t2->sgn = 1;
    t1->ln = m * gamk * 0.5; /* t1 = exp(m*gamk/2) */
    t2->ln = -t1->ln; /* t2 = exp(-m*gamk/2) */
    lnum_add(cg + k, t1, t2); /* cg[k] = 2 cosh(m*gamk/2) */
    lnum_sub(sg + k, t1, t2); /* sg[k] = 2 sinh(m*gamk/2) */
    tg[k] = lnum_get( lnum_div(t3, sg + k, cg + k) );
    //printf("k %d, gamk %g, dgamk %g, c %g, s %g, tanh %g\n",
    //    k, gamk, dgam[k], exp(cg[k].ln), exp(sg[k].ln), exp(sg[k].ln-cg[k].ln));
    printf("k %d, dgam %g, tg %g\n", k, dgam[k], tg[k]);
  }

  lnum_set(Y1, 1);
  lnum_set(Y2, 1);
  lnum_set(Y3, 1);
  lnum_set(Y4, 1);
  dlnY1 = dlnY2 = dlnY3 = dlnY4 = 0;
  ddlnY1 = ddlnY2 = ddlnY3 = ddlnY4 = 0;
  for ( k = 0; k < n; k++ ) {
    lnum_imul(Y1, &cg[2*k+1]);
    lnum_imul(Y2, &sg[2*k+1]);
    lnum_imul(Y3, &cg[2*k]);
    lnum_imul(Y4, &sg[2*k]);
    dlnY1 += dgam[2*k+1] * tg[2*k+1] * m/2;
    dlnY2 += dgam[2*k+1] / tg[2*k+1] * m/2;
    dlnY3 += dgam[2*k] * tg[2*k] * m/2;
    dlnY4 += dgam[2*k] / tg[2*k] * m/2;
    ddlnY1 += ddgam[2*k+1] * tg[2*k+1] * m/2;
    ddlnY2 += ddgam[2*k+1] / tg[2*k+1] * m/2;
    ddlnY3 += ddgam[2*k] * tg[2*k] * m/2;
    ddlnY4 += ddgam[2*k] / tg[2*k] * m/2;
  }
  //printf("dlnY %g, %g, %g, %g\n", dlnY1, dlnY2, dlnY3, dlnY4);
  lnum_add(Y, lnum_add(t2, lnum_add(t1, Y1, Y2), Y3), Y4); /* Y = Y1 + Y2 + Y3 + Y4 */
  Y->ln += (sh->ln + log2)*n*m*0.5 - log2;
  r2 = lnum_get( lnum_div(t1, Y2, Y1) ); /* r2 = Y2/Y1 */
  r3 = lnum_get( lnum_div(t1, Y3, Y1) ); /* r3 = Y3/Y1 */
  r4 = lnum_get( lnum_div(t1, Y4, Y1) ); /* r4 = Y4/Y1 */
  printf("r2 %g, %g, %g; dlnY %g %g %g %g\n", r2, r3, r4, dlnY1, dlnY2, dlnY3, dlnY4);
  xp = exp(-4*beta);
  *eav = -n*m*(1+xp)/(1-xp)-(dlnY1 + dlnY2*r2 + dlnY3*r3 + dlnY4*r4)/(1 + r2 + r3 + r4);
  //*eav = n*m*(1 - 2/(1 - xp*xp));
  //*eav += -(dlnY1 + dlnY2*r2 + dlnY3*r3 + dlnY4*r4)/(1 + r2 + r3 + r4);
  free(cg);
  free(sg);
  free(tg);
  free(dgam);
  free(ddgam);
  return Y->ln;
}


int main(int argc, char **argv)
{
  int n = 4, m = 4;
  double tp = 2.3;
  lpoly_t *p;

  if ( argc == 2 ) {
    n = m = atoi(argv[1]);
  } else if ( argc >= 3 ) {
    n = atoi(argv[1]);
    m = atoi(argv[2]);
    if ( argc >= 4 ) tp = atof(argv[3]);
  }
  {
    double lnz1, eav1, cv1, lnz2 = 0, eav2 = 0, cv2 = 0;
    lnz1 = is2exact(n, m, 1./tp, &eav1, &cv1);
    lnz2 = is2_exact(n, m, 1./tp, &eav2, &cv2);
    printf("lnz %g,%g, eav %g,%g, cv %g,%g\n", lnz1, lnz2, eav1, eav2, cv1, cv2);
  }
  p = is2dos(n, m);
  if ( n*m <= 32 ) lpoly_print(p);
  is2dos_save(p, n, m);
  lpoly_close(p);
  return 0;
}

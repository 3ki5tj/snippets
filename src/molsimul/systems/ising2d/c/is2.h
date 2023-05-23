#ifndef IS2_H__
#define IS2_H__



/* two-dimensional Ising model */



#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include "mtrand.h"



typedef struct {
  int l, n;
  int M, E;
  int *s; /* 0 or 1 */
  unsigned uproba[5]; /* temporary probability for MC transitions */
  /* for Wolff's algorithm */
  int *queue;
  int *used;
} is2_t;



/* initialize an lxl Ising model */
__inline static is2_t *is2_open(int l)
{
  int i, n;
  is2_t *is;

  if ( (is = calloc(1, sizeof(*is))) == NULL ) {
    fprintf(stderr, "no memory for is2\n");
    return NULL;
  }
  is->l = l;
  is->n = n = l * l;
  if ( (is->s = calloc(n, sizeof(*is->s))) == NULL ) {
    fprintf(stderr, "no memory for is->n\n");
    return NULL;
  }
  for (i = 0; i < n; i++) is->s[i] = -1;
  is->M = -n;
  is->E = -2*n;
  is->uproba[0] = 0xffffffff;
  if ( (is->queue = calloc(n, sizeof(*is->queue))) == NULL ) {
    fprintf(stderr, "no memory for is->queue\n");
    return NULL;
  }
  if ( (is->used = calloc(n, sizeof(*is->used))) == NULL ) {
    fprintf(stderr, "no memory for is->used\n");
    return NULL;
  }
  return is;
}



__inline static void is2_close(is2_t *is)
{
  free(is->s);
  free(is->queue);
  free(is->used);
  free(is);
}



/* set transition probability */
__inline static void is2_setuproba(double bet, unsigned *p)
{
  double x = exp(-4 * bet);
  /* the index of `p` is `h` computed from `pick()`
   * the change of the energy is given by `2*h`
   * p[h] = exp(-2*h*beta); */
  p[2] = (unsigned) ((double)(0xffffffff) * x);
  p[4] = (unsigned) ((double)(0xffffffff) * x*x);
}



/* pick a random site, count neighbors with different spins */
__inline static int is2_pick(const is2_t *is, int *h)
{
  int ix, ixp, ixm, iy, iyp, iym, id, ssn;
  int l = is->l, n = is->n;

  id = (int) ( rand01() * n );
  ix = id % l;
  iy = id - ix;
  ixp = ( ix + 1 ) % l;
  ixm = ( ix + l - 1 ) % l;
  iyp = ( iy + l ) % n;
  iym = ( iy + n - l ) % n;
  ssn = is->s[iy  + ixp]
      + is->s[iy  + ixm]
      + is->s[iyp +  ix]
      + is->s[iym +  ix];
  *h = is->s[id] * ssn; /* -(*h) is the energy before the flip */
  return id;
}



/* flip site id, with (-h) is the energy before the flip */
__inline static int is2_flip(is2_t *is, int id, int h)
{
  is->M += (is->s[id] = -is->s[id]) * 2;
  return is->E += h * 2;
}


/* faster macros for systems with fixed (upon compiling) size
 * to use them one must define IS2_LB before including
 * IS2_PICK()/IS2_PSEQ() and IS2_FLIP() */
#ifdef  IS2_LB  /* L = 2^LB, N = L*L */
#define IS2_L   (1 << IS2_LB)
#define IS2_N   (IS2_L * IS2_L)

#define IS2_GETH(is, id, h) { \
  unsigned ix, ixp, ixm, iy, iyp, iym; \
  ix = id % IS2_L; \
  iy = id - ix; \
  ixp = (ix + 1) % IS2_L; \
  ixm = (ix + (IS2_L - 1)) % IS2_L; \
  iyp = (iy + IS2_L) % IS2_N; \
  iym = (iy + (IS2_N - IS2_L)) % IS2_N; \
  h = is->s[id] * ( is->s[iy  + ixp] \
                  + is->s[iy  + ixm] \
                  + is->s[iyp + ix ] \
                  + is->s[iym + ix ] ); }
#define IS2_IRND(is, id)  id = mtrand() >> (32 - 2*IS2_LB);
/* random picking */
#define IS2_PICK(is, id, h) { IS2_IRND(is, id); IS2_GETH(is, id, h); }
#define IS2_ISEQ(is, id)  id = (id + 1) % IS2_N;
/* sequential picking */
#define IS2_PSEQ(is, id, h) { IS2_ISEQ(is, id); IS2_GETH(is, id, h); }

#define IS2_FLIP(is, id, h) { \
  is->M += (is->s[id] = -is->s[id]) * 2; \
  is->E += h * 2; }

#else

#define IS2_PICK(is, id, h)  id = is2_pick(is, &h)
#define IS2_FLIP(is, id, h)  is2_flip(is, id, h)

#endif



/* compute total energy and magnetization */
__inline static int is2_em(is2_t *is)
{
  int l, n, i, j, e, m;

  e = m = 0;
  l = is->l;
  n = l * l;
  for ( i = 0; i < n; i += l ) {
    for ( j = 0; j < l; j++ ) {
      int id = i + j;
      int idr = i + (j + 1) % l;
      int idu = (i + l) % n + j;
      int s = is->s[id];
      int su = is->s[idu];
      int sr = is->s[idr];
      m += s;
      e += s * (su + sr);
    }
  }
  is->M = m;
  return is->E = -e;
}



/* add spin j to the queue if s[j] is different from s
 * return the spin */
__inline static int is2_addtoqueue(is2_t *is, int j, int s,
    double r, int *cnt)
{
  int sj = is->s[j];

  if ( sj == s && !is->used[j] && rand01() < r ) {
    is->queue[ (*cnt)++ ] = j;
    is->used[j] = (char) 1;
  }
  return sj;
}



/* Wolff algorithm */
__inline static int is2_wolff(is2_t *is, double padd)
{
  int l = is->l, n = is->n, i, ix, iy, id, s, cnt = 0, h = 0;

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  s = is->s[id];
  is->queue[ cnt++ ] = id;
  for ( i = 0; i < n; i++ ) {
    is->used[i] = 0;
  }
  is->used[id] = (char) 1;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = is->queue[i];
    /* flip the spin to correctly compute the local field,
     * which is the total magnetization of all spins
     * surrounding the cluster.
     *
     * consider a bond id-jd, with jd being a neighbor of id
     * 1) if jd does not make it to the cluster, then it
     *    lies on the border, and it contributes
     *    s[jd] to the local field
     * 2) if s[jd] == s, and will be included in the cluster
     *    in the future, it should contribute zero to the
     *    local field.  But we let it contribute s to the
     *    local field for now.  Since jd is added to the
     *    queue, when jd is considered in this loop, or
     *    when the bond jd-id is reconsidered, it will
     *    contribute an s[id] to the local field.  But at
     *    that time, s[id] = -s due to the flip here,
     *    so the total contribution would be s + (-s) = 0.
     *  */
    is->s[id] = -s;
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    h += is2_addtoqueue(is, iy + (ix + 1) % l,     s, padd, &cnt);
    h += is2_addtoqueue(is, iy + (ix + l - 1) % l, s, padd, &cnt);
    h += is2_addtoqueue(is, (iy + l) % n + ix,     s, padd, &cnt);
    h += is2_addtoqueue(is, (iy + n - l) % n + ix, s, padd, &cnt);
  }

  is->E += 2 * s * h;
  is->M -= 2 * s * cnt;
  return 0;
}



__inline static int is2_save(const is2_t *is, const char *fname)
{
  FILE *fp;
  int i, j, l, *p;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  l = is->l;
  fprintf(fp, "2 %d %d %d\n", l, l, is->n);
  for (p = is->s, i = 0; i < l; i++) {
    for (j = 0; j < l; j++, p++)
      fprintf(fp, "%c", (*p > 0) ? '#' : ' ');
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



__inline static int is2_load(is2_t *is, const char *fname)
{
  FILE *fp;
  int i, lx, ly, n, c;
  char s[80];

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fname);
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "missing first line %s\n", fname);
    fclose(fp);
    return -1;
  }
  if (4 != sscanf(s, "%d%d%d%d", &i, &lx, &ly, &n)
      || i != 2 || lx != ly || lx != is->l || n != is->n) {
    fprintf(stderr, "bad setting: %dD, %dx%d = %d\n", i, lx, ly, n);
    return -1;
  }
  for (i = 0; i < n; i++) {
    while ((c = fgetc(fp)) != EOF && c == '\n') ;
    if (c == EOF) break;
    is->s[i] = (c == ' ') ? -1 : 1;
  }
  if (i < n)
    fprintf(stderr, "%s: data stopped at i = %d\n", fname, i);
  fclose(fp);
  is2_em(is);
  return 0;
}


const double ln0 = -10000;

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



const int blksz_ = 64;

typedef struct {
  int n, ncap;
  lnum_t *a; /* coefficients */
} lpoly_t;

__inline static lpoly_t *lpoly_open(void)
{
  lpoly_t *p;
  int i;

  p = calloc(1, sizeof(*p));
  p->n = 0;
  p->ncap = blksz_;
  p->a = calloc(p->ncap, sizeof(lnum_t));
  for ( i = 0; i < p->ncap; i++ ) {
    lnum_set(p->a + i, 0);
  }
  return p;
}

__inline static void lpoly_close(lpoly_t *p)
{
  free(p->a);
  free(p);
}

__inline static void lpoly_print(lpoly_t *p)
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
__inline static void lpoly_resize(lpoly_t *p, int n)
{
  int i;

  if ( n > p->ncap ) {
    p->ncap = ((n / blksz_) + 1) * blksz_;
    p->a = realloc(p->a, p->ncap * sizeof(lnum_t));
  }
  p->n = n;
  for ( i = 0; i < p->ncap; i++ )
    lnum_set(p->a + i, 0);
}

/* q = p */
__inline static lpoly_t *lpoly_copy(lpoly_t *q, const lpoly_t *p)
{
  int i, n = p->n;

  lpoly_resize(q, n);
  for ( i = 0; i < n; i++ )
    lnum_copy(q->a + i, p->a + i);
  return q;
}

/* p = a[0] + a[1]*x + a[2]*x^2 + ... */
__inline static void lpoly_set(lpoly_t *p, int n, ...)
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
__inline static void lpoly_mullnum(lpoly_t *p, lpoly_t *q, lnum_t *ls)
{
  int i, n = q->n;

  lpoly_resize(p, n);
  for ( i = 0; i < n; i++ ) {
    p->a[i].sgn = q->a[i].sgn * ls->sgn;
    p->a[i].ln = q->a[i].ln + ls->ln;
  }
}

/* p *= s */
__inline static void lpoly_imulnum(lpoly_t *p, double s)
{
  int i;
  lnum_t ls;

  lnum_set(&ls, s);
  for ( i = 0; i < p->n; i++ ) {
    p->a[i].sgn *= ls.sgn;
    p->a[i].ln += ls.ln;
  }
}

/* p += q */
#define lpoly_iadd(p, q, t) lpoly_copy(p, lpoly_add(t, p, q))

/* r = p + q */
__inline static lpoly_t *lpoly_add(lpoly_t *r, const lpoly_t *p, const lpoly_t *q)
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
  return r;
}

#define lpoly_sub(r, p, q) lpoly_sadd(r, p, q, -1)

/* r = p + s * q */
__inline static lpoly_t *lpoly_sadd(lpoly_t *r, lpoly_t *p, lpoly_t *q, double s)
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
  return r;
}

/* y *= x */
#define lpoly_imul(y, x, t) lpoly_copy(y, lpoly_mul(t, y, x))

/* z = x * y */
__inline static lpoly_t *lpoly_mul(lpoly_t *z, const lpoly_t *x, const lpoly_t *y)
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

/* y = x^n */
__inline static lpoly_t *lpoly_pow(lpoly_t *y, const lpoly_t *x, int n, lpoly_t *t)
{
  int i;
  lpoly_copy(y, x);
  for ( i = 1; i < n; i++ )
    lpoly_imul(y, x, t);
  return y;
}


/* return the exact logarithmic partition function
 * also compute the average energy and heat capacity */
__inline static double is2_exact(int n, int m, double beta, double *eav, double *cv)
{
  lnum_t a, b, c, d, ch, sh, ccs, dg, ddg, cg, sg, z[4];
  double mh, xp, den, th, gam, dgam, ddgam, tg, ln2, ish2;
  double dlnz[4], ddlnz[4], r[4], rt = 0, dz = 0, ddz = 0;
  int k, i;

  mh = m * 0.5;
  ln2 = log(2);
  /* sh = sinh(2*K) = exp(2*beta)(1-exp(-4*beta)/2); */
  xp = exp(-2*beta);
  lnum_set(&a, (1 - xp*xp)/2);
  lnum_imul(lnum_setln(&sh, 2*beta), &a);
  lnum_copy(&c, &sh);
  c.ln = -c.ln; /* c = 1/sinh(2*K); */
  lnum_add(&ccs, &sh, &c); /* ccs = sinh(2*K) + 1/sinh(2*K) */
  ish2 = exp(2*c.ln); /* 1/sinh^2(2*K) */
  /* dg = ccs' = 2*cosh(2*K)*(1 - 1/sinh^2(2*K)) */
  lnum_set(&d, (1 + xp*xp)*(1 - ish2));
  lnum_imul(lnum_setln(&dg, 2*beta), &d);
  /* ddg = ccs'' */
  lnum_setln(&ddg, sh.ln + ln2*2 + log(1 + ish2 + 2*ish2*ish2));
  for ( i = 0; i < 4; i++ ) {
    lnum_set(&z[i], 1);
    dlnz[i] = ddlnz[i] = 0;
  }
  for ( k = 0; k < n * 2; k++ ) {
    if ( k == 0 ) {
      gam = 2*beta + log((1 - xp)/(1 + xp));
      den = 1 - xp*xp;
      dgam = 2 + 4*xp/den;
      ddgam = -8*xp*(1+xp*xp)/(den*den);
    } else {
      lnum_add(&ch, &ccs, lnum_set(&d, -cos(M_PI*k/n))); /* ch = cosh(gam) */
      /* gam = arccosh(ch) = ln(ch + sh) + ln(ch) + ln(1+th) */
      th = sqrt(1 - exp(-2*ch.ln)); /* th = tanh(gam) */
      gam = ch.ln + log(1 + th);
      /* sinh(gam) = cosh(gam) * tanh(gam); */
      lnum_mul(&sh, &ch, lnum_set(&a, th));
      dgam = lnum_get( lnum_div(&a, &dg, &sh) ); /* dg/sinh(gam) */
      lnum_imul(lnum_set(&c, dgam*dgam), &ch);
      lnum_sub(&d, &ddg, &c); /* d = ddg - cosh(gam)*dgam^2 */
      ddgam = lnum_get( lnum_div(&b, &d, &sh) ); /* dg/sinh(gam) */
    }
    lnum_setln(&a,  m * gam * 0.5); /* a = exp(m*gam/2) */
    lnum_setln(&b, -m * gam * 0.5); /* b = exp(-m*gam/2) */
    lnum_add(&cg, &a, &b); /* cg = 2 cosh(m*gam/2) */
    lnum_sub(&sg, &a, &b); /* sg = 2 sinh(m*gam/2) */
    tg = lnum_get( lnum_div(&c, &sg, &cg) );
    i = (k % 2) * 2;
    lnum_imul(&z[i],   &cg);
    lnum_imul(&z[i+1], &sg);
    dlnz[i]   += dgam * tg * mh;
    dlnz[i+1] += dgam / tg * mh;
    ddlnz[i]   += ( dgam*dgam*4*exp(-2*cg.ln)*mh + ddgam*tg) * mh;
    ddlnz[i+1] += (-dgam*dgam*4*exp(-2*sg.ln)*mh + ddgam/tg) * mh;
  }

  lnum_add(&d, lnum_add(&b, lnum_add(&a, &z[0], &z[1]), &z[2]), &z[3]);
  lnum_setln(&a, 2*beta);
  lnum_sub(&sh, &a, lnum_setln(&b, -2*beta));
  for ( i = 0; i < 4; i++ ) {
    r[i] = lnum_get( lnum_div(&a, &z[i], &z[0]) );
    rt += r[i];
    dz += r[i] * dlnz[i];
    ddz += r[i] * (ddlnz[i] + dlnz[i] * dlnz[i]);
  }
  dz /= rt;
  *eav = n*m*(1 - 2/(1 - xp*xp)) - dz;
  *cv = beta * beta * (ddz/rt - dz * dz - 2*n*m*ish2);
  return d.ln + sh.ln*n*m*0.5 - ln2;
}

/* exact density of states */
__inline static double *is2dos(int n, int m)
{
  lpoly_t *a, *b, *c, *s, *bm, *u, *v, *w, *xp, *xm;
  lpoly_t **p, **q, *Z1, *Z2, *Z3, *Z4;
  lnum_t tmp;
  double *lndos, lnmfac, lnimax;
  int i, k;

  a = lpoly_open();
  b = lpoly_open();
  c = lpoly_open();
  s = lpoly_open();
  u = lpoly_open();
  v = lpoly_open();
  w = lpoly_open();
  bm = lpoly_open();
  xp = lpoly_open();
  xm = lpoly_open();
  Z1 = lpoly_open();
  Z2 = lpoly_open();
  Z3 = lpoly_open();
  Z4 = lpoly_open();
  p = calloc(m + 1, sizeof(*p));
  q = calloc(m + 1, sizeof(*q));
  for ( i = 0; i <= m; i++ ) {
    p[i] = lpoly_open();
    q[i] = lpoly_open();
  }
  lpoly_set(b, 4, 0.0, 2.0, 0.0, -2.0); /* beta = 2x(1-x^2) */
  lpoly_pow(u, b, m, v); /* beta^m */
  lnum_setln(&tmp, (m - 1)*log(0.5)); /* tmp = 0.5^(m-1) */
  lpoly_mullnum(bm, u, &tmp); /* bm = 0.5^(m-1) beta^m */

  for ( lnmfac = 0, i = 2; i <= m; i++ )
    lnmfac += log(i * 0.5); /* lnmfac = ln(m!/2^(m-1)) */

  lpoly_set(u, 2, 1.0, 1.0); /* u = 1 + x */
  lpoly_pow(xp, u, m, v); /* xp = (1 + x)^m */
  lpoly_set(u, 2, 1.0, -1.0); /* u = 1 - x */
  lpoly_pow(xm, u, m, v); /* xm = (1 - x)^m */
  lpoly_resize(u, m + 1);
  lnum_set(u->a + m, 1); /* u = x^m */
  lpoly_mul(v, u, xp); /* v = x^m (1 + x)^m */
  lpoly_mul(w, u, xm); /* w = x^m (1 - x)^m */
  lpoly_add(Z3, xm, v); /* Z3 = c0 = (1 - x)^m + x^m (1 + x)^m */
  lpoly_sub(Z4, xm, v); /* Z4 = s0 = (1 - x)^m - x^m (1 + x)^m */
  lpoly_add(Z1, xp, w); /* Z1 = cn = (1 + x)^m + x^m (1 - x)^m */
  lpoly_sub(Z2, xp, w); /* Z2 = sn = (1 + x)^m - x^m (1 - x)^m */

  if ( n % 2 == 0 ) { /* n is even */
    lpoly_imul(Z3, Z1, u);
    lpoly_imul(Z4, Z2, u);
    lpoly_set(Z1, 1, 1.0);
    lpoly_set(Z2, 1, 1.0);
  }

  for ( k = 1; k < n; k++ ) {
    lnum_set(&tmp, -cos(M_PI*k/n));
    lpoly_mullnum(u, b, &tmp); /* u = beta*cos(pi*k/n) */
    lpoly_set(v, 5, 1.0, 0.0, 2.0, 0.0, 1.0); /* v = (1+x^2)^2 */
    lpoly_add(a, v, u); /* ak = (1+x^2)^2 - beta*cos(pi*k/n) */

    /* p[i] = (ak^2 - beta^2)^i/(2i)! */
    lpoly_pow(u, a, 2, w);
    lpoly_pow(v, b, 2, w);
    lpoly_sub(w, u, v); /* w = ak^2 - beta^2 */
    lpoly_set(p[0], 1, 1.0);
    for ( i = 1; i <= m/2; i++ ) {
      lpoly_mul(p[i], p[i-1], w);
      lpoly_imulnum(p[i], 1.0/(2.*i*(2*i-1)));
    }
    /* q[i] = ak^i m!/i!/2^(m-1) */
    lpoly_set(q[0], 1, 1.0);
    lnum_setln(&q[0]->a[0], lnmfac);
    for ( i = 1; i <= m; i++ ) {
      lpoly_mul(q[i], q[i-1], a);
      lpoly_imulnum(q[i], 1./i);
    }

    /* ck^2 = bm + Sum_{i=0 to m/2} p[i]*q[m-2i] */
    lpoly_copy(c, bm);
    for ( i = 0; i <= m/2; i++ ) {
      lpoly_mul(u, p[i], q[m-2*i]);
      lpoly_iadd(c, u, v); /* c += u */
    }
    /* sk^2 = ck^2 - 2*(beta^m/2^(m-1)) */
    lpoly_sadd(s, c, bm, -2.0);
    if ( k % 2 ) {
      lpoly_imul(Z1, c, u);
      lpoly_imul(Z2, s, u);
    } else {
      lpoly_imul(Z3, c, u);
      lpoly_imul(Z4, s, u);
    }
  }
  lpoly_iadd(Z1, Z2, u);
  lpoly_iadd(Z1, Z3, u);
  lpoly_iadd(Z1, Z4, u);
  lpoly_imulnum(Z1, 0.5);

  /* export to a 1D array */
  lnimax = log((double) ULONG_MAX);
  lndos = calloc(n*m + 1, sizeof(*lndos));
  for ( i = 0; i < Z1->n; i+= 2 ) {
    double y = Z1->a[i].ln;
    if ( y < 0 ) y = ln0;
    else if ( y < lnimax ) /* round to nearest integer */
      y = log((double)((unsigned long) (exp(y)+.5)));
    lndos[i/2] = y;
  }

  lpoly_close(a);
  lpoly_close(b);
  lpoly_close(c);
  lpoly_close(s);
  lpoly_close(u);
  lpoly_close(v);
  lpoly_close(w);
  lpoly_close(bm);
  lpoly_close(xp);
  lpoly_close(xm);
  lpoly_close(Z1);
  lpoly_close(Z2);
  lpoly_close(Z3);
  lpoly_close(Z4);
  for ( i = 0; i <= m; i++ ) {
    lpoly_close(p[i]);
    lpoly_close(q[i]);
  }
  free(p);
  free(q);
  return lndos;
}

__inline static int is2dos_save(const double *lndos, int n, int m,
    const char *prefix)
{
  char fn[128];
  FILE *fp;
  int i, digs;

  sprintf(fn, "%sis2lndos%dx%d.dat", prefix, n, m);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  digs = (int)(-log10(lndos[n*m/2]*DBL_EPSILON));
  if ( digs < 0 ) digs = 0;
  for ( i = 0; i <= n*m; i++ ) {
    fprintf(fp, "%d %.*f\n", -2*n*m+4*i, digs, lndos[i]);
  }
  fclose(fp);
  return 0;
}

#endif /* IS2_H__ */


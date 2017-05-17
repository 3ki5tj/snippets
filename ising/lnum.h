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

/* r = p + q */
__inline static void lpoly_add(lpoly_t *r, lpoly_t *p, lpoly_t *q)
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
__inline static void lpoly_sadd(lpoly_t *r, lpoly_t *p, lpoly_t *q, double s)
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

/* y *= x */
__inline static lpoly_t *lpoly_imul(lpoly_t *y, const lpoly_t *x, lpoly_t *t)
{
  lpoly_mul(t, y, x);
  lpoly_copy(y, t);
  return y;
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



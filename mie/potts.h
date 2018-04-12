#include "mtrand.h"
#include "util.h"

enum {
  SDRT, /* direct */
  SLIN, /* linear extrapolation */
  SEXP, /* exponential extrapolation */
  SBLK,
  SBAV,
  STOT
};

/* one-dimensional Potts model chain (open boundary) */
typedef struct {
  int q;
  int n;
  int *s;
  int e;

  long tot; /* total number of samples */
  long **cnt1d; /* raw counts for 1D entropy */
  long **cnt2d; /* raw counts for 2D entropy */
  long **cnt3d; /* raw counts for 3D entropy */
  double *ent1d[STOT];
  double *ent2d[STOT];
  double *ent3d[STOT];
  double ent1[STOT], ent1r;
  double ent2[STOT], ent2r;
  double ent3[STOT], ent3r;
  double entr;

  int trjn, trji; /* trajectory length and current position */
  int trjwsz; /* size to save each frame */
  char *trj;
} potts_t;

static potts_t *potts_open(int q, int n, int trjlen)
{
  potts_t *p;
  int i, k, s;

  xnew(p, 1);
  p->q = q;
  p->n = n;
  xnew(p->s, n);
  for ( i = 0; i < n; i++ ) {
    p->s[i] = 0;
  }
  p->e = -(n - 1);

  p->tot = 0;
  xnew(p->cnt1d, n);
  for ( i = 0; i < n; i++ ) {
    xnew(p->cnt1d[i], q);
    for ( s = 0; s < q; s++ )
      p->cnt1d[i][s] = 0;
  }

  xnew(p->cnt2d, n*n);
  for ( i = 0; i < n*n; i++ ) {
    xnew(p->cnt2d[i], q*q);
    for ( s = 0; s < q*q; s++ )
      p->cnt2d[i][s] = 0;
  }

  xnew(p->cnt3d, n*n*n);
  for ( i = 0; i < n*n*n; i++ ) {
    xnew(p->cnt3d[i], q*q*q);
    for ( s = 0; s < q*q*q; s++ )
      p->cnt3d[i][s] = 0;
  }

  for ( k = 0; k < STOT; k++ ) {
    xnew(p->ent1d[k], n);
    for ( i = 0; i < n; i++ )
      p->ent1d[k][i] = 0;

    xnew(p->ent2d[k], n*n);
    for ( i = 0; i < n*n; i++ )
      p->ent2d[k][i] = 0;

    xnew(p->ent3d[k], n*n*n);
    for ( i = 0; i < n*n*n; i++ )
      p->ent3d[k][i] = 0;
  }

  /* initialize the trajectory data */
  p->trjn = trjlen;
  p->trji = 0;
  p->trj = NULL;
  p->trjwsz = n;
  if ( q >= 128 ) {
    fprintf(stderr, "traj. mode does not support %d state\n", q);
    return NULL;
  }
  //fprintf(stderr, "using trajectory of size %gM\n", trjlen*n/(1024.*1024));
  xnew(p->trj, p->trjn * p->trjwsz);
  return p;
}

static void potts_close(potts_t *p)
{
  int i, k, n = p->n;

  for ( i = 0; i < n; i++ ) free(p->cnt1d[i]);
  free(p->cnt1d);
  for ( i = 0; i < n*n; i++ ) free(p->cnt2d[i]);
  free(p->cnt2d);
  for ( i = 0; i < n*n*n; i++ ) free(p->cnt3d[i]);
  free(p->cnt3d);
  for ( k = 0; k < STOT; k++ ) {
    free(p->ent1d[k]);
    free(p->ent2d[k]);
    free(p->ent3d[k]);
  }
  free(p->s);
  free(p->trj);
  free(p);
}

static void potts_encode(potts_t *p, int trji, const int *s)
{
  int i;
  char *c = p->trj + trji * p->trjwsz;

  if ( trji >= p->trjn ) {
    fprintf(stderr, "encode: index overflow %d >= %d\n", trji, p->trjn);
    exit(1);
  }
  for ( i = 0; i < p->n; i++ )
    c[i] = (char) s[i];
}

static void potts_decode(potts_t *p, int trji, int *s)
{
  int i;
  char *c = p->trj + trji * p->trjwsz;

  if ( trji >= p->trjn ) {
    fprintf(stderr, "decode: index overflow %d >= %d\n", trji, p->trjn);
    exit(1);
  }
  for ( i = 0; i < p->n; i++ )
    s[i] = c[i];
}

/* reset counts */
static void potts_resetcounts(potts_t *p)
{
  int i, s, n = p->n, q = p->q;

  p->tot = 0;

  for ( i = 0; i < n; i++ )
    for ( s = 0; s < q; s++ )
      p->cnt1d[i][s] = 0;

  for ( i = 0; i < n*n; i++ )
    for ( s = 0; s < q*q; s++ )
      p->cnt2d[i][s] = 0;

  for ( i = 0; i < n*n*n; i++ )
    for ( s = 0; s < q*q*q; s++ )
      p->cnt3d[i][s] = 0;
}


/* register the state ps into the counts */
static void potts_reg(potts_t *p, const int *ps)
{
  int i, j, k, n = p->n, s, t, r, q = p->q;

  p->tot += 1;

  for ( i = 0; i < n; i++ ) {
    s = ps[i];
    p->cnt1d[i][s] += 1;
  }

  for ( i = 0; i < n; i++ ) {
    s = ps[i];
    for ( j = i + 1; j < n; j++ ) {
      t = ps[j];
      p->cnt2d[i*n + j][s*q + t] += 1;
    }
  }

  for ( i = 0; i < n; i++ ) {
    s = ps[i];
    for ( j = i + 1; j < n; j++ ) {
      t = ps[j];
      for ( k = j + 1; k < n; k++ ) {
        r = ps[k];
        p->cnt3d[(i*n + j)*n + k][(s*q + t)*q + r] += 1;
      }
    }
  }
}

/* add a trajectory frame */
static void potts_add(potts_t *p)
{
  potts_encode(p, p->trji, p->s);
  p->trji++;
}

/* count occurrences from frame start to frame end */
static void potts_count(potts_t *p, int start, int end)
{
  int i, *s;

  xnew(s, p->n);
  potts_resetcounts(p);
  for ( i = start; i < end; i++ ) {
    potts_decode(p, i, s);
    potts_reg(p, s);
  }
  free(s);
}

/* compute the first-order approximation of the entropy */
static double potts_ent1(potts_t *p, double *ent1d)
{
  double ent1 = 0, ent1i, tot = (double) p->tot, pr;
  int i, n = p->n, q = p->q, s;
  long c;

  for ( i = 0; i < n; i++ ) {
    ent1i = 0;
    for ( s = 0; s < q; s++ ) {
      c = p->cnt1d[i][s];
      if ( c <= 0 ) continue;
      pr = 1.0*c/tot;
      ent1i += -pr*log(pr);
    }
    ent1d[i] = ent1i;
    ent1 += ent1i;
  }
  return ent1;
}

/* compute the second-order approximation of the entropy */
static double potts_ent2(potts_t *p,
    const double *ent1d, double *ent2d)
{
  double ent2 = 0, eij, tot = (double) p->tot, pr, ds;
  int i, j, ij, n = p->n, s, t, st, q = p->q;

  for ( i = 0; i < n; i++ )
    ent2 += ent1d[i];

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      ij = i*n + j;
      eij = 0;
      for ( s = 0; s < q; s++ ) {
        for ( t = 0; t < q; t++ ) {
          st = s * q + t;
          long c = p->cnt2d[ij][st];
          if ( c <= 0 ) continue;
          pr = 1.0*c/tot;
          eij += -pr*log(pr);
        }
      }
      ent2d[ij] = eij;
      ds = eij - ent1d[i] - ent1d[j];
      /* the pair correction has to be nonpositive */
      if ( ds > 0 ) ds = 0;
      ent2 += ds;
    }
  }
  return ent2;
}

/* compute the third-order approximation of the entropy */
static double potts_ent3(potts_t *p,
    const double *ent1d, const double *ent2d, double *ent3d)
{
  double ent3 = 0, eijk, tot = (double) p->tot, pr, ds;
  int i, j, k, ij, ijk, n = p->n, s, t, r, str, q = p->q;

  for ( i = 0; i < n; i++ )
    ent3 += ent1d[i];

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      ij = i * n + j;
      ent3 += ent2d[ij] - ent1d[i] - ent1d[j];
      for ( k = j + 1; k < n; k++ ) {
        ijk = ij * n + k;
        eijk = 0;
        for ( s = 0; s < q; s++ ) {
          for ( t = 0; t < q; t++ ) {
            for ( r = 0; r < q; r++ ) {
              str = (s * q + t) * q + r;
              long c = p->cnt3d[ijk][str];
              if ( c <= 0 ) continue;
              pr = 1.0*c/tot;
              eijk += -pr*log(pr);
            }
          }
        }
        ent3d[ijk] = eijk;
        ds = eijk - ent2d[ij] - ent2d[i*n+k] - ent2d[j*n+k]
           + ent1d[i] + ent1d[j] + ent1d[k];
        ent3 += ds;
      }
    }
  }
  //printf("%g, %g\n", ent3d[n+2], ent3); getchar();
  return ent3;
}



/* linear extrapolation */
static double linext(double xt, double xb, double t, double tb)
{
  /* since we have St = S - a/t, Sp = S - a*npart/t
   * a/t = (St - Sp)/(npart-1)
   * S = (npart*St - Sp)/(npart-1); */
  return (xt*t - xb*tb)/(t - tb);
}


/* exponential extrapolation */
static double expext(double xt, double xb, double t, double tb)
{
  double ds, xp, xpmax = 0.99 * t / tb, xpc;
  ds = xt - xb;
  xp = exp(ds);
  if ( xp >= xpmax ) xp = xpmax;
  xpc = (t - xp*tb)/(t - tb);
  return xt - log(xpc);
}


static double potts_entropy(potts_t *p, int npart)
{
  int trjn = p->trji, ip, blksz;
  int i, j, k, ij, ijk, n = p->n;
  double ds;

  /* entropy estimated from trajectory block(s) */
  blksz = trjn / npart;
  for ( i = 0; i < n; i++ )
    p->ent1d[SBAV][i] = 0;
  for ( i = 0; i < n*n; i++ )
    p->ent2d[SBAV][i] = 0;
  for ( i = 0; i < n*n*n; i++ )
    p->ent3d[SBAV][i] = 0;
  p->ent1[SBAV] = 0;
  p->ent2[SBAV] = 0;
  p->ent3[SBAV] = 0;
  /* average over the blocks */
  for ( ip = 0; ip < npart; ip++ ) {
    /* compute the entropy from the ith block */
    potts_count(p, ip * blksz, (ip + 1) * blksz);
    p->ent1[SBLK] = potts_ent1(p, p->ent1d[SBLK]);
    p->ent2[SBLK] = potts_ent2(p, p->ent1d[SBLK], p->ent2d[SBLK]);
    p->ent3[SBLK] = potts_ent3(p, p->ent1d[SBLK], p->ent2d[SBLK], p->ent3d[SBLK]);
    for ( i = 0; i < n; i++ ) {
      p->ent1d[SBAV][i] += p->ent1d[SBLK][i];
      for ( j = i + 1; j < n; j++ ) {
        ij = i * n + j;
        p->ent2d[SBAV][ij] += p->ent2d[SBLK][ij];
        for ( k = j + 1; k < n; k++ ) {
          ijk = ij * n + k;
          p->ent3d[SBAV][ijk] += p->ent3d[SBLK][ijk];
        }
      }
    }
    p->ent1[SBAV] += p->ent1[SBLK];
    p->ent2[SBAV] += p->ent2[SBLK];
    p->ent3[SBAV] += p->ent3[SBLK];
  }
  /* compute the block averages from the block sums */
  for ( i = 0; i < n; i++ ) {
    p->ent1d[SBAV][i] /= npart;
    for ( j = i + 1; j < n; j++ ) {
      ij = i * n + j;
      p->ent2d[SBAV][ij] /= npart;
      for ( k = j + 1; k < n; k++ ) {
        ijk = ij * n + k;
        p->ent3d[SBAV][ijk] /= npart;
      }
    }
  }
  p->ent1[SBAV] /= npart;
  p->ent2[SBAV] /= npart;
  p->ent3[SBAV] /= npart;

  /* entropy estimate from the entire trajectory */
  potts_count(p, 0, trjn);
  p->ent1[SDRT] = potts_ent1(p, p->ent1d[SDRT]);
  p->ent2[SDRT] = potts_ent2(p, p->ent1d[SDRT], p->ent2d[SDRT]);
  p->ent3[SDRT] = potts_ent3(p, p->ent1d[SDRT], p->ent2d[SDRT], p->ent3d[SDRT]);

  /* exponential extrapolation */
  p->ent1[SLIN] = 0;
  p->ent2[SLIN] = 0;
  p->ent3[SLIN] = 0;
  p->ent1[SEXP] = 0;
  p->ent2[SEXP] = 0;
  p->ent3[SEXP] = 0;
  for ( i = 0; i < n; i++ ) {
    if ( p->ent1d[SBAV][i] > p->ent1d[SDRT][i] ) {
      fprintf(stderr, "block average %d, is greater %g > %g\n", i, p->ent1d[SBAV][i], p->ent1d[SDRT][i]);
      exit(1);
    }
    p->ent1d[SLIN][i] = linext(p->ent1d[SDRT][i], p->ent1d[SBAV][i], trjn, blksz);
    p->ent1[SLIN] += p->ent1d[SLIN][i];

    p->ent1d[SEXP][i] = expext(p->ent1d[SDRT][i], p->ent1d[SBAV][i], trjn, blksz);
    p->ent1[SEXP] += p->ent1d[SEXP][i];
    for ( j = i + 1; j < n; j++ ) {
      ij = i * n + j;
      if ( p->ent2d[SBAV][ij] > p->ent2d[SDRT][ij] ) {
        fprintf(stderr, "block average %d %d, is greater %g > %g\n", i, j, p->ent2d[SBAV][ij], p->ent2d[SDRT][ij]);
        exit(1);
      }
      p->ent2d[SLIN][ij] = linext(p->ent2d[SDRT][ij], p->ent2d[SBAV][ij], trjn, blksz);
      ds = p->ent2d[SLIN][ij] - p->ent1d[SLIN][i] - p->ent1d[SLIN][j];
      if ( ds > 0 ) ds = 0;
      p->ent2[SLIN] += ds;

      p->ent2d[SEXP][ij] = expext(p->ent2d[SDRT][ij], p->ent2d[SBAV][ij], trjn, blksz);
      ds = p->ent2d[SEXP][ij] - p->ent1d[SEXP][i] - p->ent1d[SEXP][j];
      if ( ds > 0 ) ds = 0;
      p->ent2[SEXP] += ds;

      for ( k = j + 1; k < n; k++ ) {
        ijk = ij * n + k;
        p->ent3d[SLIN][ijk] = linext(p->ent3d[SDRT][ijk], p->ent3d[SBAV][ijk], trjn, blksz);
        ds = p->ent3d[SLIN][ijk] - p->ent2d[SLIN][ij]
           - p->ent2d[SLIN][i*n + k] - p->ent2d[SLIN][j*n + k]
           + p->ent1d[SLIN][i] + p->ent1d[SLIN][j] + p->ent1d[SLIN][k];
        p->ent3[SLIN] += ds;

        p->ent3d[SEXP][ijk] = expext(p->ent3d[SDRT][ijk], p->ent3d[SBAV][ijk], trjn, blksz);
        ds = p->ent3d[SEXP][ijk] - p->ent2d[SEXP][ij]
           - p->ent2d[SEXP][i*n + k] - p->ent2d[SEXP][j*n + k]
           + p->ent1d[SEXP][i] + p->ent1d[SEXP][j] + p->ent1d[SEXP][k];
        p->ent3[SEXP] += ds;
      }
    }
  }
  p->ent2[SLIN] += p->ent1[SLIN];
  p->ent2[SEXP] += p->ent1[SEXP];
  p->ent3[SLIN] += p->ent2[SLIN];
  p->ent3[SEXP] += p->ent2[SEXP];

  //fprintf(stderr, "%g, %g\n", p->ent1[SLIN], p->ent2[SLIN]); getchar();
  //p->ent1[SEXP] = expext(p->ent1[SDRT], p->ent1[SBAV], trjn, blksz);
  //p->ent2[SEXP] = expext(p->ent2[SDRT], p->ent2[SBAV], trjn, blksz);

  return p->ent3[SEXP];
}

__inline static int potts_energy(potts_t *p)
{
  int i, e = 0;
  for ( i = 0; i < p->n - 1; i++ ) {
    e -= (p->s[i] == p->s[i+1]);
  }
  return e;
}

static int potts_metro(potts_t *p, double beta)
{
  int i, s, s1, q = p->q, de, acc;

  i = (int) (p->n * rand01());
  s = p->s[i];
  s1 = (s + 1 + (int)((q - 1) * rand01())) % q;
  de = 0;
  if ( i > 0 ) {
    if ( p->s[i-1] == s ) de += 1;
    else if ( p->s[i-1] == s1 ) de -= 1;
  }
  if ( i < p->n - 1 ) {
    if ( p->s[i+1] == s ) de += 1;
    else if ( p->s[i+1] == s1 ) de -= 1;
  }
  //printf("i %d, de %d, s %d, s1 %d\n", i, de, s, s1);
  acc = 0;
  if ( de <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp(-beta*de) );
  }
  if ( acc ) {
    p->s[i] = s1;
    p->e += de;
  }
  return acc;
}

/* compute the reference entropy, 1st order */
__inline static double potts_ent1ref(potts_t *p)
{
  return p->n * log(p->q);
}

const double LN0 = -1000000;

/* return log(exp(a) + exp(b)) */
__inline static double lnadd(double a, double b)
{
  if ( a > b ) {
    return a + log(1 + exp(b-a));
  } else {
    return b + log(1 + exp(a-b));
  }
}

/* A = B * C */
__inline static void lnmatmul(double *a, double *b, double *c, int q)
{
  int s, t, u;
  double x;

  for ( s = 0; s < q; s++ ) {
    for ( t = 0; t < q; t++ ) {
      x = LN0;
      for ( u = 0; u < q; u++ ) {
        /* exp(x) += exp(b(s,u)) * exp(c(u,t)) */
        x = lnadd(x, b[s*q + u] + c[u*q + t]);
      }
      a[s*q + t] = x;
    }
  }
}

/* compute the reference entropy */
__inline static double potts_entref(potts_t *p, double beta)
{
  int q = p->q, k, n = p->n;
  double f = exp(beta) - 1, lna, lnb, lnz, lnp1, lnp2, sk, dsk;

  p->ent1r = n * log(q);

  /* compute the correction from pair information */
  p->ent2r = 0;
  for ( k = 1; k < n; k++ ) {
    /* for the transition matrix of a pair k sites apart */
    lna = k * log(f);
    /* ai = f^k, bi = ((f+q)^k - f^k)/q */
    lnb = k * log(f + q) - log(q) + log(1 - pow(f/(f+q), k));
    lnz = log(q) + k * log(f + q); /* partition function */
    lnp1 = lnadd(lna, lnb) - lnz; /* diagonal element of the transition matrix */
    lnp2 = lnb - lnz; /* non-diagonal element of the transition matrix */
    /* entropy of the pair of k sites apart
       sk = -(p2*ln(p2)*(q^2 - q) + p1*ln(p1)*q) */
    sk = -(lnp2*exp(lnp2)*(q*q - q) + lnp1*exp(lnp1)*q);
    dsk = sk - 2 * log(q);
    p->ent2r += dsk * (n - k);
    //fprintf(stderr, "k %d: e^s %g\n", k, exp(sk));
    p->ent3r += -dsk * (n - k) * (k - 1);
  }
  //fprintf(stderr, "ent2r %g\n", p->ent2r);
  p->ent2r += p->ent1r;
  p->ent3r += p->ent2r;
  p->entr = log(q) + (n - 1) * (log(f + q) - beta * (f + 1) / (f + q));

  return p->ent2r;
}

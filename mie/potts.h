#include "mtrand.h"
#include "util.h"

/* one-dimensional Potts model chain (open boundary) */
typedef struct {
  int q;
  int n;
  int *s;
  int e;

  long tot; /* total number of samples */
  long **cnt1; /* raw counts for 1D entropy */
  long **cnt2; /* raw counts for 2D entropy */
  double *enti;
  double ent1, ent1c, ent1r;
  double ent2, ent2c, ent2r;
  double entr;

  int trjn, trji; /* trajectory length and current position */
  int trjwsz; /* size to save each frame */
  char *trj;
} potts_t;

static potts_t *potts_open(int q, int n, int trjlen)
{
  potts_t *p;
  int i, s;

  xnew(p, 1);
  p->q = q;
  p->n = n;
  xnew(p->s, n);
  for ( i = 0; i < n; i++ ) {
    p->s[i] = 0;
  }
  p->e = -(n - 1);

  p->tot = 0;
  xnew(p->cnt1, n);
  for ( i = 0; i < n; i++ ) {
    xnew(p->cnt1[i], q);
    for ( s = 0; s < q; s++ )
      p->cnt1[i][s] = 0;
  }
  xnew(p->enti, n);

  xnew(p->cnt2, n*n);
  for ( i = 0; i < n*n; i++ ) {
    xnew(p->cnt2[i], q*q);
    for ( s = 0; s < q*q; s++ )
      p->cnt2[i][s] = 0;
  }

  /* initialize the trajectory data */
  p->trjn = trjlen;
  p->trji = 0;
  p->trj = NULL;
  p->trjwsz = n;
  if ( trjlen > 0 ) {
    if ( q >= 128 ) {
      fprintf(stderr, "traj. mode does not support %d state\n", q);
      return NULL;
    }
    //fprintf(stderr, "using trajectory of size %gM\n", trjlen*n/(1024.*1024));
    xnew(p->trj, p->trjn * p->trjwsz);
  }
  return p;
}

static void potts_close(potts_t *p)
{
  int i, n = p->n;

  for ( i = 0; i < n; i++ ) {
    free(p->cnt1[i]);
  }
  free(p->cnt1);
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

/* empty counts */
static void potts_empty(potts_t *p)
{
  int i, s, n = p->n, q = p->q;

  p->tot = 0;

  for ( i = 0; i < n; i++ )
    for ( s = 0; s < q; s++ )
      p->cnt1[i][s] = 0;

  for ( i = 0; i < n*n; i++ )
    for ( s = 0; s < q*q; s++ )
      p->cnt2[i][s] = 0;
}


/* register the state ps into the counts */
static void potts_reg(potts_t *p, const int *ps)
{
  int i, j, n = p->n, s, t, q = p->q;

  p->tot += 1;

  for ( i = 0; i < n; i++ ) {
    s = ps[i];
    p->cnt1[i][s] += 1;
  }

  for ( i = 0; i < n; i++ ) {
    s = ps[i];
    for ( j = i + 1; j < n; j++ ) {
      t = ps[j];
      p->cnt2[i*n + j][s*q + t] += 1;
    }
  }
}

/* add a trajectory frame */
static void potts_add(potts_t *p)
{
  if ( p->trjn > 0 ) {
    potts_encode(p, p->trji, p->s);
    p->trji++;
  } else { /* register directly */
    potts_reg(p, p->s);
  }
}

/* count occurrences from frame start to frame end */
static void potts_count(potts_t *p, int start, int end)
{
  int i, *s;

  xnew(s, p->n);
  potts_empty(p);
  for ( i = start; i < end; i++ ) {
    potts_decode(p, i, s);
    potts_reg(p, s);
  }
  free(s);
}

/* compute the first-order approximation of the entropy */
static double potts_ent1(potts_t *p)
{
  double ent1 = 0, tot = (double) p->tot, pr;
  int i, n = p->n, s;

  for ( i = 0; i < n; i++ ) {
    p->enti[i] = 0;
    for ( s = 0; s < p->q; s++ ) {
      long c = p->cnt1[i][s];
      if ( c <= 0 ) continue;
      pr = 1.0*c/tot;
      p->enti[i] += -pr*log(pr);
    }
    ent1 += p->enti[i];
  }
  p->ent1 = ent1;
  return ent1;
}

/* compute the second-order approximation of the entropy */
static double potts_ent2(potts_t *p)
{
  double ent2 = 0, eij, tot = (double) p->tot, pr;
  int i, j, ij, n = p->n, s, t, st, q = p->q;

  ent2 = potts_ent1(p);
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      ij = i*n + j;
      eij = 0;
      for ( s = 0; s < q; s++ ) {
        for ( t = 0; t < q; t++ ) {
          st = s * q + t;
          long c = p->cnt2[ij][st];
          if ( c <= 0 ) continue;
          pr = 1.0*c/tot;
          eij += -pr*log(pr);
        }
      }
      ent2 += eij - p->enti[i] - p->enti[j];
    }
  }
  p->ent2 = ent2;
  return ent2;
}


static double potts_entropy(potts_t *p, int npart, int blk1st)
{
  if ( p->trjn == 0 ) {
    potts_ent2(p);
    p->ent1c = p->ent1;
    p->ent2c = p->ent2;
  } else {
    int trjn = p->trji, ip, blksz;
    double ent1p = 0, ent1t;
    double ent2p = 0, ent2t;

    /* entropy estimated from trajectory block(s) */
    blksz = trjn / npart;
    if ( blk1st ) { /* just from the first block */
      for ( ip = 0; ip < npart; ip++ ) {
        potts_count(p, ip * blksz, (ip + 1) * blksz);
        potts_ent2(p);
        ent1p += p->ent1;
        ent2p += p->ent2;
      }
      ent1p /= npart;
      ent2p /= npart;
    } else { /* average over several blocks */
      potts_count(p, 0, blksz);
      ent2p = potts_ent2(p);
      ent1p = p->ent1;
    }

    /* entropy estimate from the entire trajectory */
    potts_count(p, 0, trjn);
    potts_ent2(p);
    ent1t = p->ent1;
    ent2t = p->ent2;

    /* since we have St = S - a/t, Sp = S - a*npart/t
     * a/t = (St - Sp)/(npart-1)
     * S = (npart*St - Sp)/(npart-1); */
    p->ent1 = ent1t;
    p->ent1c = (ent1t*npart - ent1p)/(npart - 1);
    p->ent2 = ent2t;
    p->ent2c = (ent2t*npart - ent2p)/(npart - 1);
  }
  return p->ent2;
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
__inline static double potts_ent2ref(potts_t *p, double beta)
{
  int q = p->q, k, n = p->n;
  double f = exp(beta) - 1, lna, lnb, lnz, lnp1, lnp2, sk;

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
    p->ent2r += (sk - 2*log(q)) * (n - k);
    //fprintf(stderr, "k %d: e^s %g\n", k, exp(sk));
  }
  //fprintf(stderr, "ent2r %g\n", p->ent2r);
  p->ent1r = n * log(q);
  p->ent2r += p->ent1r;
  p->entr = log(q) + (n - 1) * (log(f + q) - beta * (f + 1) / (f + q));

  return p->ent2r;
}

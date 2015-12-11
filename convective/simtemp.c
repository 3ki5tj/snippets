/* simulated tempering */
#include "tconv.h"

#ifndef IS2_LB
#define IS2_LB 5
#endif
#define L (1 << IS2_LB)
#define N ((L) * (L))
#include "is2.h"

#define BATCHSIZE N


typedef struct {
  int type;

  int nbeta;
  double *beta;
  double *hist;
  double *esum;
  double *eesum;

  /* Ising model parameters */
  double *lnz;
  unsigned (*uproba)[5];
  double *eref;
  double *cvref;

  double *p;
  int *idmap;
  double *ps;
  double *cp;
} simtemp_t;


static simtemp_t *simtemp_open(int type,
    double tp0, double dtp, int ntp)
{
  simtemp_t *st;
  int i;
  double tp;

  xnew(st, 1);

  st->type = type;

  st->nbeta = ntp;

  xnew(st->beta, st->nbeta);
  xnew(st->hist, st->nbeta);
  xnew(st->esum, st->nbeta);
  xnew(st->eesum, st->nbeta);

  xnew(st->lnz, st->nbeta);
  xnew(st->uproba, st->nbeta);
  xnew(st->eref, st->nbeta);
  xnew(st->cvref, st->nbeta);
  for ( i = 0; i < st->nbeta; i++ ) {
    tp = tp0 + dtp * i;
    st->beta[i] = 1 / tp;
    st->lnz[i] = is2_exact(L, L, st->beta[i],
        &st->eref[i], &st->cvref[i]);
    is2_setuproba(st->beta[i], st->uproba[i]);
    printf("tp %7.4f, lnz %8.3f, eav %8.3f, cv %8.3f\n",
        tp, st->lnz[i], st->eref[i], st->cvref[i]);
  }

  for ( i = 0; i < st->nbeta; i++ )
    st->hist[i] = 0;

  xnew(st->p, st->nbeta);
  xnew(st->idmap, st->nbeta);
  xnew(st->ps, st->nbeta);
  xnew(st->cp, st->nbeta + 1);
  return st;
}


static void simtemp_close(simtemp_t *st)
{
  free(st->beta);
  free(st->hist);
  free(st->esum);
  free(st->eesum);

  free(st->lnz);
  free(st->uproba);
  free(st->eref);
  free(st->cvref);

  free(st->p);
  free(st->idmap);
  free(st->ps);
  free(st->cp);
  free(st);
}


/* choose a temperature according the current energy */
static int simtemp_choose(simtemp_t *st, double E, int ibet)
{
  int i;
  double max = -DBL_MAX;

  /* compute the probabilities */
  for ( i = 0; i < st->nbeta; i++ ) {
    st->p[i] = -st->beta[i] * E - st->lnz[i];
    if ( st->p[i] > max ) max = st->p[i];
  }
  /* no need to normalize to sum 1 */
  for ( i = 0; i < st->nbeta; i++ ) {
    st->p[i] = exp( st->p[i] - max );
  }

  return tc_select(st->type, st->nbeta, ibet, st->p,
      st->idmap, st->ps, st->cp);
}

static void simtemp(int type)
{
  simtemp_t *st;
  is2_t *is;
  int id, h, it, ibet;
  unsigned long t, nsteps = 1000000000L;
  double flatness, htot;

  printf("L %d\n", L);
  is = is2_open(L);
  st = simtemp_open(type, 4.0, 0.1, 5);

  ibet = 0;
  for ( t = 0; t < nsteps; t+= BATCHSIZE ) {
    /* run a batch of MC steps */
    for ( it = 0; it < BATCHSIZE; it++ ) {
      IS2_PICK(is, id, h);
      if ( h <= 0 || mtrand() <= st->uproba[ibet][h] ) {
        IS2_FLIP(is, id, h);
      }
    }
    ibet = simtemp_choose(st, is->E, ibet);
    if ( t < N * 100 ) continue;
    st->hist[ibet] += 1;
    st->esum[ibet] += is->E;
    st->eesum[ibet] += 1. * is->E * is->E;
  }

  for ( htot = 0, ibet = 0; ibet < st->nbeta; ibet++ ) {
    htot += st->hist[ibet];
  }

  flatness = 0;
  for ( ibet = 0; ibet < st->nbeta; ibet++ ) {
    double eav, cv;

    eav = st->esum[ibet] / st->hist[ibet];
    cv = st->eesum[ibet] / st->hist[ibet] - eav * eav;
    cv *= st->beta[ibet] * st->beta[ibet];
    flatness -= log(st->hist[ibet]*st->nbeta/htot)/st->nbeta;

    printf("H %8.5f, eav %8.3f/%8.3f, cv %8.3f/%8.3f\n",
        st->hist[ibet]/htot, eav, st->eref[ibet],
        cv, st->cvref[ibet]);
  }
  printf("flatness %g\n", flatness);
  is2_close(is);
  simtemp_close(st);
}


int main(int argc, char **argv)
{
  int type = 1;

  if ( argc > 1 ) {
    type = atoi( argv[1] );
  }

  mtscramble( time(NULL) );

  simtemp(type);
  return 0;
}

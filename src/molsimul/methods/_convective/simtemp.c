/* simulated tempering for Ising model */
#include "simtemp.h"



#ifndef IS2_LB
#define IS2_LB 5
#endif
#define L (1 << IS2_LB)
#define N ((L) * (L))
#include "is2.h"



/* Ising model parameters */
typedef struct {
  int n;
  double *beta;
  unsigned (*uproba)[5];
  double *lnzref;
  double *eref;
  double *cvref;
} is2param_t;



static is2param_t *is2param_open(int l, simtemp_t *st)
{
  is2param_t *isp;
  int i;

  xnew(isp, 1);
  isp->n = st->nbeta;
  /* link the temperature array */
  isp->beta = st->beta;
  xnew(isp->uproba, isp->n);
  xnew(isp->lnzref, isp->n);
  xnew(isp->eref,   isp->n);
  xnew(isp->cvref,  isp->n);
  for ( i = 0; i <  isp->n; i++ ) {
    isp->lnzref[i] = is2_exact(l, l, st->beta[i],
        &isp->eref[i], &isp->cvref[i]);
    is2_setuproba(isp->beta[i], isp->uproba[i]);
    printf("tp %7.4f, lnz %8.3f, eav %8.3f, cv %8.3f\n",
        1/isp->beta[i], isp->lnzref[i], isp->eref[i], isp->cvref[i]);
  }

  return isp;
}



static void is2param_close(is2param_t *isp)
{
  free(isp->uproba);
  free(isp->lnzref);
  free(isp->eref);
  free(isp->cvref);
  free(isp);
}



#define BATCHSIZE N



static void simtemp_is2(int type)
{
  simtemp_t *st;
  is2_t *is;
  is2param_t *isp;
  int id, h, it, ibet;
  unsigned long t, nsteps = 1000000000L;
  double flatness, htot;

  st = simtemp_open(type, 4.0, 0.1, 5);
  is = is2_open(L);
  isp = is2param_open(L, st);

  /* set the exact partition function array
   * for the weights of simulated tempering */
  simtemp_setlnz(st, isp->lnzref);

  printf("L %d, type %d, %s\n", L,
      (type & TC_MATRIXMASK),
      ((type & TC_FREQ) ? "freq" : "prob"));

  ibet = 0;
  for ( t = 0; t < nsteps; t+= BATCHSIZE ) {
    /* run a batch of MC steps */
    for ( it = 0; it < BATCHSIZE; it++ ) {
      IS2_PICK(is, id, h);
      if ( h <= 0 || mtrand() <= isp->uproba[ibet][h] ) {
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
        st->hist[ibet]/htot,
        eav, isp->eref[ibet],
        cv, isp->cvref[ibet]);
  }
  printf("flatness %g\n", flatness);

  is2param_close(isp);
  is2_close(is);
  simtemp_close(st);
}



int main(int argc, char **argv)
{
  int type = 1;

  if ( argc > 1 ) {
    type = atoi( argv[1] );
  }

  if ( argc > 2 ) { /* frequency-based sampling */
    if ( atoi( argv[2] ) > 0 ) {
      type |= TC_FREQ;
    }
  }

  mtscramble( time(NULL) );

  simtemp_is2(type);
  return 0;
}

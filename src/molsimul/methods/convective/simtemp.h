/* simulated tempering */
#ifndef SIMTEMP_H__
#define SIMTEMP_H__



#include "tconv.h"



typedef struct {
  int type;

  int nbeta;
  double *beta;
  double *hist;
  double *esum;
  double *eesum;

  double *lnz;

  double *p;
  int *idmap;
  double *sp;  /* sorted version of `p` */
  double *cp;  /* cumulative probabilities */

  /* pseudo-counts for frequency-based sampling */
  double *cnt;  /* unsorted */
  double *scnt; /* sorted */
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
  for ( i = 0; i < st->nbeta; i++ ) {
    tp = tp0 + dtp * i;
    st->beta[i] = 1 / tp;
  }

  for ( i = 0; i < st->nbeta; i++ )
    st->hist[i] = 0;

  xnew(st->p, st->nbeta);
  xnew(st->idmap, st->nbeta);
  xnew(st->sp, st->nbeta);
  xnew(st->cp, st->nbeta + 1);

  xnew(st->cnt, st->nbeta);
  xnew(st->scnt, st->nbeta);

  return st;
}



static void simtemp_close(simtemp_t *st)
{
  free(st->beta);
  free(st->hist);
  free(st->esum);
  free(st->eesum);

  free(st->lnz);

  free(st->p);
  free(st->idmap);
  free(st->sp);
  free(st->cp);

  free(st->cnt);
  free(st->scnt);

  free(st);
}



/* set the array for the partition function */
static void simtemp_setlnz(simtemp_t *st, const double *lnz)
{
  int i;

  for ( i = 0; i < st->nbeta; i++ ) {
    st->lnz[i] = lnz[i];
  }
}



/* choose a temperature according the current energy */
static int simtemp_choose(simtemp_t *st, double E, int ibet)
{
  int i;
  double max = -DBL_MAX;

  double psum = 0;

  /* compute the probabilities */
  for ( i = 0; i < st->nbeta; i++ ) {
    st->p[i] = -st->beta[i] * E - st->lnz[i];
    if ( st->p[i] > max ) max = st->p[i];
  }

  for ( i = 0; i < st->nbeta; i++ ) {
    st->p[i] = exp( st->p[i] - max );
    psum += st->p[i];
  }

  for ( i = 0; i < st->nbeta; i++ ) {
    st->p[i] /= psum;
  }

  if ( st->type & TC_FREQ ) {
    ibet = tc_next(st->type, st->nbeta, ibet, st->p,
          st->idmap, st->sp, st->cp, st->cnt, st->scnt);
  } else {
    ibet = tc_select(st->type, st->nbeta, ibet, st->p,
        st->idmap, st->sp, st->cp);
  }

  return ibet;
}



#endif /* SIMTEMP_H__ */


#ifndef WHAM_H__
#define WHAM_H__



/* Weighted histogram analysis method
 * this module requires `hist.h` */



#include "hist.h"
#include "mdiis.h"



typedef struct {
  double *beta; /* temperature array, reference */
  double *lnz; /* partition function, reference */
  double *lnz1;
  double *res;
  double *lndos; /* density of states */
  double *lntot;
  hist_t *hist; /* histograms, reference */
} wham_t;



static wham_t *wham_open(double *beta, double *lnz, hist_t *hist)
{
  wham_t *w;
  int i, j, nbeta = hist->rows, n = hist->n;

  xnew(w, 1);
  w->beta = beta;
  w->lnz = lnz;
  w->hist = hist;
  xnew(w->lntot, hist->rows);
  xnew(w->lnz1, hist->rows);
  xnew(w->res, hist->rows);
  xnew(w->lndos, hist->n);

  /* compute the total */
  for ( j = 0; j < nbeta; j++ ) {
    double x = 0;
    for ( i = 0; i < n; i++ )
      x += hist->arr[j*n+i];
    w->lntot[j] = log(x);
  }

  return w;
}



static void wham_close(wham_t *w)
{
  free(w->lnz1);
  free(w->res);
  free(w->lndos);
  free(w->lntot);
  free(w);
}



#define LOG0 -1e9



/* log(exp(a) + exp(b)) */
__inline static double wham_lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* save the density of states to file `fn` */
static int wham_savelndos(wham_t *w, const char *fn)
{
  FILE *fp;
  int i, n = w->hist->n;
  double emin = w->hist->xmin, de = w->hist->dx;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ )
    if ( w->lndos[i] > LOG0 )
      fprintf(fp, "%g %g\n", emin + (i+.5)*de, w->lndos[i]);
  fclose(fp);
  return 0;
}



/* compute the energy and heat capacity from the WHAM */
static void wham_getav(wham_t *w, const char *fn)
{
  hist_t *hist = w->hist;
  double T, lnZ, lnE, lnE2, Eav, Evar, enei, y, loge;
  double de = hist->dx, emin = hist->xmin;
  int i, j, n = hist->n, nbeta = hist->rows;
  FILE *fp = NULL;

  if ( fn != NULL && (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
  }
  for ( j = 0; j < nbeta; j++ ) {
    T = 1./w->beta[j];
    lnZ = lnE = lnE2 = LOG0;
    for ( i = 0; i < n; i++ ) {
      if (w->lndos[i] <= LOG0) continue;
      /* note: we do not add emin here for it may lead to
       * a negative energy whose logarithm is undefined */
      enei = (i + .5) * de;
      y = w->lndos[i] - w->beta[j] * enei;
      loge = log(enei);
      lnZ = wham_lnadd(lnZ, y);
      lnE = wham_lnadd(lnE, loge + y);
      lnE2 = wham_lnadd(lnE2, 2*loge + y);
    }
    Eav = exp(lnE - lnZ);
    Evar = exp(lnE2 - lnZ) - Eav * Eav;
    Eav += emin;
    if ( fp != NULL )
      fprintf(fp, "%g %g %g %g\n", T, Eav, Evar/(T*T), lnZ - w->beta[j] * emin);
  }
  if ( fp != NULL ) fclose(fp);
}



/* estimate the partition function using the single histogram method */
static void wham_estimatelnz(wham_t *w, double *lnz)
{
  hist_t *hist = w->hist;
  int i, j, n = hist->n, nbeta = hist->rows;
  double dbet, enei, hc, s, logsx;

  lnz[0] = 0;
  for ( j = 0; j < nbeta - 1; j++ ) {
    dbet = w->beta[j+1] - w->beta[j];
    s = 0;
    logsx = LOG0;
    for ( i = 0; i < n; i++ ) {
      hc = hist->arr[j*n + i];
      if ( hc <= 0 ) continue;
      enei = hist->xmin + (i + .5) * hist->dx; // no offset
      s += hc;
      logsx = wham_lnadd(logsx, log(hc) - dbet * enei);
    }
    if ( s <= 0 ) {
      lnz[j+1] = lnz[j];
      continue;
    }
    lnz[j+1] = lnz[j] + logsx - log(s);
    //printf("j %d, dlnZ %g, logsx %g, s %g\n", j, lnZ[j+1]-lnZ[j], logsx, s);
  }
}



static double wham_step(wham_t *w, double *lnz, double *res, int update)
{
  hist_t *hist = w->hist;
  int i, j, imin;
  int n = hist->n, nbeta = hist->rows;
  double x, num, den, enei, emin = hist->xmin, de = hist->dx, err;

  imin = -1;
  for ( i = 0; i < n; i++ ) {
    num = 0;
    den = LOG0;
    enei = emin + (i + .5) * de;
    /*        num           Sum_j arr[j, i]
     * dos = ----- = ------------------------------------------
     *        den     Sum_j tot[j] exp(-beta[j] * enei) / Z[j]
     * */
    for ( j = 0; j < nbeta; j++ ) {
      x = hist->arr[j*n + i];
      if ( x <= 0 ) continue;
      num += x;
      den = wham_lnadd(den, w->lntot[j] - w->beta[j] * enei - lnz[j]);
    }
    if ( num > 0 ) {
      w->lndos[i] = log(num) - den;
      if ( imin < 0 ) imin = i;
    } else {
      w->lndos[i] = LOG0;
    }
  }

  /* shift the origin of the density of states */
  for ( x = w->lndos[imin], i = n - 1; i >= 0; i-- )
    if ( w->lndos[i] > LOG0 )
      w->lndos[i] -= x;

  /* refresh the partition function */
  for ( j = 0; j < nbeta; j++ ) {
    for ( w->lnz1[j] = LOG0, i = 0; i < n; i++ ) {
      if ( w->lndos[i] <= LOG0) continue;
      enei = hist->xmin + (i + .5) * hist->dx;
      w->lnz1[j] = wham_lnadd(w->lnz1[j], w->lndos[i] - w->beta[j] * enei);
    }
  }
  for ( x = w->lnz1[0], j = 0; j < nbeta; j++ )
    w->lnz1[j] -= x; /* shift the baseline */

  for ( err = 0, j = 0; j < nbeta; j++ ) {
    res[j] = w->lnz1[j] - lnz[j];
    if ( fabs(res[j]) > err ) err = fabs(res[j]);
    if ( update ) lnz[j] = w->lnz1[j];
  }

  return err;
}



static double wham_getres(void *w, double *lnz, double *res)
{
  return wham_step((wham_t *) w, lnz, res, 0);
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double wham_getlndos(wham_t *w, double *lnz, int itmax, double tol)
{
  hist_t *hist = w->hist;
  int j, iter, nbeta = hist->rows;
  double err;

  for ( j = 0; j < nbeta; j++ )
    if ( fabs(lnz[j]) > DBL_MIN ) break;
  if ( j == nbeta )
    wham_estimatelnz(w, lnz);

  for ( iter = 1; iter <= itmax; iter++ ) {
    err = wham_step(w, lnz, w->res, 1);
    if (iter % 1000 == 0 )
      fprintf(stderr, "iter %d, err = %g\n", iter, err);
    if (err < tol) break;
  }
  fprintf(stderr, "partition function converged at step %d, error %g\n", iter, err);
  return err;
}



/* weighted histogram analysis method */
static void wham(hist_t *hist, double *beta, double *lnz,
    int itmax, double tol,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, lnz, hist);
  wham_getlndos(w, lnz, itmax, tol);
  if ( fnlndos ) wham_savelndos(w, fnlndos);
  wham_getav(w, fneav);
  wham_close(w);
}



static void wham_mdiis(hist_t *hist, double *beta, double *lnz,
    int nbases, double damp, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, lnz, hist);
  iter_mdiis(lnz, hist->rows, wham_getres, w, nbases, damp,
     itmax, tol, verbose);
  if ( fnlndos ) wham_savelndos(w, fnlndos);
  wham_getav(w, fneav);
  wham_close(w);
}



#endif /* WHAM_H__ */


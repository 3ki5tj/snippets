#ifndef WHAM_H__
#define WHAM_H__



/* Weighted histogram analysis method
 * this module requires `hist.h` */



#include "hist.h"



static int wham_savelndos(double *lndos,
    int n, double emin, double de, const char *fn);
static void wham_getav(double *lndos, int n, double emin, double de,
    double *beta, int nt, const char *fn);



#define LOG0 -1e9



/* log(exp(a) + exp(b)) */
__inline static double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* save the density of states to file `fn` */
static int wham_savelndos(double *lndos,
    int n, double emin, double de, const char *fn)
{
  FILE *fp;
  int i;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ )
    if (lndos[i] > LOG0)
      fprintf(fp, "%g %g\n", emin + (i+.5)*de, lndos[i]);
  fclose(fp);
  return 0;
}



/* compute the energy and heat capacity from the WHAM */
static void wham_getav(double *lndos, int n, double emin, double de,
    double *beta, int nt, const char *fn)
{
  double T, lnZ, lnE, lnE2, Eav, Evar, enei, y, loge;
  int i, j;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    exit(1);
  }
  for ( j = 0; j < nt; j++ ) {
    T = 1./beta[j];
    lnZ = lnE = lnE2 = LOG0;
    for ( i = 0; i < n; i++ ) {
      if (lndos[i] <= LOG0) continue;
      /* note: we do not add emin here for it may lead to
       * a negative energy whose logarithm is undefined */
      enei = (i + .5) * de;
      y = lndos[i] - beta[j] * enei;
      loge = log(enei);
      lnZ = lnadd(lnZ, y);
      lnE = lnadd(lnE, loge + y);
      lnE2 = lnadd(lnE2, 2*loge + y);
    }
    Eav = exp(lnE - lnZ);
    Evar = exp(lnE2 - lnZ) - Eav * Eav;
    Eav += emin;
    fprintf(fp, "%g %g %g %g\n", T, Eav, Evar/(T*T), lnZ - beta[j] * emin);
  }
  fclose(fp);
}



/* estimate the partition function using the single histogram method */
static void wham_estimatelnZ(hist_t *hist, double *beta, double *lnZ)
{
  int i, j, n = hist->n, nbeta = hist->rows;
  double dbet, enei, hc, s, logsx;

  lnZ[0] = 0;
  for ( j = 0; j < nbeta - 1; j++ ) {
    dbet = beta[j+1] - beta[j];
    s = 0;
    logsx = LOG0;
    for ( i = 0; i < n; i++ ) {
      hc = hist->arr[j*n + i];
      if ( hc <= 0 ) continue;
      enei = hist->xmin + (i + .5) * hist->dx; // no offset
      s += hc;
      logsx = lnadd(logsx, log(hc) - dbet * enei);
    }
    if ( s <= 0 ) {
      lnZ[j+1] = lnZ[j];
      continue;
    }
    lnZ[j+1] = lnZ[j] + logsx - log(s);
    //printf("j %d, dlnZ %g, logsx %g, s %g\n", j, lnZ[j+1]-lnZ[j], logsx, s);
  }
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double *wham_getlndos(hist_t *hist, double *beta, double *lnZ)
{
  double *lndos, *lntot, *lnZ0 = lnZ, *lnZ1, num, den, enei;
  double err, errmax, x;
  int i, imin, j, n = hist->n, nbeta = hist->rows;
  int iter;

  xnew(lndos, n);
  xnew(lntot, nbeta);
  if (lnZ == NULL)
    xnew(lnZ, nbeta);
  xnew(lnZ1, nbeta);

  for ( j = 0; j < nbeta; j++ ) {
    if ( fabs(lnZ[j]) > 1e-30 ) break;
  }
  if ( j == nbeta )
    wham_estimatelnZ(hist, beta, lnZ);

  /* compute the total */
  for ( j = 0; j < nbeta; j++ ) {
    num = 0;
    for ( i = 0; i < n; i++ )
      num += hist->arr[j*n+i];
    lntot[j] = log(num);
  }

#ifndef WHAM_ITERMAX
#define WHAM_ITERMAX 100000
#endif
#ifndef WHAM_LNZTOL
#define WHAM_LNZTOL 0.1
#endif

  for ( iter = 1; iter <= WHAM_ITERMAX; iter++ ) {
    imin = -1;
    for ( i = 0; i < n; i++ ) {
      num = 0;
      den = LOG0;
      enei = hist->xmin + (i + .5) * hist->dx;
      /*        num           Sum_j arr[j, i]
       * dos = ----- = ------------------------------------------
       *        den     Sum_j tot[j] exp(-beta[j] * enei) / Z[j]
       * */
      for ( j = 0; j < nbeta; j++ ) {
        x = hist->arr[j*n + i];
        if ( x <= 0 ) continue;
        num += x;
        den = lnadd(den, lntot[j] - beta[j] * enei - lnZ[j]);
      }
      if ( num > 0 ) {
        lndos[i] = log(num) - den;
        if ( imin < 0 ) imin = i;
      } else {
        lndos[i] = LOG0;
      }
    }

    /* shift the origin of the density of states */
    for ( i = n - 1; i >= 0; i-- )
      if (lndos[i] > LOG0)
        lndos[i] -= lndos[imin];

    /* refresh the partition function */
    for ( j = 0; j < nbeta; j++ ) {
      for ( lnZ1[j] = LOG0, i = 0; i < n; i++ ) {
        if (lndos[i] <= LOG0) continue;
        enei = hist->xmin + (i + .5) * hist->dx;
        lnZ1[j] = lnadd(lnZ1[j], lndos[i] - beta[j] * enei);
      }
    }
    errmax = 0;
    for ( i = -1, j = nbeta - 1; j >= 0; j-- ) {
      lnZ1[j] -= lnZ1[0]; /* shift the base */
      err = fabs(lnZ1[j] - lnZ[j]);
      if (err > errmax) { i = j; errmax = err; }
      lnZ[j] = lnZ1[j];
    }

    if (iter % 1000 == 0 ) {
      fprintf(stderr, "iter %d, err = %g, i %d\n", iter, errmax, i);
      //wham_savelndos(lndos, n, hist->xmin, hist->dx, "alndos.dat");
      //wham_getav(lndos, n, hist->xmin, hist->dx, beta, nbeta, "aeav.dat");
      //getchar();
    }
    if (errmax < WHAM_LNZTOL) break;
  }
  fprintf(stderr, "partition function converged at step %d, error %g\n", iter, errmax);
  free(lntot);
  if (lnZ0 == NULL) free(lnZ);
  free(lnZ1);
  return lndos;
}



/* weighted histogram analysis method */
static void wham(hist_t *hs, double *beta, double *lnZ,
    const char *fnlndos, const char *fneav)
{
  double *lndos = NULL;

  /* estimate the density of states */
  lndos = wham_getlndos(hs, beta, lnZ);
  if (fnlndos)
    wham_savelndos(lndos, hs->n, hs->xmin, hs->dx, fnlndos);
  wham_getav(lndos, hs->n, hs->xmin, hs->dx,
      beta, hs->rows, fneav);
  free(lndos);
}



#endif /* WHAM_H__ */


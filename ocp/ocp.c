/* one component plasma
 * To compile
 *  gcc ocp.c -lfftw3 -lm
 * To run the program
 *  ./a.out
 * To view the output
 *  gnuplot "cr.dat" u 1:2 t "c(r)", "" u 1:3 t "t(r)", "" u 1:($2+$3) t "h(r)"
 **/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <fftw3.h>
#define PI M_PI
#define xnew(x, n) { \
  if ((x = calloc((size_t)(n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %d\n", #x, (int) (n)); \
    exit(1); } }



double beta = 1.0;
int npt = 1024;
double rmax = 10.24;
double rho = 0.3;
double a = 1.0;
double damp = 1;
double tol = 1e-7;
int itmax = 1000;



/* compute
 *    out(k) = 2*fac/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, double *in, double *out, double fac,
    fftw_plan p, double *arr, double *ri, double *ki)
{
  int i;

  for ( i = 0; i < npt; i++ ) /* form in(x) * x */
    arr[i] = in[i] * ri[i];
  fftw_execute(p);
  for ( i = 0; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / ki[i];
}



/* solve the integral equation */
static void iter(int npt, double rho, double *ri, double *ki,
    double *dcr, double *dtr, double *dck, double *dtk,
    double *bphi, double *t0r, double *t0k,
    double facr2k, double fack2r,
    fftw_plan plan, double *arr, int itmax)
{
  int i, it;
  double x, err, errmax, ck;

  for ( it = 0; it < itmax; it++ ) {
    sphr(npt, dcr, dck, facr2k, plan, arr, ri, ki); /* dc(r) --> dc(k) */
    for ( i = 0; i < npt; i++ ) {
      ck = dck[i] - t0k[i]; /* add back the long range part */
      dtk[i] = rho * ck * ck / (1 - rho * ck) - t0k[i];
    }
    sphr(npt, dtk, dtr, fack2r, plan, arr, ki, ri); /* dt(k) --> dt(r) */
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = exp(-bphi[i] + t0r[i] + dtr[i]) - (1 + dtr[i]);
      if ((err = fabs(dcr[i] - x)) > errmax) errmax = err;
      dcr[i] += damp * (x - dcr[i]);
    }
    if ( errmax < tol ) break;
  }
}



int main(void)
{
  double dr, dk, facr2k, fack2r, surfr, surfk;
  double *bphi, *dcr, *dtr, *dck, *dtk, *t0r, *t0k;
  double *arr, *ri, *ki, *ri2, *ki2;
  int i;
  fftw_plan plan = NULL;
  FILE *fp;

  xnew(arr, npt);
  xnew(ri, npt);
  xnew(ki, npt);
  xnew(ri2, npt);
  xnew(ki2, npt);
  xnew(bphi, npt);
  xnew(dcr, npt);
  xnew(dtr, npt);
  xnew(dck, npt);
  xnew(dtk, npt);
  xnew(t0r, npt);
  xnew(t0k, npt);
  plan = fftw_plan_r2r_1d(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);

  dr = rmax / npt;
  dk = PI / (dr * npt);
  facr2k = PI*2 * dr;
  fack2r = pow(PI*2, -2) * dk;
  surfr = PI*4;
  surfk = surfr * pow(PI*2, -3);
  for ( i = 0; i < npt; i++ ) {
    ri[i] = dr * (i * 2 + 1) / 2;
    ki[i] = dk * (i * 2 + 1) / 2;
    ri2[i] = surfr * ri[i] * ri[i] * dr;
    ki2[i] = surfk * ki[i] * ki[i] * dk;
  }

  for ( i = 0; i < npt; i++ ) {
    bphi[i] = beta / ri[i];
    t0r[i] = bphi[i] * erf(ri[i]/sqrt(2)/a);
  }
  sphr(npt, t0r, t0k, facr2k, plan, arr, ri, ki); /* t0(r) --> t0(k) */

  for ( i = 0; i < npt; i++ ) /* initialize dc(r) for iteration */
    dcr[i] = -bphi[i] + t0r[i];

  iter(npt, rho, ri, ki, dcr, dtr, dck, dtk, bphi, t0r, t0k,
      facr2k, fack2r, plan, arr, itmax);

  /* output */
  if ((fp = fopen("cr.dat", "w")) != NULL) {
    for ( i = 0; i < npt; i++ ) {
      fprintf(fp, "%g %g %g %g %g %g %g %g\n",
          ri[i], dcr[i] - t0r[i], dtr[i] + t0r[i], t0r[i],
          ki[i], dck[i] - t0k[i], dtk[i] + t0k[i], t0k[i]);
    }
    fclose(fp);
  }

  free(bphi);
  free(dcr);
  free(dtr);
  free(dck);
  free(dtk);
  free(t0r);
  free(t0k);
  free(arr);
  free(ri);
  free(ki);
  free(ri2);
  fftw_cleanup();
  return 0;
}



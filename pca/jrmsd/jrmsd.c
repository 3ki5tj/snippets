#include "mat.h"
#include "eig.h"
#include "mtrand.h"


#if 1 /* dimer system */

#define N 2
double xref[][D] = {{-1, 0}, {1, 0}}; /* reference structure */
double x1[][D] = {{-1, 0}, {1, 0}}; /* test structure 1 */
double x2[][D] = {{-2, 0}, {2, 0}}; /* test structure 2 */
double x3[][D] = {{0, -4}, {0, 4}}; /* test structure 3 */

#else /* trimer system */

#define N 3
double xref[N][D] = {{-1, 0}, {0, 0}, {1, 0}};
double x1[N][D] = {{-1, 0}, {0, 0}, {1, 0}};
double x2[N][D] = {{-2, 0}, {0, 0}, {2, 0}};
double x3[N][D] = {{0, 1}, {0, 0}, {1, 0}};

#endif


int ntrials = 1000000;
double sigma = 0.1;



/* compute the Jacobian of fitting `z0` to `x0` */
static double jrmsd(int n, double x0[][D], double z0[][D],
    int nt, double sig)
{
  double (*x)[D], (*y)[D], (*y0)[D], *dy;
  double *mat, *eval, *evec;
  double rmsd0, rmsd, lnJ = 0;
  int i, d, nd = n * D, ii, jj, t;

  x = calloc(n, sizeof(x[0]));
  y = calloc(n, sizeof(y[0]));
  y0 = calloc(n, sizeof(y0[0]));
  dy = calloc(nd, sizeof(dy[0]));

  mat = calloc(nd * nd, sizeof(mat[0]));
  eval = calloc(nd, sizeof(eval[0]));
  evec = calloc(nd * nd, sizeof(evec[0]));

  /* first fit z0 to x0 */
  rmsd0 = vrmsd(z0, y0, x0, NULL, n, 0, NULL, NULL);

  for ( t = 0; t < nt; t++ ) {
    /* perturb x0 a bit to make x */
    for ( i = 0; i < n; i++ ) {
      for ( d = 0; d < D; d++ ) {
        x[i][d] = x0[i][d] + sig * randgaus();
      }
    }

    /* try to fit `z0` to x */
    rmsd = vrmsd(z0, y, x, NULL, n, 0, NULL, NULL);
    //printf("RMSD %g, %g\n", rmsd, rmsd0);

    /* compute the deviation */
    ii = 0;
    for ( d = 0; d < D; d++ ) {
      for ( i = 0; i < n; i++ ) {
        dy[ii++] = y[i][d] - y0[i][d];
      }
    }

    /* accumulate the variance matrix */
    for ( ii = 0; ii < nd; ii++ ) {
      for ( jj = 0; jj < nd; jj++ ) {
        mat[ii*nd + jj] += dy[ii] * dy[jj];
      }
    }
  }

  printf("covariance matrix:\n");
  for ( ii = 0; ii < nd; ii++ ) {
    for ( jj = 0; jj < nd; jj++ ) {
      mat[ii*nd + jj] /= nt;
      printf("%8.5f ", mat[ii*nd + jj]);
    }
    printf("\n");
  }

  /* compute the eigenvalues of mat */
  eigsym(mat, eval, evec, nd);

  /* compute the Jacobian from the nonzero eigenvalues */
  lnJ = -nd * log(sigma);
  for ( ii = 0; ii < nd; ii++ ) {
    if ( fabs(eval[ii]) > 1e-12 ) {
      lnJ += log(eval[ii]) * 0.5;
      printf("eigval %d: %g\n", ii, eval[ii]);
    }
  }
  printf("%d trials, RMSD %g, log Jacob %g\n\n",
      ntrials, rmsd0, lnJ);

  free(x);
  free(y);
  free(y0);
  free(dy);
  free(mat);
  free(eval);
  free(evec);

  return lnJ;
}



int main(void)
{
  double lnJ1, lnJ2, lnJ3;

  lnJ1 = jrmsd(N, x1, xref, ntrials, sigma);
  lnJ2 = jrmsd(N, x2, xref, ntrials, sigma);
  lnJ3 = jrmsd(N, x3, xref, ntrials, sigma);

  printf("lnJ1 %g, lnJ2 %g, lnJ3 %g\n", lnJ1, lnJ2, lnJ3);
  printf("ln(J2/J1)/ln(2) = %g\n", (lnJ2-lnJ1)/log(2));
  return 0;
}

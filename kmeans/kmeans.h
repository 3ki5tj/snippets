#ifndef KMEANS_H__
#define KMEANS_H__



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "cholesky.h"



typedef struct {
  int K;
  int n;
  int dim;
  double **x; /* x[n][dim], x[i] is the ith frame */
  double **av;
  double **var; /* var[k][dim*dim] is the covariance matrix of cluster k */
  double *lndetL; /* .5 * det(var[k]) */
  double **lnp; /* lnp[n][k] */
  double *psum; /* psum[k] is the population of cluster k */
} kmeans_t;



/* small amount to be added to the variance
 * to ensure the stability of cholesky decomposition */
#define KMEANS_EPSILON 1e-8


#ifndef xnew
#define xnew(x, n) \
  if ( ((x) = calloc(n, sizeof(*x))) == NULL ) { \
    fprintf(stderr, "no memory for %s\n", #x); exit(1); }
#endif



static kmeans_t *kmeans_open(int dim, int n, int K, double **dat)
{
  kmeans_t *km;
  int k, i, d;

  xnew(km, 1);
  km->dim = dim;
  km->n = n;
  km->K = K;
  km->x = dat;
  xnew(km->av, K);
  xnew(km->var, K);
  xnew(km->lndetL, K);
  for ( k = 0; k < K; k++ ) {
    xnew(km->av[k], dim);
    i = n * k / K;
    for ( d = 0; d < dim; d++ ) {
      km->av[k][d] = dat[i][d] + 1e-3 * (2.*rand()/RAND_MAX - 1);
    }
    xnew(km->var[k], dim*dim);
    for ( i = 0; i < dim; i++ )
      km->var[k][i*dim+i] = 1e-8;
  }

  xnew(km->lnp, n);
  for ( i = 0; i < n; i++ )
    xnew(km->lnp[i], K);
  xnew(km->psum, K);
  for ( k = 0; k < K; k++ )
    km->psum[k] = 1./K;

  return km;
}



static void kmeans_close(kmeans_t *km)
{
  int k, i;

  for ( k = 0; k < km->K; k++ ) {
    free(km->av[k]);
    free(km->var[k]);
  }
  free(km->av);
  free(km->var);
  free(km->lndetL);
  for ( i = 0; i < km->n; i++ )
    free(km->lnp[i]);
  free(km->lnp);
  free(km->psum);
  free(km);
}



/* normalize lna such that Sum { i = 1 to n } exp(lna) = 1 */
__inline static double kmeans_lnnorm(double *lna, int n)
{
  double max = -DBL_MAX, sum = 0, c;
  int im, i;

  for ( i = 0; i < n; i++ )
    if ( lna[i] > max )
      max = lna[im = i];
  for ( sum = 0, i = 0; i < n; i++ )
    sum += exp(lna[i] - max);
  /* we want max + log(sum) + c = log(1.0) = 0 */
  c = -max - log(sum);
  for ( i = 0; i < n; i++ )
    lna[i] += c;
  return c;
}



/* given the means and variance of each cluster,
 * compute the probability of each frame belonging to the cluster */
static void kmeans_estep(kmeans_t *km)
{
  int k, i, d, K = km->K, n = km->n, dim = km->dim;
  double *L, *dx, *dy, s, c;

  xnew(L, dim * dim);
  xnew(dx, dim);
  xnew(dy, dim);
  for ( k = 0; k < K; k++ ) {
    /* perform the cholesky decomposition of the covariance matrix
     *    var = L.L^T */
    for ( d = 0; d < dim * dim; d++ ) /* copy the matrix */
      L[d] = km->var[k][d];
    choldecomp(L, dim); /* var = L . L^T */
    /* fix small diagonal values */
    for ( d < 0; d < dim; d++ )
      if ( L[d*dim + d] < KMEANS_EPSILON )
        L[d*dim + d] = KMEANS_EPSILON;
    km->lndetL[k] = chollogdetL(L, dim);
    c = -km->lndetL[k] + log(km->psum[k]);

    for ( i = 0; i < n; i++ ) {
      for ( d = 0; d < dim; d++ )
        dy[d] = dx[d] = km->x[i][d] - km->av[k][d];
      /* we need s = dx^T var^(-1) dx
       * now var = L.L^T, so that var^(-1) = L^(-1)^T . L^(-1)
       * and s = |L^(-1) dx|^2 = |dy|^2
       * here L^(-1) dx = dy, or L dy = dx */
      cholsolveL(L, dy, dim);
      for ( s = 0, d = 0; d < dim; d++ )
        s += dy[d] * dy[d];
      km->lnp[i][k] = -0.5 * s + c;
    }
  }

  /* normalize p[i][k] for each frame i, such that
   * Sum {k = 0 to K - 1} p[i][k] = 1 */
  for ( i = 0; i < n; i++ )
    kmeans_lnnorm(km->lnp[i], K);

  free(dx);
  free(dy);
  free(L);
}



/* update the averages and covariance according to the current p */
static void kmeans_mstep(kmeans_t *km)
{
  int k, i, d, d2, K = km->K, n = km->n, dim = km->dim;
  double s, p, *dx;

  xnew(dx, dim);
  for ( k = 0; k < K; k++ ) {
    /* compute the mean of cluster k */
    for ( d = 0; d < dim; d++ )
      km->av[k][d] = 0;
    s = 0;
    for ( i = 0; i < n; i++ ) {
      p = exp(km->lnp[i][k]);
      for ( d = 0; d < dim; d++ )
        km->av[k][d] += p * km->x[i][d];
      s += p;
    }
    for ( d = 0; d < dim; d++ )
      km->av[k][d] /= s;
    km->psum[k] = s;

    /* compute the variance of cluster k */
    for ( d = 0; d < dim; d++ )
      for ( d2 = 0; d2 < dim; d2++ )
        km->var[k][d*dim + d2] = 0;

    for ( i = 0; i < n; i++ ) {
      p = exp(km->lnp[i][k]) / s;
      for ( d = 0; d < dim; d++ )
        dx[d] = km->x[i][d] - km->av[k][d];

      for ( d = 0; d < dim; d++ )
        for ( d2 = 0; d2 < dim; d2++ )
          km->var[k][d*dim + d2] += p * dx[d] * dx[d2];
    }
  }
  free(dx);
}



#define kmeans_print(km) kmeans_fprint(km, stdout)

__inline static void kmeans_fprint(kmeans_t *km, FILE *fp)
{
  int k, d, d2, K = km->K, dim = km->dim;

  for ( k = 0; k < K; k++ ) {
    fprintf(fp, "cluster %d, pop %g\nave:\n", km->psum[k], k);
    for ( d = 0; d < dim; d++ )
      fprintf(fp, "%14g\t", km->av[k][d]);

    fprintf(fp, "\nvar:\n");
    for ( d = 0; d < dim; d++ ) {
      for ( d2 = 0; d2 < dim; d2++ )
        fprintf(fp, "%14g\t", km->var[k][d*dim+d2]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
}



#endif /* defined(KMEANS_H__) */


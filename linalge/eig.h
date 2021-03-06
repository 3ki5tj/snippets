#ifndef EIG_H__
#define EIG_H__



/* diagonalization of a symmetric real matrix */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>



/* return sqrt(a*a + b*b) */
__inline static double dblhypot(double x, double y)
{
  return sqrt(x*x + y*y);
}


/* To reduce a real symmetric matrix 'm' to tridiagonal by Householder transformations.
 * The diagonal elements are saved in vector 'd' and off-diagonal elements 'e'.  */
static void tridiag(double *m, double d[], double e[], int n)
{
  int i, j, k;
  double H, sigma, p, K, *x;

  /* use d[i] to indicate if the i'th Householder transformation is performed */
  for (i = 0; i < n; i++) d[i] = 0;

  /* n-2 Householder transformations */
  for (i = 0; i < n - 2; i++) {
    x = m + i*n; /* alias x[k] == m[i*n+k] */

    for (H = 0, k = i + 1; k < n; k++) H += x[k]*x[k];
    sigma = x[i + 1] > 0 ? sqrt(H) : -sqrt(H); /* sigma = sgn(x1) |x| */
    e[i] = -sigma; /* P x = - sigma e1 */
    H += sigma*x[i + 1]; /* H= (1/2) |u|^2 = |x|^2 + sigma x1 */

    /* To avoid singularity due to (partially) diagonal matrix as input */
    p = fabs(sigma);
    if (p <= fabs(m[i*n + i]) * DBL_EPSILON || p <= DBL_MIN) {
      e[i] = m[i*n + i + 1];
      continue;
    }

    x[i + 1] += sigma;  /* u = x + sigma e1, we now switch to 'u' */
    for (j = i + 1; j < n; j++) m[j*n + i] = x[j]/H; /* save u/H in column i */

    /* calculate P A P */
    K = 0;
    for (j = i + 1; j < n; j++) {
      /* calculate p=A u /H, we only use the up triangle */
      for (p = 0, k = i + 1; k <= j; k++) p += m[k*n + j]*x[k];
      for (k = j + 1; k < n; k++) p += m[j*n + k]*x[k];
      e[j] = (p /= H); /* save p temporarily to e[j], notice e[i+1..n-1] are not used yet.*/
      K += x[j]*p; /* K = u' p / (2H) */
    }
    K /= (2*H); /* K = u' p / (2H) */
    for (j = i + 1; j < n; j++) e[j] -= K*x[j];  /* form  q = p - K u */
    for (j = i + 1; j < n; j++) /* calculate A' = A - q u' - u q' (only right-top triangle) */
      for (k = j; k < n; k++)
        m[j*n + k] -= e[j]*x[k] + x[j]*e[k];

    d[i] = 1; /* indicate that the transformation is performed */
  }
  if ( n > 1 ) e[n - 2] = m[(n - 2)*n + n - 1]; /* for i == n-2 */
  e[n - 1] = 0;

  /* if only eigenvalues are required, enable the above line and ignore the rest */

  /* To form Q = P1 ... Pn-2 */
  /* copy last two eigenvalues */
  if ( n > 1 )
    d[n - 2] = m[(n - 2)*n + n - 2];
  d[n - 1] = m[(n - 1)*n + n - 1];
  /* initialize the right-bottom corner */
  if ( n > 1 ) {
    m[(n - 2)*n + n - 2] = 1;
    m[(n - 2)*n + n - 1] = 0;
    m[(n - 1)*n + n - 2] = 0;
  }
  m[(n - 1)*n + n - 1] = 1;

  /* P Q = (1 - u u'/H) Q = Q - (u/H) (u' Q) */
  for (i = n - 3; i >= 0; i--) { /* for each P */
    x = m + i*n; /* alias x[k] == m[i*n+k] */

    /* form eigenvector, ONLY if i'th transformation is performed */
    if (fabs(d[i]) > DBL_MIN) {
      for (j = i + 1; j < n; j++) {
        /* form K = u'Q */
        for (K = 0, k = i + 1; k < n; k++) K += x[k]*m[k*n + j];
        /* Q = Q - K (u/H)  */
        for (k = i + 1; k < n; k++) m[k*n + j] -= K*m[k*n + i];
      }
    }
    /* copy the diagonal element and proceed */
    d[i] = m[i*n + i];
    m[i*n + i] = 1;
    for (j = i + 1; j < n; j++) m[i*n + j] = m[j*n + i] = 0.;
  }
}


/* diagonalize the tridiagonal matrix by QR algorithm,
 * whose diagonal is d[0..n-1], off-diagonal is e[0..n-2];
 * the eigenvalues are saved in d[0..n-1] on return */
static void eigtriqr(double d[], double e[], int n, double *mat)
{
  const int itermax = 1000;
  int i, j, k, m, iter, sgn;
  double ks = 0, r, c, s, delta, f = 0, t1, t2, tol;

  e[n - 1] = 0;
  tol = 1e-12f;
  /* reduction starts from the left-top corner to the right-bottom corner */
  for ( i = 0; i < n; i++ ) {
    /* for each eigenvalue */
    for ( iter = 0; iter < itermax; iter++ ) {
      /* Look for a single small subdiagonal element to split the matrix */
      for (m = i; m < n - 1; m++) {
        double d1 = fabs(d[m + 1]), d2 = fabs(d[m]);
        if (fabs(e[m]) < (d1 + d2) * tol)
          break;
      }

      /* I have isolated d[i] from other matrix elements
       * so that d[i] is the eigenvalue.
       * stop iteration and look for next(i+1) eigenvalue  */
      if (m == i) break;

      /* form shift ks */
      delta = d[m] - d[m - 1];
      sgn = ((delta > 0) ? 1: -1);
      delta /= e[m - 1];
      r = dblhypot(delta, 1);
      ks = d[m] + sgn*e[m - 1]/(r + fabs(delta));

      /* rotations */
      for (j = i; j <= m - 1; j++) {
        /* calculate c and s */
        if (j == i) {
          /* First rotation */
          r = dblhypot(d[i] - ks, e[i]);
          c = (d[i] - ks)/r;
          s = e[i]/r;
        } else {
          /* Givens rotations */
          r = dblhypot(e[j - 1], f);
          c = e[j - 1]/r;
          s = f/r;
          e[j - 1] = r;
        }

        /* update the diagonal and off-diagonal elements */
        r = s*(d[j + 1] - d[j]) + 2*c*e[j];
        d[j]     += s*r;
        d[j + 1] -= s*r;
        e[j]      = c*r - e[j];
        f         = s*e[j + 1];
        e[j + 1] *= c;

        /* update eigenvectors */
        for (k = 0; k < n; k++) {
          t1 = mat[k*n + j];
          t2 = mat[k*n + j + 1];
          mat[k*n + j]   = c*t1 + s*t2;
          mat[k*n + j + 1] = -s*t1 + c*t2;
        }
      } /* end of rotations */
    } /* end for iteration */
    /*printf("Iterate %d times for %d'th eigenvalue.\n", iter, i);*/
  }/* end for each eigenvalue */
}



/* sort eigenvalues and eigenvectors in ascending order */
static void eigsort(double *d, double *v, int n)
{
  int i, j, im;
  double max, tmp;

  for (i = 0; i < n - 1; i++) {
    /* find the largest eigenvalue */
    for (max = d[i], im = i, j = i + 1; j < n; j++) {
      if (d[j] > max) max = d[im = j];
    }
    if (im != i) { /* swap columns im and i */
      tmp = d[i], d[i] = d[im], d[im] = tmp;
      for (j = 0; j < n; j++)
        tmp = v[j*n + i], v[j*n + i] = v[j*n + im], v[j*n + im] = tmp;
    }
  }
}



/* solve eigensystem of a real symmetric matrix `mat',
 * eigenvalues saved to `d', eigenvectors to v */
static int eigsym(double *mat, double *d, double *v, int n)
{
  double *e;
  int i;

  if ( (e = calloc(n, sizeof(*e))) == NULL ) {
    exit(1);
  }
  for ( i = 0; i < n * n; i++ ) {
    v[i] = mat[i];
  }
  tridiag(v, d, e, n);
  eigtriqr(d, e, n, v);
  eigsort(d, v, n);
  free(e);
  return 0;
}



#endif /* EIG_H__ */


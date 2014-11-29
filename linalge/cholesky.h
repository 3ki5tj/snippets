#ifndef CHOLESKY_H__
#define CHOLESKY_H__



/* Cholesky decomposition */



/* compute the Cholesky decomposition
 * the input matrix `a' should be positive definite
 * on return the left-bottom triangle is filled by the matrix L */
int choldecomp(double *a, int n)
{
  int i, j, k;
  double y;

  for ( i = 0; i < n; i ++ ) {
    for ( j = 0; j <= i; j++ ) {
      /* A(i, j) = L(i, 0) L(j, 0) + L(i, 1) L(j, 1) + ... + L(i, j) L(j, j)
       *         = sum {k = 0 to j - 1} L(i, k) L(j, k) + L(i, j) L(j, j) */
      for ( y = a[i*n + j], k = 0; k < j; k++ )
        y -= a[i*n + k] * a[j*n + k];
      if ( i == j ) {
        if ( y < 0 ) {
          fprintf(stderr, "cholesky: negative element on %d\n", i);
          return -1;
        }
        a[i*n + i] = sqrt(y);
      } else {
        a[j*n + i] = a[i*n + j] = y / a[j*n + j];
      }
    }
  }
  return 0;
}



/* solve L x = b, with a = L.L^T by Cholesky decomposition
 * on return, `x' is saved in `b' */
__inline static void cholsolveL(double *L, double *b, int n)
{
  int i, j;
  double y;

  /* sum {j = 0 to i - 1} L(i, j) y(j) + L(i, i) y(i) = b(i) */
  for ( i = 0; i < n; i++ ) {
    for ( y = b[i], j = 0; j < i; j++ )
      y -= L[i*n + j] * b[j];
    b[i] = y / L[i*n + i];
  }
}



/* solve L^T x = b
 * on return, `x' is saved in `b' */
__inline static void cholsolveLT(double *L, double *b, int n)
{
  int i, j;
  double y;

  for ( i = n - 1; i >= 0; i-- ) {
    for ( y = b[i], j = i + 1; j < n; j++ )
      y -= L[i*n + j] * b[j]; /* a(i, j) == a(j, i) */
    b[i] = y / L[i*n + i];
  }
}



/* solve a x = b by Cholesky decomposition
 * on return, `x' is saved in `b' */
int cholsolve(double *a, double *b, int n)
{
  if ( choldecomp(a, n) != 0 ) return -1;
  cholsolveL(a, b, n); /* solve L y = b */
  cholsolveLT(a, b, n); /* solve L^T x = y */
  return 0;
}



/* inverse the matrix `a', b = a^(-1) by Cholesky decomposition */
__inline static int cholinv(double *a, double *b, int n)
{
  int i, j, k;
  double y;

  if ( choldecomp(a, n) != 0 ) return -1;
  for ( k = 0; k < n; k++ ) {
    /* solve L y = b_k, with b_k(i) = delta(k, i)
     * sum {j = 0 to i - 1} L(i, j) y(j) + L(i, i) y(i) = b(i) */
    for ( i = 0; i < n; i++ ) {
      y = (i == k) ? 1 : 0;
      for ( j = 0; j < i; j++ )
        y -= a[i*n + j] * b[j*n + k];
      b[i*n + k] = y / a[i*n + i];
    }
    /* solve L^T x = y */
    for ( i = n - 1; i >= 0; i-- ) {
      for ( y = b[i*n + k], j = i + 1; j < n; j++ )
        y -= a[i*n + j] * b[j*n + k]; /* a(i, j) == a(j, i) */
      b[i*n + k] = y / a[i*n + i];
    }
  }
  return 0;
}



/* return determinant of L
 * note that the determinant of the original matrix is twice as large */
__inline static double chollogdetL(double *L, int n)
{
  int i;
  double y = 0;

  for ( i = 0; i < n; i++ ) y += log(L[i*n+i]);
  return y;
}

#endif

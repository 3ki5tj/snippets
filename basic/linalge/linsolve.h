#include <math.h>
#include <float.h>


/* solve the linear equation: a x = b */
static int linsolve(double *a, double *x, double *b, int n)
{
  int i, j, k, ip;
  double y;

  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    /* 1. select the pivot of the ith column
     * pivot: the maximal element a(j = i..n-1, i)  */
    y = fabs( a[(ip = i)*n + i] );
    for ( j = i + 1; j < n; j++ )
      if ( fabs( a[j*n + i] ) > y )
        y = fabs( a[(ip = j)*n + i] );

    /* 2. swap the pivot (ip'th) row with the ith row */
    if ( ip != i ) {
      y = b[ip], b[ip] = b[i], b[i] = y;
      for ( j = i; j < n; j++ )
        y = a[i*n + j], a[i*n + j] = a[ip*n + j], a[ip*n + j] = y;
    }
    y = a[i*n + i];
    if ( fabs(y) < DBL_MIN ) {
      fprintf(stderr, "singular matrix on %dth row\n", i);
      return -1;
    }

    /* 3. normalize the ith row */
    b[i] /= y;
    for ( k = i; k < n; k++ )
      a[i*n + k] /= y;

    /* 4. use the pivot row to eliminate the following rows */
    for ( j = i + 1; j < n; j++ ) { /* for rows */
      y = a[j*n + i];
      b[j] -= y * b[i];
      for ( k = i; k < n; k++ )
        a[j*n + k] -= y * a[i*n + k];
    }
  }

  /* 5. now that the matrix is upper-triangular
   *    solve for x */
  for ( i = n - 1; i >= 0; i-- ) {
    x[i] = b[i] / a[i*n + i];
    for ( j = 0; j < i; j++ ) {
      b[j] -= a[j*n + i] * x[i];
    }
  }
  return 0;
}


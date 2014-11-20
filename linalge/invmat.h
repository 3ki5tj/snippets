/* compute the inverse matrix b = a^(-1) by Gaussian elimination */
static int invmat(double *a, double *b, int n)
{
  int i, j, k, ip;
  double x;

  /* initialize b as the identity matrix */
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ )
      b[i*n + j] = (i == j);

  /* Gaussian elimination */
  for ( i = 0; i < n; i++ ) {
    /* choose the pivot as the largest element of column i */
    x = fabs( a[(ip = i)*n + i] );
    for ( k = ip + 1; k < n; k++ )
      if ( fabs( a[k*n + i] ) > x )
        x = fabs( a[(ip = k)*n + i] );

    /* swap the pivot (ip'th) row with the present row i */
    for ( k = i; k < n; k++ )
      x = a[i*n + k], a[i*n + k] = a[ip*n + k], a[ip*n + k] = x;
    for ( k = 0; k < n; k++ )
      x = b[i*n + k], b[i*n + k] = b[ip*n + k], b[ip*n + k] = x;

    /* normalize this row */
    x = a[i*n + i];
    if ( fabs(x) < DBL_MIN ) {
      fprintf(stderr, "Error: singular matrix of %dx%d\n", n, n);
      return -1;
    }
    for ( k = i; k < n; k++ ) a[i*n + k] /= x;
    for ( k = 0; k < n; k++ ) b[i*n + k] /= x;

    /* use the pivot row to zero the rest rows */
    for ( j = i + 1; j < n; j++ ) {
      x = a[j*n + i];
      for ( k = i; k < n; k++ )
        a[j*n + k] -= x * a[i*n + k];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= x * b[i*n + k];
    }
  }

  /* now that the matrix is upper triangular
   * make it diagonal */
  for ( i = n - 1; i >= 0; i-- ) {
    /* note a[i*n + i] should be 1 now */
    for ( j = 0; j < i; j++ ) {
      x = a[j*n + i];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= b[i*n + k] * x;
    }
  }
  return 0;
}





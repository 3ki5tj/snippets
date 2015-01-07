/* Cholesky decomposition */



"use strict";



/* compute the Cholesky decomposition
 * the input matrix `a' should be positive definite
 * on return the left-bottom triangle is filled by the matrix L */
function choldecomp(a, n)
{
  for ( var i = 0; i < n; i ++ ) {
    for ( var j = 0; j <= i; j++ ) {
      var y = a[i*n + j];
      // A(i, j) = L(i, 0) L(j, 0) + L(i, 1) L(j, 1) + ... + L(i, j) L(j, j)
      //         = sum {k = 0 to j - 1} L(i, k) L(j, k) + L(i, j) L(j, j)
      for ( var k = 0; k < j; k++ )
        y -= a[i*n + k] * a[j*n + k];
      if ( i == j ) {
        if ( y < 0 ) {
          console.log("cholesky: negative element, i", i, "y", y);
          y = 0;
        }
        a[i*n + i] = Math.sqrt(y);
      } else {
        a[j*n + i] = a[i*n + j] = y / a[j*n + j];
      }
    }
  }
  return 0;
}



/* solve L x = b, with a = L.L^T by Cholesky decomposition
 * on return, `x' is saved in `b' */
function cholsolveL(L, b, n)
{
  // sum {j = 0 to i - 1} L(i, j) y(j) + L(i, i) y(i) = b(i)
  for ( var i = 0; i < n; i++ ) {
    var y = b[i];
    for ( var j = 0; j < i; j++ )
      y -= L[i*n + j] * b[j];
    b[i] = y / L[i*n + i];
  }
}



/* solve L^T x = b
 * on return, `x' is saved in `b' */
function cholsolveLT(L, b, n)
{
  for ( var i = n - 1; i >= 0; i-- ) {
    var y = b[i];
    for ( var j = i + 1; j < n; j++ )
      y -= L[i*n + j] * b[j]; // a(i, j) == a(j, i)
    b[i] = y / L[i*n + i];
  }
}



/* solve a x = b by Cholesky decomposition
 * on return, `x' is saved in `b' */
function cholsolve(a, b, n)
{
  if ( choldecomp(a, n) != 0 ) return -1;
  cholsolveL(a, b, n); // solve L y = b
  cholsolveLT(a, b, n); // solve L^T x = y
  return 0;
}



/* inverse the matrix `a', b = a^(-1) by Cholesky decomposition */
function cholinv(a, b, n)
{
  var i, j, k;
  var y;

  if ( choldecomp(a, n) != 0 ) return -1;
  for ( k = 0; k < n; k++ ) {
    // solve L y = b_k, with b_k(i) = delta(k, i)
    // sum {j = 0 to i - 1} L(i, j) y(j) + L(i, i) y(i) = b(i)
    for ( i = 0; i < n; i++ ) {
      y = (i == k) ? 1 : 0;
      for ( j = 0; j < i; j++ )
        y -= a[i*n + j] * b[j*n + k];
      b[i*n + k] = y / a[i*n + i];
    }
    // solve L^T x = y
    for ( i = n - 1; i >= 0; i-- ) {
      for ( y = b[i*n + k], j = i + 1; j < n; j++ )
        y -= a[i*n + j] * b[j*n + k]; // a(i, j) == a(j, i)
      b[i*n + k] = y / a[i*n + i];
    }
  }
  return 0;
}



/* return determinant of L
 * note that the determinant of the original matrix is twice as large */
function chollogdetL(L, n, eps)
{
  var y = 0;

  for ( var i = 0; i < n; i++ )
    y += Math.log( Math.max(L[i*n+i], eps) );
  return y;
}


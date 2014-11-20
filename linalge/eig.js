/* return sqrt(a*a + b*b) */
function dblhypot(x, y)
{
  return Math.sqrt(x*x + y*y);
}


/* To reduce a real symmetric matrix 'm' to tridiagonal by Householder transformations.
 * The diagonal elements are saved in vector 'd' and off-diagonal elements 'e' */
function tridiag(m, d, e, n)
{
  var i, j, k;
  var H, sigma, p, K;

  // use d[i] to indicate if the i'th Householder transformation is performed
  for (i = 0; i < n; i++) d[i] = 0;

  // n-2 Householder transformations
  for (i = 0; i < n - 2; i++) {
    for (H = 0, k = i + 1; k < n; k++)
      H += m[i*n+k]*m[i*n+k];
    sigma = m[i*n+i+1] > 0 ? Math.sqrt(H) : -Math.sqrt(H); // sigma = sgn(x1) |m[i]|
    e[i] = -sigma; // P m[i] = - sigma e1
    H += sigma*m[i*n+i+1]; // H= (1/2) |u|^2 = |m[i]|^2 + sigma x1

    // To avoid singularity due to (partially) diagonal matrix as input
    if (Math.abs(sigma) <= Math.abs(m[i*n+i])*1e-14) {
      e[i] = m[i*n+i+1];
      continue;
    }

    m[i*n+i+1] += sigma;  // u = m[i] + sigma e1, we now switch to 'u'
    for (j = i + 1; j < n; j++) m[j*n + i] = m[i*n + j]/H; // save u/H in column i

    // calculate P A P
    K = 0;
    for (j = i + 1; j < n; j++) {
      // calculate p=A u /H, we only use the up triangle
      for (p = 0, k = i + 1; k <= j; k++)
        p += m[k*n + j]*m[i*n + k];
      for (k = j + 1; k < n; k++)
        p += m[j*n + k]*m[i*n + k];
      e[j] = (p /= H); // save p temporarily to e[j], notice e[i+1..n-1] are not used yet
      K += m[i*n + j]*p; // K = u' p / (2H)
    }
    K /= (2*H); // K = u' p / (2H)
    for (j = i + 1; j < n; j++)
      e[j] -= K*m[i*n+j];  // form  q = p - K u
    // calculate A' = A - q u' - u q' (only right-top triangle)
    for (j = i + 1; j < n; j++)
      for (k = j; k < n; k++)
        m[j*n + k] -= e[j]*m[i*n + k] + m[i*n + j]*e[k];

    d[i] = 1; // indicate that the transformation is performed
  }
  if ( n > 1 ) e[n - 2] = m[(n - 2)*n + n - 1]; // for i == n-2
  e[n - 1] = 0;

  // if only eigenvalues are required, enable the above line and ignore the rest

  // To form Q = P1 ... Pn-2
  // copy last two eigenvalues
  if ( n > 1 )
    d[n - 2] = m[(n - 2)*n + n - 2];
  d[n - 1] = m[(n - 1)*n + n - 1];
  // initialize the right-bottom corner
  if ( n > 1 ) {
    m[(n - 2)*n + n - 2] = 1;
    m[(n - 2)*n + n - 1] = 0;
    m[(n - 1)*n + n - 2] = 0;
  }
  m[(n - 1)*n + n - 1] = 1;

  // P Q = (1 - u u'/H) Q = Q - (u/H) (u' Q)
  for (i = n - 3; i >= 0; i--) { // for each P
    // form eigenvector, ONLY if i'th transformation is performed
    if (Math.abs(d[i]) > 1e-30) {
      for (j = i + 1; j < n; j++) {
        // form K = u'Q
        for (K = 0, k = i + 1; k < n; k++)
          K += m[i*n + k]*m[k*n + j];
        // Q = Q - K (u/H)
        for (k = i + 1; k < n; k++)
          m[k*n + j] -= K*m[k*n + i];
      }
    }
    // copy the diagonal element and proceed
    d[i] = m[i*n + i];
    m[i*n + i] = 1;
    for (j = i + 1; j < n; j++) m[i*n + j] = m[j*n + i] = 0.;
  }
}


/* diagonalize the tridiagonal matrix by QR algorithm,
 * whose diagonal is d[0..n-1], off-diagonal is e[0..n-2];
 * reduce from the left-top to right-left */
function eigtriqr(d, e, n, mat)
{
  var itermax = 100;
  var i, j, k, m, iter, sgn;
  var ks = 0, r, c, s, delta, f = 0, t1, t2, tol;

  e[n - 1] = 0;
  tol = 1e-12;
  for (i = 0; i < n; i++) {
    // for each eigenvalue
    for (iter = 0; iter < itermax; iter++) {
      // Look for a single small subdiagonal element to split the matrix
      for (m = i; m < n - 1; m++) {
        var d1 = Math.abs(d[m + 1]);
        var d2 = Math.abs(d[m]);
        if (Math.abs(e[m]) < (d1 + d2) * tol)
          break;
      }

      // I have isolated d[i] from other matrix elements
      // so that d[i] is the eigenvalue.
      // stop iteration and look for next(i+1) eigenvalue
      if (m == i) break;

      // form shift ks
      delta = d[m] - d[m - 1];
      sgn = ((delta > 0) ? 1: -1);
      delta /= e[m - 1];
      r = dblhypot(delta, 1);
      ks = d[m] + sgn*e[m - 1]/(r + Math.abs(delta));

      // rotations
      for (j = i; j <= m - 1; j++) {
        // calculate c and s
        if (j == i) {
          // first rotation
          r = dblhypot(d[i] - ks, e[i]);
          c = (d[i] - ks)/r;
          s = e[i]/r;
        } else {
          // Givens rotations
          r = dblhypot(e[j - 1], f);
          c = e[j - 1]/r;
          s = f/r;
          e[j - 1] = r;
        }

        // update the diagonal and off-diagonal elements
        r = s*(d[j + 1] - d[j]) + 2*c*e[j];
        d[j]     += s*r;
        d[j + 1] -= s*r;
        e[j]      = c*r - e[j];
        f         = s*e[j + 1];
        e[j + 1] *= c;

        // update eigenvectors
        for (k = 0; k < n; k++) {
          t1 = mat[k*n + j];
          t2 = mat[k*n + j + 1];
          mat[k*n + j]   = c*t1 + s*t2;
          mat[k*n + j + 1] = -s*t1 + c*t2;
        }
      } // end of rotations
    } // end for iteration
    if ( iter >= itermax )
      console.log("failed to converge after " + iter + " iterations");
  }// end for each eigenvalue
}



/* sort eigenvalues and eigenvectors in ascending order */
function eigsort(d, v, n)
{
  var i, j, im;
  var max, tmp;

  for (i = 0; i < n - 1; i++) {
    /* search the maximal eigenvalue */
    for (max = d[i], im = i, j = i + 1; j < n; j++) {
      if (d[j] > max) max = d[im = j];
    }
    if (im != i) { /* change column im and i */
      tmp = d[i], d[i] = d[im], d[im] = tmp;
      for (j = 0; j < n; j++)
        tmp = v[j*n + i], v[j*n + i] = v[j*n + im], v[j*n + im] = tmp;
    }
  }
}



/* solve eigensystem of a real symmetric matrix `mat',
 * eigenvalues saved to `d', eigenvectors to v */
function eigsym(mat, d, v, n)
{
  var e = new Array(n);
  var i;

  for (i = 0; i < n*n; i++) v[i] = mat[i];
  tridiag(v, d, e, n);
  eigtriqr(d, e, n, v);
  eigsort(d, v, n);
}




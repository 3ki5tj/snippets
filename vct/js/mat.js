


"use strict";



/* return a unit matrix */
function munit(a)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      a[i][j] = (i === j) ? 1.0 : 0.0;
    }
  }
}



/* a = b */
function mcopy(a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      a[i][j] = b[i][j];
    }
  }
}



/* a = b^T */
function mtrans(a)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = i + 1; j < D; j++ ) {
      var x = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = x;
    }
  }
}



/* c = a^T b */
function mvtxv(c, a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      c[i][j] = a[i] * b[j];
    }
  }
}



/* c = a b */
function vmxv(c, a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      var s = 0;
      for ( var k = 0; k < D; k++ ) {
        s += a[i][k] * b[k];
      }
      c[i] = s;
    }
  }
}



/* c = a b */
function mmxm(c, a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      var s = 0;
      for ( var k = 0; k < D; k++ ) {
        s += a[i][k] * b[k][j];
      }
      c[i][j] = s;
    }
  }
}



/* c = a^T b */
function mmtxm(c, a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      var s = 0;
      for ( var k = 0; k < D; k++ ) {
        s += a[k][i] * b[k][j];
      }
      c[i][j] = s;
    }
  }
}



/* c = a b^T */
function mmxmt(c, a, b)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      c[i][j] = vdot(a[i], b[j]);
    }
  }
}



function mxrot3d(m, theta)
{
  theta *= Math.PI / 180;
  var c = Math.cos(theta);
  var s = Math.sin(theta);
  var m2 = newarr2d(3, 3);
  var d;
  for ( d = 0; d < 3; d++ ) {
    m2[0][d] = m[0][d];
    m2[1][d] = c * m[1][d] - s * m[2][d];
    m2[2][d] = s * m[1][d] + c * m[2][d];
  }
  return m2;
}



function myrot3d(m, theta)
{
  theta *= Math.PI / 180;
  var c = Math.cos(theta);
  var s = Math.sin(theta);
  var m2 = newarr2d(3, 3);
  var d;
  for ( d = 0; d < 3; d++ ) {
    m2[1][d] = m[1][d];
    m2[2][d] = c * m[2][d] - s * m[0][d];
    m2[0][d] = s * m[2][d] + c * m[0][d];
  }
  return m2;
}



function mzrot3d(m, theta)
{
  theta *= Math.PI / 180;
  var c = Math.cos(theta);
  var s = Math.sin(theta);
  var m2 = newarr2d(3, 3);
  var d;
  for ( d = 0; d < 3; d++ ) {
    m2[2][d] = m[2][d];
    m2[0][d] = c * m[0][d] - s * m[1][d];
    m2[1][d] = s * m[0][d] + c * m[1][d];
  }
  return m2;
}



function m2str(m)
{
  var s = "[";
  for ( var d = 0; d < D; d++ ) {
    if ( d > 0 ) {
      s += ", ";
    }
    s += m[d].toString();
  }
  s += "]";
  return s;
}



/* a += b * s */
function msinc(a, b, s)
{
  for ( var i = 0; i < D; i++ ) {
    for ( var j = 0; j < D; j++ ) {
      a[i][j] += b[i][j] * s;
    }
  }
}



/* compute the inverse matrix b = a^(-1), by Gaussian elimination */
function minv(b, a)
{
  var i, j, k, ip, x;

  // initialize b as the identity matrix
  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      b[i][j] = (i == j);
    }
  }

  // Gaussian elimination
  for ( i = 0; i < D; i++ ) {
    // choose the pivot as the largest element of column i
    x = Math.abs( a[i][i] );
    for ( ip = i, k = ip + 1; k < D; k++ ) {
      if ( Math.abs( a[k][i] ) > x ) {
        ip = k;
        x = Math.abs( a[k][i] );
      }
    }

    // swap the pivot (ip'th) row with the present row i
    for ( k = i; k < D; k++ ) {
      x = a[i][k];
      a[i][k] = a[ip][k];
      a[ip][k] = x;
    }
    for ( k = 0; k < D; k++ ) {
      x = b[i][k];
      b[i][k] = b[ip][k];
      b[ip][k] = x;
    }

    // normalize this row
    x = a[i][i];
    if ( Math.abs(x) < 1e-300 ) {
      console.log("Error: singular matrix");
      return -1;
    }
    for ( k = i; k < D; k++ ) {
      a[i][k] /= x;
    }
    for ( k = 0; k < D; k++ ) {
      b[i][k] /= x;
    }

    // use the pivot row to zero the rest rows
    for ( j = i + 1; j < D; j++ ) {
      x = a[j][i];
      for ( k = i; k < D; k++ ) {
        a[j][k] -= x * a[i][k];
      }
      for ( k = 0; k < D; k++ ) {
        b[j][k] -= x * b[i][k];
      }
    }
  }

  // now that the matrix is upper triangular, make it diagonal
  for ( i = D - 1; i >= 0; i-- ) {
    // note a[i][i] should be 1 now
    for ( j = 0; j < i; j++ ) {
      x = a[j][i];
      for ( k = 0; k < D; k++ ) {
        b[j][k] -= b[i][k] * x;
      }
    }
  }
  return 0;
}




/* full pivot
 * return the pivot row r and column c, starting from (r0, c0)
 * cmap[r0] registers the actual column index */
function mpivotf_(m, r0, cmap, sgn)
{
  var i, j, r = r0, c = r0;
  var tmp, max, t;

  // 1. find the pivot row and column
  max = -1;
  for ( j = r0; j < D; j++ ) {
    for ( i = r0; i < D; i++ ) {
      if ( (tmp = Math.abs(m[i][j])) > max ) {
        r = i;
        c = j;
        max = tmp;
      }
    }
  }

  // 2. put the pivot to the left-top corner
  // swap rows r and r0, which doesn't affect the solution
  if ( r != r0 ) {
    vswap(m[r], m[r0]);
    if ( sgn ) {
      sgn *= -1;
    }
  }

  if ( c != r0 ) { // swap columns c and r0
    if ( sgn ) {
      sgn *= -1;
    }
    for ( i = 0; i < D; i++ ) { // must be from row 0
      t = m[i][c];
      m[i][c] = m[i][r0];
      m[i][r0] = t;
    }
    if ( cmap ) {
      i = cmap[c];
      cmap[c] = cmap[r0];
      cmap[r0] = i;
    }
  }

  return sgn;
}



/* return the determinant */
function mdet(m)
{
  var y, det = 1.0, a = newarr2d(D, D);
  var i, j, k, sgn = 1;

  mcopy(a, m);
  for ( i = 0; i < D; i++ ) {
    // find the pivot, the largest element in the matrix
    sgn = mpivotf_(a, i, null, sgn);
    if ( Math.abs(a[i][i]) <= 0 ) {
      break;
    }

    det *= a[i][i];

    // normalize the row i
    y = a[i][i];
    for ( k = i; k < D; k++ ) {
      a[i][k] /= y;
    }

    // use the pivot to simplify the matrix
    for ( j = i + 1; j < D; j++ ) { // for rows j
      for ( y = a[j][i], k = i; k < D; k++ ) { // for columns k >= i
        a[j][k] -= y * a[i][k];
      }
    }
  }

  return i < D ? 0 : det * sgn;
}



var msolvezero_reltol = 1e-12;

/* Solve matrix equation a x = 0 by Gaussian elimination (full-pivot)
 * The matrix 'a' is destroyed, solutions are saved as *row* vectors in 'x'
 * return the number of solutions */
function msolvezero(a, x, reltol)
{
  var tol = 0, y;
  var i, j, k, cmap = newarr(D), sgn = 1;

  if ( !reltol ) {
    reltol = msolvezero_reltol;
  }

  for ( i = 0; i < D; i++ ) {
    cmap[i] = i;
  }

  for ( i = 0; i < D; i++ ) {
    // find the pivot, the largest element in the matrix
    sgn = mpivotf_(a, i, cmap, sgn);
    y = Math.abs(a[i][i]);
    if ( y <= tol ) { // we have D - i solutions
      break;
    }
    if ( i === 0 ) {
      tol = y * reltol;
    }

    // normalize the row i
    y = a[i][i];
    for ( k = i; k < D; k++ ) {
      a[i][k] /= y;
    }

    // use the pivot to simplify the matrix
    for ( j = 0; j < D; j++ ) { // for rows j
      if ( j != i ) {
        for ( y = a[j][i], k = i; k < D; k++ ) { // for columns k >= i
          a[j][k] -= y * a[i][k];
        }
      }
    }
  }

  // solve the D - i solutions
  for ( j = 0; j < D - i; j++ ) {
    var v = newarr(D);
    vzero( v );
    for ( k = 0; k < i; k++ ) {
      v[ cmap[k] ] = -a[k][i + j];
    }
    v[ cmap[i + j] ] = 1.0;
    vnormalize( v );
    x.push( v );
  }

  return D - i;
}



function meigvecs_low(vecs, mat, val, reltol)
{
  var m = newarr2d(D, D);

  mcopy(m, mat); // make a matrix
  for ( var d = 0; d < D; d++ ) {
    m[d][d] -= val;
  }
  return msolvezero(m, vecs, reltol);
}



/* maximal acceptable relative tolerance for solving eigenvectors */
var meig_reltol = 1e-6;

/* given an eigenvalue, return the corresponding eigenvectors
 * Note: there might be multiple eigenvectors for the eigenvalue */
function meigvecs(vecs, mat, val, reltol)
{
  var rtol, i = 0;

  if ( !reltol ) {
    reltol = meig_reltol;
  }

  // increase the tolerance, until a solution is found
  for ( rtol = 1e-14; rtol < reltol; rtol *= 10 ) {
    if ( (i = meigvecs_low(vecs, mat, val, rtol)) > 0 )
      break;
  }
  return i;
}



/* sort `s' to descending order, order `u' and `v' correspondingly */
function msort2(s, u, v)
{
  for ( var i = 0; i < D; i ++ ) {
    for ( var k = i, j = i + 1; j < D; j++ ) {
      if ( s[j] > s[k] ) {
        k = j;
      }
    }

    if ( k != i ) {
      var t = s[i]; s[i] = s[k]; s[k] = t;
      if ( u ) {
        vswap( u[i], u[k] );
      }
      if ( v ) {
        vswap( v[i], v[k] );
      }
    }
  }
}



/* compute eigenvalues of a 3x3 matrix
 * by solving a cubic equation */
function meigval(v, a)
{
  var m, p, q, pr, pr3, a00, a11, a22;

  m = (a[0][0] + a[1][1] + a[2][2])/3;
  a00 = a[0][0] - m;
  a11 = a[1][1] - m;
  a22 = a[2][2] - m;
  q = ( a00 * (a11*a22 - a[1][2]*a[2][1])
      + a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a22)
      + a[0][2] * (a[1][0]*a[2][1] - a11*a[2][0]) ) / 2.0;
  p = (a00*a00 + a11*a11 + a22*a22) / 6.0
    + (a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]) / 3.0;
  // solve x^3 - 3 p x  - 2 q = 0
  pr = Math.sqrt(p);
  pr3 = p * pr;
  if ( pr3 <= Math.abs(q) ) {
    if ( q < 0.0 ) { // choose phi = pi/3
      v[1] = v[0] = m + pr;
      v[2] = m - 2.0 * pr;
    } else { // phi = 0
      v[0] = m + 2.0 * pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    var phi = Math.acos(q / pr3)/3.0; // 0 < phi < pi/3

    v[0] = m + 2.0 * pr * Math.cos(phi);  // largest
    v[1] = m + 2.0 * pr * Math.cos(phi - 2*Math.PI/3); // second largest
    v[2] = m + 2.0 * pr * Math.cos(phi + 2*Math.PI/3); // smallest
  }
  return v;
}



/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs'
 * ideally, eigenvalues are sorted in descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
function meigsys(v, vecs, mat, nt, reltol)
{
  var vs = [];
  var n = 0, i = 0;

  if ( !reltol ) {
    reltol = meig_reltol;
  }

  meigval(v, mat);

  for ( i = 0; i < 3; i++ ) {
    n = meigvecs_low(vs, mat, v[vs.length], reltol);
    if ( n === 0 ) {
      console.log("meigsys failed: increase msolvezero_reltol", i, vs.length);
      return -1;
    }
    if ( vs.length >= 3 ) {
      break;
    }
  }

  mcopy(vecs, vs);
  msort2(v, vecs, null);

  if ( !nt ) {
    mtrans(vecs);
  }
  return 0;
}



var msvd_reltol = 1e-6;

/* SVD decomposition of a matrix A = U S V^T */
function msvd(a, u, s, v)
{
  var i, rank;
  var ata = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var us = [[0, 0, 0],[0, 0, 0], [0, 0, 0]];

  // A^T A = V S^2 V^T, so (A^T A) V = V S^2

  // 1. compute A^T A and its eigenvectors, which is V
  mmtxm(ata, a, a);
  meigsys(s, v, ata, 1);

  // 2. U^T = S^{-1} V^T A^T, and each row of U^T is an eigenvector
  // since eigenvectors are to be normalized, S^{-1} is unnecessary
  if (s[0] <= 0.0) {
    rank = 0;
    mcopy(u, v);
  } else {
    var tol = msvd_reltol;

    // the test i = 1 + (s[1] > s[0]*tol) + (s[2] > s[0]*tol);
    mmxmt(u, v, a);
    for ( i = 0; i < 3; i++ ) {
      vcopy(us[i], u[i]); // save a copy of V^T A^T before normalizing it
      s[i] = vnorm(u[i]);
      if ( s[i] > 0 ) {
        vsmul(u[i], 1/s[i]);
      }
    }
    rank = 1;
    rank += (Math.abs( vdot(u[0], u[1]) ) < tol && s[1] > tol);
    rank += (Math.abs( vdot(u[0], u[2]) ) < tol
          && Math.abs( vdot(u[1], u[2]) ) < tol && s[2] > tol);
    if ( rank <= 2 ) {
      if ( rank == 1 ) {
        var z = [0, 0, 0], w, tmp;

        w = Math.abs( u[0][i = 0] );
        if ( (tmp = Math.abs(u[0][1])) < w ) {
          w = tmp;
          i = 1;
        }
        if ( (tmp = Math.abs(u[0][2])) < w ) {
          i = 2;
        }
        z[i] = 1.0; // select the smallest element in u[0] as z
        vnormalize( vcross3d(u[1], z, u[0]) );
        s[1] = vdot(u[1], us[1]); // S = U^T (V^T A^T)^T is more accurate than sqrt(A^T A)
        if (s[1] < 0) {
          s[1] = -s[1];
          vneg( u[1] );
        } // make sure s[1] > 0
      }
      vnormalize( vcross3d(u[2], u[0], u[1]) );
      s[2] = vdot(u[2], us[2]);
      if ( s[2] < 0 ) {
        s[2] = -s[2];
        vneg( u[2] );
      }
    }
    msort2(s, u, v);
  }
  mtrans(v);
  mtrans(u);
}



/* Fit x to y by rotation and translation of the `x'
 * If `refl', reflection can also be used.
 * The best-fit structure is saved to `xf', if not null */
function vrmsd(x, xf, y, w, n, refl, r, t)
{
  var i;
  var wi, wtot = 0, sq, dev = 0, dev0, detm;
  var xc = [0, 0, 0], yc = [0, 0, 0], xs = [0, 0, 0], ys = [0, 0, 0];
  var sig = [0, 0, 0], t_ = [0, 0, 0];
  var u = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var v = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var s = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xy = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var r_ = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xfi = [0, 0, 0], dx = [0, 0, 0];


  if ( !r ) {
    r = r_;
  }
  if ( !t ) {
    t = t_;
  }

  // 1. compute the centers
  vzero( xc );
  vzero( yc );
  if ( !w ) {
    for ( i = 0; i < n; i++ ) {
      vinc(xc, x[i]);
      vinc(yc, y[i]);
    }
    wtot = n;
  } else {
    for ( wtot = 0.0, i = 0; i < n; i++ ) {
      vsinc(xc, x[i], w[i]);
      vsinc(yc, y[i], w[i]);
      wtot += w[i];
    }
  }
  vsmul(xc, 1.0/wtot);
  vsmul(yc, 1.0/wtot);

  // 2. compute the asymmetric covariance matrix S = (x-xc) (y-yc)^T
  for ( i = 0; i < n; i++ ) {
    vdiff(xs, x[i], xc); // shift to the center avoid the translation
    vdiff(ys, y[i], yc);
    mvtxv(xy, xs, ys);
    sq  = vsqr(xs);
    sq += vsqr(ys);
    if ( w ) {
      wi = w[i];
    } else {
      wi = 1.0;
    }
    msinc(s, xy, wi);
    dev += wi * sq; // Tr(x^T x + y^T y)
  }
  dev0 = dev;

  // 3. SVD decompose S = u sig v^T
  msvd(s, u, sig, v);

  // 4. compute R = v u^T
  mmxmt(r, v, u);
  detm = mdet(r);

  if ( detm < 0 && !refl ) { // to avoid a reflection
    mtrans(u);
    vneg(u[2]); // flip the last eigenvector
    mmxm(r, v, u);
    dev -= 2*(sig[0] + sig[1] - sig[2]);
    detm = mdet(r);
  } else {
    dev -= 2 * (sig[0] + sig[1] + sig[2]); // -2 Tr(R x y^T)
  }
  if ( dev < 0 ) {
    dev = 0;
  }
  vmxv(xs, r, xc); // xs = R xc
  vdiff(t, yc, xs); // t = yc - R xc

  // 5. compute the rotated structure
  if ( xf || dev < dev0 * 0.01 ) { // if there's a large cancellation recompute the deviation
    for ( dev = 0, i = 0; i < n; i++ ) {
      vmxv(xs, r, x[i]); // xs = R x
      vadd(xfi, xs, t); // xfi = R x + t
      sq = vsqr( vdiff(dx, y[i], xfi) );
      if ( xf ) {
        vcopy(xf[i], xfi);
      }
      dev += (w ? w[i]*sq : sq); // recompute the deviation
    }
  }
  return Math.sqrt(dev/wtot);
}




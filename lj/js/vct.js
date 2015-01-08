/* vector routines */



"use strict";



var D = 2;



function vzero(x)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] = 0;
  }
  return x;
}



function vcopy(x, y)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] = y[d];
  }
  return x;
}



function vinc(x, dx)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] += dx[d];
  }
  return x;
}



function vdec(x, dx)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] -= dx[d];
  }
  return x;
}



function vsinc(x, dx, s)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] += dx[d] * s;
  }
  return x;
}



function vdiff(c, a, b)
{
  for ( var d = 0; d < D; d++ ) {
    c[d] = a[d] - b[d];
  }
  return c;
}



function vsmul(x, s)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] *= s;
  }
  return x;
}



function vsqr(x)
{
  for ( var s = 0, d = 0; d < D; d++ ) {
    s += x[d] * x[d];
  }
  return s;
}



function vdot(x, y)
{
  for ( var s = 0, d = 0; d < D; d++ ) {
    s += x[d] * y[d];
  }
  return s;
}



function vcross2d(x, y)
{
  return x[0]*y[1] - x[1]*y[0];
}



function vcross3d(z, x, y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}



/* inverse matrix a^(-1) */
function rm3_inv(a)
{
  var d00 = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  var d01 = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  var d02 = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  var detm = a[0][0]*d00 + a[0][1]*d01 + a[0][2]*d02;
  var dmin = 1e-20;

  if (detm < dmin && detm > -dmin) {
    detm = (detm < 0) ? -dmin: dmin;
  }
  var b = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  b[0][0] = d00/detm;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detm;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/detm;
  b[1][0] = d01/detm;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/detm;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/detm;
  b[2][0] = d02/detm;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/detm;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detm;
  return b;
}



function vpbc(v, l, invl)
{
  for ( var d = 0; d < D; d++ ) {
    v[d] -= (Math.floor(v[d]*invl + 1000.5) - 1000.0) * l;
  }
  return v;
}



function vwrap(v, l)
{
  for ( var d = 0; d < D; d++ ) {
    while ( v[d] < 0 ) {
      v[d] += l;
    }
    while ( v[d] > l ) {
      v[d] -= l;
    }
  }
}



/* return a unit matrix */
function munit(a)
{
  for ( var d1 = 0; d1 < D; d1++ ) {
    for ( var d2 = 0; d2 < D; d2++ ) {
      a[d1][d2] = (d1 === d2) ? 1.0 : 0.0;
    }
  }
}



/* copy a matrix, a = b*/
function mcopy(a, b)
{
  for ( var d1 = 0; d1 < D; d1++ ) {
    for ( var d2 = 0; d2 < D; d2++ ) {
      a[d1][d2] = b[d1][d2];
    }
  }
}



/* return the product of matrix m and vector v */
function mmulv(m, v)
{
  var u = newarr2d(D);

  for ( var d = 0; d < D; d++ ) {
    u[d] = 0;
    for ( var d2 = 0; d2 < D; d2++ ) {
      u[d] += m[d][d2] * v[d2];
    }
  }
  return u;
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


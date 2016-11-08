/* vector routines */



"use strict";



var D = 3;



function vzero(x)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] = 0.0;
  }
  return x;
}



function vneg(x)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] = -x[d];
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



function vswap(x, y)
{
  for ( var d = 0; d < D; d++ ) {
    var z = x[d];
    x[d] = y[d];
    y[d] = z;
  }
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



function vadd(c, a, b)
{
  for ( var d = 0; d < D; d++ ) {
    c[d] = a[d] + b[d];
  }
  return c;
}



function vsadd(c, a, b, s)
{
  for ( var d = 0; d < D; d++ ) {
    c[d] = a[d] + s * b[d];
  }
  return c;
}



function vdiff(c, a, b)
{
  for ( var d = 0; d < D; d++ ) {
    c[d] = a[d] - b[d];
  }
  return c;
}



function vnadd(c, a, b)
{
  for ( var d = 0; d < D; d++ ) {
    c[d] = - a[d] - b[d];
  }
  return c;
}



function vlincomb2(z, x, y, c, s)
{
  for ( var d = 0; d < D; d++ ) {
    z[d] = x[d] * c + y[d] * s;
  }
  return z;
}



function vsmul(x, s)
{
  for ( var d = 0; d < D; d++ ) {
    x[d] *= s;
  }
  return x;
}



function vsmul2(y, x, s)
{
  for ( var d = 0; d < D; d++ ) {
    y[d] = x[d] * s;
  }
  return y;
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


// bond angle interaction
function vang(xi, xj, xk, gi, gj, gk)
{
  var xij= [0,0,0], xkj = [0,0,0], ri, rk, dot, ang;

  ri = Math.sqrt( vsqr( vdiff(xij, xi, xj) ) );
  vsmul(xij, 1.0/ri);

  rk = Math.sqrt( vsqr( vdiff(xkj, xk, xj) ) );
  vsmul(xkj, 1.0/rk);

  dot = Math.max( Math.min( vdot(xij, xkj), 1.0 ), -1.0);
  ang = Math.acos( dot );

  if ( gi && gj && gk ) {
    var sn, gij, gkj;
    var d;
    sn = -1.0 / Math.sqrt(1 - dot * dot); // -1.0/sin(phi)
    for ( d = 0; d < D; d++ ) {
      gij = sn * (xkj[d] - xij[d]*dot) / ri;
      gkj = sn * (xij[d] - xkj[d]*dot) / rk;
      gi[d] = gij;
      gk[d] = gkj;
      gj[d] = -(gij + gkj);
    }
  }
  return ang;
}



function vdih(xi, xj, xk, xl, gi, gj, gk, gl)
{
  var tol, phi, cosphi = 1;
  var nxkj, nxkj2, m2, n2;
  var xij = [0,0,0], xkj = [0,0,0], xkl = [0,0,0];
  var uvec = [0,0,0], vvec = [0,0,0], svec = [0,0,0];
  var m = [0,0,0], n = [0,0,0]; // the planar vector of xij x xkj, and xkj x xkj

  vdiff(xij, xi, xj);
  vdiff(xkj, xk, xj);
  vdiff(xkl, xk, xl);
  nxkj2 = vsqr(xkj);
  nxkj = Math.sqrt(nxkj2);
  tol = nxkj2 * 1e-16;

  m2 = vsqr( vcross3d(m, xij, xkj) );
  n2 = vsqr( vcross3d(n, xkj, xkl) );
  if (m2 > tol && n2 > tol) {
    cosphi = Math.max( Math.min(
          vdot(m, n) / Math.sqrt(m2 * n2),
          1), -1);
  }
  phi = Math.acos(cosphi);
  if (vdot(n, xij) < 0.0) {
    phi = -phi;
  }

  // optionally calculate the gradient
  if ( gi && gj && gk && gl ) {
    if ( m2 > tol && n2 > tol ) {
      vsmul2(gi, m, nxkj/m2);
      vsmul2(gl, n, -nxkj/n2);
      vsmul2(uvec, gi, vdot(xij, xkj)/nxkj2);
      vsmul2(vvec, gl, vdot(xkl, xkj)/nxkj2);
      vdiff(svec, uvec, vvec);
      vdiff(gj, svec, gi);
      vnadd(gk, svec, gl);
    } else { /* clear the gradients */
      vzero(gi);
      vzero(gj);
      vzero(gk);
      vzero(gl);
    }
  }
  return phi;
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



/* return the norm the vector */
function vnorm(a)
{
  return Math.sqrt( vsqr(a) );
}



/* return the distance */
function vdistx(dx, a, b)
{
  return vnorm( vdiff(dx, a, b) );
}



/* return the distance */
function vdist(a, b)
{
  var dx = newarr(D);
  return vnorm( vdiff(dx, a, b) );
}



/* normalize the vector */
function vnormalize(v)
{
  var s = Math.sqrt( vsqr(v) );
  return vsmul(v, 1.0 / s);
}



/* perpendicular component of `x' to `y', normalized */
function vperpen(z, x, y)
{
  return vnormalize( vsadd(z, x, y, -vdot(x, y)/vsqr(y)) );
}



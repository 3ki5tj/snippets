



"use strict";



/* bond length interaction */
function ebondlen(r0, k, xi, xj, fi, fj)
{
  var dx = newarr(D);
  vdiff(dx, xi, xj);
  var r2 = vsqr(dx);
  var r = Math.sqrt(r2);
  var dr = r - r0;
  if ( fi && fj ) {
    var amp = 2 * k * dr / r;
    vsinc(fi, dx, -amp);
    vsinc(fj, dx,  amp);
  }
  return k * dr * dr;
}



/* bond angle interaction */
function ebondang(ang0, k, xi, xj, xk, fi, fj, fk)
{
  var xij = newarr(D), xkj = newarr(D), ri, rk;
  var dot, ang, dang;

  // compute the cosine of the angle
  ri = Math.sqrt( vsqr( vdiff(xij, xi, xj) ) );
  vsmul(xij, 1.0/ri); // unit vector from j to i

  rk = Math.sqrt( vsqr( vdiff(xkj, xk, xj) ) );
  vsmul(xkj, 1.0/rk); // unit vector from j to k

  dot = Math.max( Math.min( vdot(xij, xkj), 1.0 ), -1.0 );
  ang = Math.acos(dot);
  dang = ang - ang0;
  //console.log(ang, dang, dot, xi, xj, xk, xij, xkj, ri, rk);

  if ( fi && fj && fk ) {
    var sn = -1.0 / Math.sqrt( 1 - dot * dot );
    var d, amp = 2 * k * dang, gi, gk;
    for ( d = 0; d < D; d++ ) {
      gi = sn * (xkj[d] - xij[d]*dot) / ri;
      gk = sn * (xij[d] - xkj[d]*dot) / rk;
      fi[d] -= amp * gi;
      fk[d] -= amp * gk;
      fj[d] += amp * (gi + gk);
    }
  }
  return k * dang * dang;
}



/* WCA interaction */
function ewca(sig2, eps, xi, xj, fi, fj)
{
  var dx = newarr(D), dr2, ir2, ir6, fs;

  dr2 = vsqr( vdiff(dx, xi, xj) );
  if ( dr2 < sig2 ) {
    ir2 = sig2 / dr2;
    ir6 = ir2 * ir2 * ir2;
    if ( fi && fj ) {
      fs = eps * ir6 * (12 * ir6 - 12); // f.r
      fs /= dr2; // f.r / r^2
      vsinc(fi, dx, fs);
      vsinc(fj, dx, -fs);
    }
    return eps * (ir6 * (ir6 - 2) + 1);
  } else {
    return 0;
  }
}





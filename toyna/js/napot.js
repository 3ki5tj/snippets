



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



/* stack interaction */
function estack(r0, phi10, phi20, kr, kphi, ust0,
    xp1, xs1, xb1, xp2, xs2, xb2, xp3,
    fp1, fs1, fb1, fp2, fs2, fb2, fp3)
{
  var dxbb = [0,0,0], rbb, drbb = 0;
  var phi1 = 0, dphi1 = 0, gi1 = [0,0,0], gj1 = [0,0,0], gk1 = [0,0,0], gl1 = [0,0,0];
  var phi2 = 0, dphi2 = 0, gi2 = [0,0,0], gj2 = [0,0,0], gk2 = [0,0,0], gl2 = [0,0,0];
  var den;

  rbb = Math.sqrt( vsqr( vdiff(dxbb, xb1, xb2) ) );
  drbb = rbb - r0;

  phi1 = vdih(xp1, xs1, xp2, xs2, gi1, gj1, gk1, gl1);
  dphi1 = phi1 - phi10;

  if ( xp3 ) {
    phi2 = vdih(xs1, xp2, xs2, xp3, gi2, gj2, gk2, gl2);
    dphi2 = phi2 - phi20;
  }

  den = 1 + kr * drbb * drbb + kphi * dphi1 * dphi1 + kphi * dphi2 * dphi2;

  if ( fp1 && fs1 && fb1 && fp2 && fs2 && fb2 ) {
    var fs, den2 = den * den;
    
    fs = 2 * kr * drbb * ust0 / den2 / rbb;
    vsinc(fb1, dxbb,  fs);
    vsinc(fb2, dxbb, -fs);

    fs = 2 * kphi * dphi1 * ust0 / den2;
    vsinc(fp1, gi1, fs);
    vsinc(fs1, gj1, fs);
    vsinc(fp2, gk1, fs);
    vsinc(fs2, gl1, fs);

    if ( xp3 ) {
      fs = 2 * kphi * dphi2 * ust0 / den2;
      vsinc(fs1, gi2, fs);
      vsinc(fp2, gj2, fs);
      vsinc(fs2, gk2, fs);
      vsinc(fp3, gl2, fs);
    }

  }
  return ust0 / den;
}



/* dielectric constant of water
 * Eq. (12) of Denesyuk 2013 */
function getdielecwater(tp)
{
  tp -= 273.15;
  return 87.740 - 0.4008*tp + 9.398e-4*tp*tp - 1.410e-6*tp*tp*tp;
}



/* return the Bjerrum length
 * Eq. (11) of Denesyuk 2013 */
function getBjerrumlen(tp)
{
  var eps = getdielecwater(tp);
  return KE2 / (eps * BOLTZK * tp);
}



/* return the charge reduction factor Q
 * Eq. (10) of Denesyuk 2013 */
function getchargeQ(tp)
{
  var b = 4.4;
  return b / getBjerrumlen(tp);
}



/* electrostatic interaction under Debye-Huckel approximation */
function echargeDH(QQ ,debyel, xi, xj, fi, fj)
{
  var dx = [0,0,0], dr, xp, fs;

  dr = Math.sqrt( vsqr( vdiff(dx, xi, xj) ) );
  xp = Math.exp(-dr/debyel);
  if ( fi && fj ) {
    fs = xp * (1/debyel + 1/dr) / (dr * dr);
    vsinc(fi, dx, fs);
    vsinc(fj, dx, -fs);
  }
  return QQ * xp / dr;
}



/* get the Debye screening length */
function getDebyel(q, conc, n, tp)
{
  var i;
  var s = 0;

  for ( i = 0; i < n; i++ ) {
    s += q[i] * q[i] * conc[i] * (AVOGADRO * 1e-27);
  }
  return Math.sqrt( getdielecwater(tp) * BOLTZK * tp / (s * 4 * PI * KE2) );
}





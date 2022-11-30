#ifndef NAPOT_H__
#define NAPOT_H__



#include "vct.h"



/* bond length interaction */
__inline static double ebondlen(double r0, double k,
    const double *xi, const double *xj,
    double *fi, double *fj)
{
  double dx[D], r, dr, amp;

  r = sqrt( vsqr( vdiff(dx, xi, xj) ) );
  dr = r - r0;
  if ( fi && fj ) {
    amp = 2 * k * dr / r;
    vsinc(fi, dx, -amp);
    vsinc(fj, dx,  amp);
  }
  return k * dr * dr;
}



/* bond angle interaction */
__inline static double ebondang(double ang0, double k,
    const double *xi, const double *xj, const double *xk,
    double *fi, double *fj, double *fk)
{
  double xij[D], xkj[D], ri, rk, dot, ang, dang;

  ri = sqrt( vsqr( vdiff(xij, xi, xj) ) );
  vsmul(xij, 1.0/ri);

  rk = sqrt( vsqr( vdiff(xkj, xk, xj) ) );
  vsmul(xkj, 1.0/rk);

  dot = vdot(xij, xkj);
  if ( dot > 1.0 ) dot = 1.0;
  else if ( dot < -1.0 ) dot = -1.0;
  ang = acos( dot );
  dang = ang - ang0;

  if ( fi && fj && fk ) {
    double sn, amp, gi, gk;
    int d;
    sn = -1.0 / sqrt(1 - dot * dot); /* -1.0/sin(phi) */
    amp = 2 * k * dang * sn;
    for ( d = 0; d < D; d++ ) {
      gi = (xkj[d] - xij[d]*dot) / ri;
      gk = (xij[d] - xkj[d]*dot) / rk;
      fi[d] -= amp * gi;
      fk[d] -= amp * gk;
      fj[d] += amp * (gi + gk);
    }
  }
  return k * dang * dang;
}



/* WCA interaction */
__inline static double ewca(double sig2, double eps,
    const double *xi, const double *xj, double *fi, double *fj)
{
  double dx[D], dr2, ir2, ir6, fs;

  dr2 = vsqr( vdiff(dx, xi, xj) );
  if ( dr2 < sig2 ) {
    ir2 = sig2 / dr2;
    ir6 = ir2 * ir2 * ir2;
    if ( fi && fj ) {
      fs = eps * ir6 * (12 * ir6 - 12); /* f.r */
      fs /= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(fj, dx, -fs);
    }
    return eps * (ir6 * (ir6 - 2) + 1);
  } else {
    return 0;
  }
}



/* stack interaction */
__inline static double estack(double r0, double phi10, double phi20,
    double kr, double kphi, double ust0,
    const double *xp1, const double *xs1, const double *xb1,
    const double *xp2, const double *xs2, const double *xb2,
    const double *xp3,
    double *fp1, double *fs1, double *fb1,
    double *fp2, double *fs2, double *fb2,
    double *fp3)
{
  double dxbb[D], rbb, drbb = 0;
  double phi1 = 0, dphi1 = 0, gi1[D], gj1[D], gk1[D], gl1[D];
  double phi2 = 0, dphi2 = 0, gi2[D] = {0}, gj2[D] = {0}, gk2[D] = {0}, gl2[D] = {0};
  double den;

  rbb = sqrt( vsqr( vdiff(dxbb, xb1, xb2) ) );
  drbb = rbb - r0;

  phi1 = vdih(xp1, xs1, xp2, xs2, gi1, gj1, gk1, gl1);
  dphi1 = phi1 - phi10;

  if ( xp3 ) {
    phi2 = vdih(xs1, xp2, xs2, xp3, gi2, gj2, gk2, gl2);
    dphi2 = phi2 - phi20;
  }

  den = 1 + kr * drbb * drbb + kphi * dphi1 * dphi1 + kphi * dphi2 * dphi2;

  if ( fp1 && fs1 && fb1 && fp2 && fs2 && fb2 ) {
    double fs, den2 = den * den;

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



/* hydrogen-bond interaction */
__inline static double ehbond(double r0, double th10, double th20,
    double psi0, double psi10, double psi20,
    double kr, double kth, double kpsi, double uhb0,
    const double *x1, const double *x2, const double *x3,
    const double *x4, const double *x5, const double *x6,
    double *f1, double *f2, double *f3,
    double *f4, double *f5, double *f6)
{
  double dx[D], r, dr = 0;
  double th1  = 0, dth1  = 0, ga1[D], gb1[D], gc1[D];
  double th2  = 0, dth2  = 0, ga2[D], gb2[D], gc2[D];
  double psi  = 0, dpsi  = 0, gi [D], gj [D], gk [D], gl [D];
  double psi1 = 0, dpsi1 = 0, gi1[D], gj1[D], gk1[D], gl1[D];
  double psi2 = 0, dpsi2 = 0, gi2[D], gj2[D], gk2[D], gl2[D];
  double den;

  r = sqrt( vsqr( vdiff(dx, x3, x4) ) );
  dr = r - r0;

  th1 = vang(x2, x3, x4, ga1, gb1, gc1);
  dth1 = th1 - th10;

  th2 = vang(x3, x4, x5, ga2, gb2, gc2);
  dth2 = th2 - th20;

  psi = vdih(x2, x3, x4, x5, gi, gj, gk, gl);
  dpsi = psi - psi0;

  psi1 = vdih(x1, x2, x3, x4, gi1, gj1, gk1, gl1);
  dpsi1 = psi1 - psi10;

  psi2 = vdih(x3, x4, x5, x6, gi2, gj2, gk2, gl2);
  dpsi2 = psi2 - psi20;

  den = 1 + kr * dr * dr + kth * (dth1 * dth1 + dth2 * dth2)
    + kpsi * (dpsi * dpsi + dpsi1 * dpsi1 + dpsi2 * dpsi2);

  if ( f1 && f2 && f3 && f4 && f5 && f6 ) {
    double fs, den2 = den * den;

    fs = 2 * kr * dr * uhb0 / den2 / r;
    vsinc(f3, dx,  fs);
    vsinc(f4, dx, -fs);

    fs = 2 * kth * dth1 * uhb0 / den2;
    vsinc(f2, ga1, fs);
    vsinc(f3, gb1, fs);
    vsinc(f4, gc1, fs);

    fs = 2 * kth * dth2 * uhb0 / den2;
    vsinc(f3, ga2, fs);
    vsinc(f4, gb2, fs);
    vsinc(f5, gc2, fs);

    fs = 2 * kpsi * dpsi * uhb0 / den2;
    vsinc(f2, gi, fs);
    vsinc(f3, gj, fs);
    vsinc(f4, gk, fs);
    vsinc(f5, gl, fs);

    fs = 2 * kpsi * dpsi1 * uhb0 / den2;
    vsinc(f1, gi1, fs);
    vsinc(f2, gj1, fs);
    vsinc(f3, gk1, fs);
    vsinc(f4, gl1, fs);

    fs = 2 * kpsi * dpsi2 * uhb0 / den2;
    vsinc(f3, gi2, fs);
    vsinc(f4, gj2, fs);
    vsinc(f5, gk2, fs);
    vsinc(f6, gl2, fs);
  }
  return uhb0 / den;
}



/* dielectric constant of water
 * Eq. (12) of Denesyuk 2013 */
__inline static double getdielecwater(double tp)
{
  tp -= 273.15;
  return 87.740 - 0.4008*tp + 9.398e-4*tp*tp - 1.410e-6*tp*tp*tp;
}



/* return the Bjerrum length
 * Eq. (11) of Denesyuk 2013 */
__inline static double getBjerrumlen(double tp, double *eps)
{
  *eps = getdielecwater(tp);
  return KE2 / (*eps * BOLTZK * tp);
}



/* return the charge reduction factor Q
 * Eq. (10) of Denesyuk 2013 */
__inline static double getchargeQ(double tp, double *eps)
{
  const double b = 4.4;
  return b / getBjerrumlen(tp, eps);
}



/* electrostatic interaction under Debye-Huckel approximation */
__inline static double echargeDH(double QQ, double debyel,
    const double *xi, const double *xj, double *fi, double *fj)
{
  double dx[D], dr, xp, fs;

  dr = sqrt( vsqr( vdiff(dx, xi, xj) ) );
  xp = exp(-dr/debyel);
  if ( fi && fj ) {
    fs = xp * (1/debyel + 1/dr) / (dr * dr);
    vsinc(fi, dx, fs);
    vsinc(fj, dx, -fs);
  }
  return QQ * xp / dr;
}



/* get the Debye screening length */
__inline static double getDebyel(const double *q, const double *conc,
    int n, double tp)
{
  int i;
  double s = 0;

  for ( i = 0; i < n; i++ ) {
    s += q[i] * q[i] * conc[i] * (AVOGADRO * 1e-27);
  }
  return sqrt( getdielecwater(tp) * BOLTZK * tp / (s * 4 * PI * KE2) );
}



#endif /* POTENTIAL_H__ */

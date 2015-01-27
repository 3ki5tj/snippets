#ifndef POTENTIAL_H__
#define POTENTIAL_H__



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



#endif /* POTENTIAL_H__ */

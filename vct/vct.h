#ifndef VCT_H__
#define VCT_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>



#ifndef PI
#define PI 3.141592653589793
#endif



typedef double vct[D];



__inline static void vzero(double *x)
{
  int d;
  for ( d = 0; d < D; d++ ) x[d] = 0;
}



__inline static void vcopy(double *x, const double *y)
{
  int d;
  for ( d = 0; d < D; d++ ) x[d] = y[d];
}



#define vinc(x, dx) vsinc(x, dx, 1)
#define vdec(x, dx) vsinc(x, dx, -1)

__inline static double *vsinc(double *x, const double *dx, double s)
{
  int d;
  for ( d = 0; d < D; d++ ) x[d] += dx[d] * s;
  return x;
}



__inline static double *vadd(double *c, const double *a, const double *b)
{
  int d;
  for ( d = 0; d < D; d++ ) c[d] = a[d] + b[d];
  return c;
}



__inline static double *vdiff(double *c, const double *a, const double *b)
{
  int d;
  for ( d = 0; d < D; d++ ) c[d] = a[d] - b[d];
  return c;
}



__inline static double *vnadd(double *c, const double *a, const double *b)
{
  int d;
  for ( d = 0; d < D; d++ ) c[d] = -a[d] - b[d];
  return c;
}



__inline static double *vsmul(double *x, double s)
{
  int d;
  for ( d = 0; d < D; d++ ) x[d] *= s;
  return x;
}



__inline static double *vsmul2(double *y, const double *x, double s)
{
  int d;
  for ( d = 0; d < D; d++ ) y[d] = x[d] * s;
  return y;
}



#define vsqr(x) vdot(x, x)

__inline static double vdot(const double *x, const double *y)
{
  int d;
  double s = 0;
  for ( d = 0; d < D; d++ ) s += x[d] * y[d];
  return s;
}



#if D == 2


__inline static double vcross(double *x, double *y)
{
  return x[0]*y[1] - x[1]*y[0];
}



#elif D == 3



__inline static double *vcross(double *z, const double *x, const double *y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}



/* inverse matrix b = a^(-1) */
__inline static void rm3_inv(double b[3][3], double a[3][3])
{
  double d00 = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  double d01 = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  double d02 = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  double detm = a[0][0]*d00 + a[0][1]*d01 + a[0][2]*d02;
  const double dmin = 1e-20;

  if (detm < dmin && detm > -dmin)
    detm = (detm < 0) ? -dmin: dmin;
  b[0][0] = d00/detm;
  b[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detm;
  b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/detm;
  b[1][0] = d01/detm;
  b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2])/detm;
  b[1][2] = (a[0][2]*a[1][0] - a[1][2]*a[0][0])/detm;
  b[2][0] = d02/detm;
  b[2][1] = (a[2][0]*a[0][1] - a[2][1]*a[0][0])/detm;
  b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detm;
}



/* bond angle interaction */
__inline static double vang(const double *xi, const double *xj, const double *xk,
    double *gi, double *gj, double *gk)
{
  double xij[D], xkj[D], ri, rk, dot, ang;

  ri = sqrt( vsqr( vdiff(xij, xi, xj) ) );
  vsmul(xij, 1.0/ri);

  rk = sqrt( vsqr( vdiff(xkj, xk, xj) ) );
  vsmul(xkj, 1.0/rk);

  dot = vdot(xij, xkj);
  if ( dot > 1.0 ) dot = 1.0;
  else if ( dot < -1.0 ) dot = -1.0;
  ang = acos( dot );

  if ( gi && gj && gk ) {
    double sn, gij, gkj;
    int d;
    sn = -1.0 / sqrt(1 - dot * dot); /* -1.0/sin(phi) */
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



__inline static double vdih(const double *xi, const double *xj,
    const double *xk, const double *xl,
    double *gi, double *gj, double *gk, double *gl)
{
  double tol, phi, cosphi = 1;
  double nxkj, nxkj2, m2, n2;
  double xij[3], xkj[3], xkl[3], uvec[3], vvec[3], svec[3];
  double m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */

  vdiff(xij, xi, xj);
  vdiff(xkj, xk, xj);
  vdiff(xkl, xk, xl);
  nxkj2 = vsqr(xkj);
  nxkj = sqrt(nxkj2);
  tol = nxkj2 * DBL_EPSILON;

  vcross(m, xij, xkj);
  m2 = vsqr(m);
  vcross(n, xkj, xkl);
  n2 = vsqr(n);
  if (m2 > tol && n2 > tol) {
    cosphi = vdot(m, n);
    cosphi /= sqrt(m2 * n2);
    if (cosphi >= 1) cosphi = 1;
    else if (cosphi < -1) cosphi = -1;
  }
  phi = acos(cosphi);
  if (vdot(n, xij) < 0.0) phi = -phi;

  /* optionally calculate the gradient */
  if (gi != NULL) {
    if (m2 > tol && n2 > tol) {
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



#endif /* D == 3 */



#endif /* VCT_H__ */


#ifndef UTIL_H__
#define UTIL_H__



#include "mat.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



/* return the larger of a and b */
__inline double dblmax(double a, double b)
{
  return a > b ? a : b;
}



/* return the smaller of a and b */
__inline double dblmin(double a, double b)
{
  return a < b ? a : b;
}



/* remove the center of mass motion */
static void rmcom(double (*x)[D], const double *m, int n)
{
  int i;
  double xc[D] = {0}, mtot = 0, wt;

  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



#if D == 2


/* annihilate the total angular momentum */
static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double am, r2, xc[D] = {0, 0}, xi[D];
  double mtot = 0, wt;

  /* determine the center of mass */
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

  am = r2 = 0.0;
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    am += wt * vcross(xi, v[i]);
    r2 += wt * vsqr(x[i]);
  }

  am = -am / r2;
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am * xi[1];
    v[i][1] +=  am * xi[0];
  }
}



#else



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double xc[D] = {0, 0, 0}, xi[D], ang[D], am[D] = {0, 0, 0};
  double dv[D], mat[D][D], inv[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  double mtot = 0, wt;

  /* determine the center of mass */
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    vcross(ang, xi, v[i]);
    vsinc(am, ang, wt);

    xx += wt * xi[0] * xi[0];
    yy += wt * xi[1] * xi[1];
    zz += wt * xi[2] * xi[2];
    xy += wt * xi[0] * xi[1];
    yz += wt * xi[1] * xi[2];
    zx += wt * xi[2] * xi[0];
  }
  mat[0][0] = yy + zz;
  mat[1][1] = xx + zz;
  mat[2][2] = xx + yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  minv(inv, mat);

  /* ang is the solution of M^(-1) * L */
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vcross(dv, ang, xi);
    vinc(v[i], dv);
  }
}



#endif /* D == 3 */




#endif


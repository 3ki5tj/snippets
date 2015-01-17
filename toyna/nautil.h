#ifndef UTIL_H__
#define UTIL_H__



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



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



typedef double vec[D];



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



__inline static double *vdiff(double *c, const double *a, const double *b)
{
  int d;
  for ( d = 0; d < D; d++ ) c[d] = a[d] - b[d];
  return c;
}



__inline static double *vsmul(double *x, double s)
{
  int d;
  for ( d = 0; d < D; d++ ) x[d] *= s;
  return x;
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



__inline static double vcross(const double *x, const double *y)
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


#endif /* D == 3 */



#endif


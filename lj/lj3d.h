#ifndef LJ3D_H__
#define LJ3D_H__



/* Three-dimensional Lennard-Jones fluid */



#define D 3
#include "ljcore.h"



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  double a;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l/n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        lj->x[id][0] = (i + .5) * a;
        lj->x[id][1] = (j + .5) * a;
        lj->x[id][2] = (k + .5) * a;
        id++;
      }
}



/* get the tail correction */
static double lj_gettail(lj_t *lj, double rho, int n, double *ptail)
{
  double irc, irc3, irc6, utail;

  irc = 1/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = ( 8*M_PI*rho*n/9*(irc6 - 3)*irc3 );
  if (ptail != NULL)
    *ptail = ( 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3 );
  return utail;
}



/* annihilate the total angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
static void lj_shiftang(double (*x)[D], double (*v)[D], int n)
{
  int i;
  double xc[D] = {0, 0, 0}, xi[D], ang[D], am[D] = {0, 0, 0}, dv[D], mat[D][D], inv[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) vinc(xc, x[i]);
  vsmul(xc, 1.f/n);
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(ang, xi, v[i]);
    vinc(am, ang);
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  rm3_inv(inv, mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  /* ang is the solution of M^(-1) * I */
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(dv, ang, xi);
    vinc(v[i], dv);
  }
}



#endif /* LJ_H__ */

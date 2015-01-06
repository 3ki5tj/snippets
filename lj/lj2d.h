#ifndef LJ2D_H__
#define LJ2D_H__



/* Two dimensional Lennard-Jons fluid */



#define D 2
#include "ljcore.h"



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++) {
      if ((i+j) % 2 != 0) continue;
      /* add some noise to prevent two atoms happened to
       * be separated precisely by the cutoff distance,
       * which might be half of the box */
      lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
      lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
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
  utail = M_PI*rho*n*(.4*irc6 - 1)*irc3*irc;
  if (ptail != NULL)
    *ptail = M_PI*rho*rho*(2.4*irc6 - 3)*irc3*irc;
  return utail;
}



/* annihilate the total angular momentum */
static void lj_shiftang(double (*x)[D], double (*v)[D], int n)
{
  int i;
  double am, r2, xc[D] = {0, 0}, xi[D];

  for (i = 0; i < n; i++) vinc(xc, x[i]);
  vsmul(xc, 1.f/n);
  for (am = r2 = 0.f, i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    am += vcross(xi, v[i]);
    r2 += vsqr(x[i]);
  }
  am = -am/r2;
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am*xi[1];
    v[i][1] +=  am*xi[0];
  }
}



#endif /* LJ2D_H__ */

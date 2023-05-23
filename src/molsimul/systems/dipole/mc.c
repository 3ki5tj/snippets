/* Monte Carlo for a dipole system */
#include "dipole.h"

int n = 16;
double l = 2.0;
double disp = 0.1;
double sig = 0.9;
double eps = 1.0;
double amp = 0.1;
double tp = 3;
double alpha = 10.0;
int km = 6;
long nsteps = 10000;

int main(void)
{
  dipole_t *dp;
  long t, nacc = 0;

  dp = dipole_open(n, l, disp, sig, eps, alpha, km);
  dp->ene = dipole_force(dp, dp->x, dp->f, &dp->ene_lj, &dp->ene_el);
  //dipole_testforce(dp);
  for ( t = 0; t < nsteps; t++ ) {
    nacc += dipole_metro(dp, amp, 1/tp);
    //double ene = dipole_force(dp, dp->x, NULL, NULL, NULL);
    //printf("t %ld, %g, %g\n", t, ene, dp->ene);
  }
  printf("acc %g\n", 1.*nacc/nsteps);
  dipole_savepos(dp, "dp.pos");
  dipole_close(dp);
  return 0;
}

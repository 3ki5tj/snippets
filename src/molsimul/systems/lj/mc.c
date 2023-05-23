/* basic Monte Carlo simulation */
#include "lj.h"



int n = 108;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.7;
double tp = 1.5;
double rcdef = 1e9;
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos = "lj.pos";
int dopr = 0;



int main(void)
{
  int t, acc;
  lj_t *lj;
  double epsm = 0, accsm = 0;

  lj = lj_open(n, rho, rcdef, dopr);
  lj_energy(lj);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    acc = lj_metro(lj, amp, 1/tp);
    if ( t <= nequil ) continue;
    epsm += lj->epot;
    accsm += acc;
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  printf("rho %g, tp %g, ep %g, acc %g%%\n",
      rho, tp, epsm/nsteps/n, 100.*accsm/nsteps);
  lj_close(lj);
  return 0;
}


/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#ifndef D
#define D 3
#endif
#include "ljcore.h"



int n = 108;
int nequil = 10000;
int nsteps = 100000;
double rho = 0.7;
double tp = 1.5;
double rcdef = 2.5;
double dt = 0.002;
double thdt = 0.02;
const char *fnpos = "lj.pos";



int main(void)
{
  int t, tstat = 1;
  lj_t *lj;
  double epsm = 0;

  lj = lj_open(n, rho, rcdef);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_vv(lj, dt);
    if ( tstat ) {
      lj->ekin = lj_vrescale(lj, tp, thdt);
    } else { /* test energy conservation */
      lj->ekin = lj_ekin(lj->v, n);
      if ( t % 1000 == 0 ) printf("%d, ep %g, ek %g, e %g\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
    }
    if ( t <= nequil ) continue;
    epsm += lj->epot;
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}


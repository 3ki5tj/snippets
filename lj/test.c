/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif



int main(void)
{
  int n = 108, t, nequil = 10000, nsteps = 100000, tstat = 0;
  double rho = 0.7, rcdef = 2.5, dt = 0.002, thdt = 0.02, tp = 1.5;
  lj_t *lj;
  double epsm = 0;

  lj = lj_open(n, rho, rcdef);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_vv(lj, dt);
    if ( tstat ) {
      lj->ekin = lj_vrescale(lj, tp, thdt);
    } else {
      lj->ekin = lj_ekin(lj->v, n);
      if ( t % 1000 == 0 ) printf("%d, ep %g, ek %g, e %g\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
    }
    if ( t <= nequil ) continue;
    epsm += lj->epot;
  }
  lj_writepos(lj, lj->x, lj->v, "a.pos");
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}


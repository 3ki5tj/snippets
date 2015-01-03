#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif



int n = 108;
int nequil = 10000;
int nsteps = 100000;
double rho = 0.65;
double tp = 1.25;
double pres = 1.0;
double rcdef = 1e9;
double dt = 0.002;
double thdt = 0.02;
double pdt = 1e-5;
const char *fnpos = "lj.pos";



int main(void)
{
  int t;
  double usm = 0, psm = 0, rhosm = 0;
  lj_t *lj;

  lj = lj_open(n, rho, rcdef);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_vv(lj, dt);
    lj_rmcom(lj->v, n);
    lj->ekin = lj_vrescale(lj, tp, thdt);
    if ( t % 5 == 0 ) lj_langp0(lj, pdt, tp, pres, 0);
    if ( t % 1000 == 0 ) printf("t %d, T %g, rho %g, p %g, u %g\n", t, tp, lj->rho, lj_calcp(lj, tp), lj->epot);
    if ( t <= nequil ) continue;
    usm += lj->epot;
    psm += lj_calcp(lj, tp);
    rhosm += lj->rho;
  }
  lj_writepos(lj, lj->x, lj->v, "a.pos");
  lj_close(lj);
  printf("rho %g, tp %g, ep %g, p %g\n", rhosm/nsteps, tp, usm/nsteps/n, psm/nsteps);
  return 0;
}


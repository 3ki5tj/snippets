/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#include "lj.h"
#include "ljeos.h"



int n = 108;
int nequil = 10000;
int nsteps = 100000;
double rho = 0.7;
double tp = 1.5;
double rcdef = 2.5;
double dt = 0.002;
double thdamp = 10;
enum { NONE, VRESCALE, NHCHAIN, LANGEVIN };
int tstat = NHCHAIN;
#define NHCLEN 3
double zeta[]  = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double zmass[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
const char *fnpos = "lj.pos";
int dopr = 0;


void thermostat(lj_t *lj, double dt)
{
  if ( tstat == VRESCALE ) {
    lj->ekin = lj_vrescale(lj, tp, dt * thdamp);
  } else if ( tstat == NHCHAIN ) {
    lj->ekin = lj_nhchain(lj, tp, dt, NHCLEN, zeta, zmass);
  } else if ( tstat == LANGEVIN ) {
    lj->ekin = lj_langevin(lj, tp, dt * thdamp);
  } else { /* test energy conservation */
    lj->ekin = lj_ekin(lj->v, n);
  }
}


int main(void)
{
  int t;
  lj_t *lj;
  double epsm = 0, eksm = 0;

  lj = lj_open(n, rho, rcdef, dopr);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    thermostat(lj, dt * 0.5);
    lj_vv(lj, dt);
    thermostat(lj, dt * 0.5);
    if ( t % 1000 == 0 && !tstat )
      printf("%d, ep %g, ek %g, e %g\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
    if ( t <= nequil ) continue;
    epsm += lj->epot;
    eksm += lj->ekin;
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g(%g), ep %g (ref. %g)\n", rho,
      tp, 2*eksm/nsteps/lj->dof, epsm/nsteps/n,
      ljeos3d_get(rho, tp, NULL, NULL, NULL));
  return 0;
}


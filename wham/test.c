#include "wham.h"
#include "../lj/lj3d.h"


int n = 108;
double rho = 0.3;
double rcdef = 2.5;
double dt = 0.002; /* MD time step */
double thdt = 0.02; /* thermostat time step */
int nequil = 10000;
int nsteps = 100000;

int ntp = 10;
double xmin, xmax, dx;



int main(void)
{
  hist_t *hs;
  double *beta, *lnz, *epot;
  lj_t **lj;
  int itp, istep;

  xnew(lj, ntp);
  xnew(beta, ntp);
  xnew(lnz, ntp);
  xnew(epot, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    lj[itp] = lj_open(n, rho, rcdef);
    beta[itp] = 1./(1 + .1 * itp);
    lnz[itp] = epot[itp] = 0;
  }
  hs = hist_open(ntp, xmin, xmax, dx);

  /* do the simulations */
  for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
    for ( itp = 0; itp < ntp; itp++ ) {
      lj_vv(lj[itp], dt);
      lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], thdt);
      epot[itp] = lj[itp]->epot;
    }
    if ( istep <= nequil ) continue;
    hist_add(hs, epot, 1, 0);
  }

  wham(hs, beta, lnz, NULL, NULL);
  hist_close(hs);
  for ( itp = 0; itp < ntp; itp++ )
    lj_close( lj[itp] );
  free(beta);
  free(lnz);
  free(epot);
  return 0;
}


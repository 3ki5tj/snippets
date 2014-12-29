/* test program for the two-dimensional WHAM */
#define WHAM2_MDIIS
#include "wham2.h"
#include "../lj/lj3d.h"


int n = 108;
double rho = 0.7;
double rcdef = 2.5;
double dt = 0.002; /* MD time step */
double thdt = 0.02; /* thermostat time step */
int nequil = 4000;
int nsteps = 40000;

int ntp = 5;
double emin = -500, emax =   0, de = 1.0;
double vmin =  200, vmax = 500, dv = 1.0;
int itmax = 100000;
double tol = 1e-7;
int nbases = 5;

const char *fnlndos = "lndos2.dat";
const char *fneav = "eav2.dat";
const char *fnhist = "hist2.dat";

enum { METHOD_DIRECT = 0, METHOD_MDIIS = 1 };
int method = METHOD_DIRECT;



int main(int argc, char **argv)
{
  hist2_t *hs;
  double *beta, *bp, *lnz, *epot, *vols;
  lj_t **lj;
  int itp, istep;

  if ( argc > 1 ) {
    method = atoi(argv[1]);
  }

  xnew(lj, ntp);
  xnew(beta, ntp);
  xnew(bp, ntp);
  xnew(lnz, ntp);
  xnew(epot, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    lj[itp] = lj_open(n, rho, rcdef);
    beta[itp] = 1./(0.8 + .3 * itp);
    lnz[itp] = epot[itp] = 0;
  }
  hs = hist2_open(ntp, emin, emax, de, vmin, vmax, dv);

  /* try to load the histogram, if it fails, do simulations */
  if ( 0 != hist2_load(hs, fnhist, HIST_VERBOSE) ) {
    /* do the simulations */
    for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) {
        lj_vv(lj[itp], dt);
        lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], thdt);
        epot[itp] = lj[itp]->epot;
        vols[itp] = lj[itp]->vol;
      }
      if ( istep <= nequil ) continue;
      hist2_add(hs, epot, vols, 1, 1.0, 0);
    }

    fprintf(stderr, "simulation ended, doing WHAM\n");
    hist2_save(hs, fnhist, HIST_ADDAHALF);
  }

  if ( method == METHOD_DIRECT ) {
    wham2(hs, beta, bp, lnz, itmax, tol, fnlndos, fneav);
  } else {
    wham2_mdiis(hs, beta, bp, lnz, nbases, 1.0, itmax, tol, 0, fnlndos, fneav);
  }
  hist2_close(hs);
  for ( itp = 0; itp < ntp; itp++ )
    lj_close( lj[itp] );
  free(beta);
  free(lnz);
  free(epot);
  return 0;
}


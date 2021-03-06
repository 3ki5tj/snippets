/* test program for WHAM */
#define ENABLE_MDIIS
#include "../wham.h"
#define IS2_LB 6
#include "is2.h"
#include <time.h>



int nequil = 100000;
int nsteps = 10000000;
int verbose = 0;

int ntp = 80;
double xmin = -2*IS2_N - 4, xmax = 4, dx = 4;
int method = 0;
int itmax = 100000;
double tol = 1e-7;
int nbases = 5;

const char *fnlndos = "lndos.dat";
const char *fneav = "eav.dat";
const char *fnhist = "hist.dat";



int main(int argc, char **argv)
{
  hist_t *hs;
  double *beta, *lnz, *epot;
  unsigned (*uproba)[5];
  is2_t **is;
  int itp, istep;

  if ( argc > 1 ) {
    method = atoi(argv[1]);
  }

  xnew(is, ntp);
  xnew(beta, ntp);
  xnew(lnz, ntp);
  xnew(uproba, ntp);
  xnew(epot, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    is[itp] = is2_open(IS2_L);
    beta[itp] = 1./(1.5 + 0.02 * itp);
    is2_setuproba(beta[itp], uproba[itp]);
    lnz[itp] = epot[itp] = 0;
  }

  mtscramble( time(NULL) );
  hs = hist_open(ntp, xmin, xmax, dx);

  /* try to load the histogram, if it fails, do simulations */
  if ( 0 != hist_load(hs, fnhist, HIST_VERBOSE) ) {
    int id, h;

    /* do the simulations */
    for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) {
        IS2_PICK(is[itp], id, h);
        if ( h < 0 || mtrand() <= uproba[itp][h] ) {
          IS2_FLIP(is[itp], id, h);
        }
        epot[itp] = is[itp]->E;
      }
      if ( istep <= nequil ) continue;
      //for ( itp = 0; itp < ntp; itp++ ) printf("itp %d, ep %d\n", itp, is[itp]->E);
      //printf("hs->xmin %g\n", hs->xmin); getchar();
      hist_add(hs, epot, 1, 0);
    }

    hist_save(hs, fnhist, 0);
    fprintf(stderr, "simulation ended, doing WHAM\n");
  }

  whamx(hs, beta, lnz, 0, NULL, 1.0, nbases, 0, 10.0,
      0, itmax, tol, verbose, fnlndos, fneav, method);

  hist_close(hs);
  for ( itp = 0; itp < ntp; itp++ )
    is2_close( is[itp] );
  free(beta);
  free(lnz);
  free(uproba);
  free(epot);
  return 0;
}


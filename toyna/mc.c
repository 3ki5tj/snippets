/* basic Monte Carlo simulation */
#include "nacore.h"



const char seq[] = "ACGGUUCAGCU";
int nequil = 100000;
int nsteps = 1000000;
double tp = 300.0;
double debyel;
double amp = 0.2; /* Monte Carlo move size */
double qv[1] = {1.0}; /* ion charge of the solvent */
double conc[1] = {1.0}; /* solvent concentration */
const char *fnpos = "na.pos";

double mctot = 1e-30, mcacc = 0;



int main(void)
{
  int t;
  na_t *na;
  double epsm = 0;

  debyel = getDebyel(qv, conc, 1, tp);
  na = na_open(seq, 3, tp, debyel);
  //na_energy(na);
  //for ( t = 1; t <= nequil + nsteps; t++ ) {
  //  mctot += 1;
  //  mcacc += na_metro(na, amp, 1/(BOLTZ*tp));
  //  if ( t <= nequil ) continue;
  //  epsm += na->epot;
  //}
  na_writepos(na, na->x, na->v, fnpos);
  na_close(na);
  printf("nr %d, tp %g, ep %g, acc %g%%\n",
      na->nr, tp, epsm/nsteps, 100.*mcacc/mctot);
  return 0;
}


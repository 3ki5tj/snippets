/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#ifndef D
#define D 3
#endif

#include "nacore.h"



const char seq[] = "ACGGUUCAGCU";
int nequil = 10000;
int nsteps = 100000;
double tp = 300.0;
double debyel;
double dt = 0.002;
double thdt = 0.02;
double qv[1] = {1.0};
double conc[1] = {1.0};
const char *fnpos = "na.pos";



int main(void)
{
  int t, tstat = 1;
  na_t *na;
  double epsm = 0;

  debyel = getDebyel(qv, conc, 1, tp);
  na = na_open(seq, 3, tp, debyel);
  double Q, eps;
  Q = getchargeQ(tp, &eps);
  printf("Debyel %g, eps %g, chargeq %g\n", debyel, eps, Q);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    na_vv(na, dt);
    if ( tstat ) {
      na->ekin = na_vrescale(na, tp, thdt);
    } else { /* test energy conservation */
      na->ekin = na_ekin(na->v, na->m, na->na);
    }
    if ( t % 1000 == 0 ) printf("%d, ep %g, ek %g, e %g\n", t, na->epot, na->ekin, na->epot + na->ekin);
    if ( t <= nequil ) continue;
    epsm += na->epot;
  }
  na_writepos(na, na->x, na->v, fnpos);
  na_close(na);
  printf("nr %d, tp %g, ep %g\n", na->nr, tp, epsm/nsteps);
  return 0;
}


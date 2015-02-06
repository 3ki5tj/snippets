/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#ifndef D
#define D 3
#endif

#include "nacore.h"



int model = NA_TIS;
char *seq = "ACGGUUCAGCU";
double tp = 300.0;
double debyel;
double uhb0 = 2.43;
double dt = 0.002;
double thdt = 0.02;
double qv[1] = {1.0};
double conc[1] = {1.0};
const char *fnpos = "na.pos";
const char *fninp = "refs/1RNK.tis.pdb";

int nequil = 10000;
int nsteps = 1000000;
int nstrep = 10000;



int main(int argc, char **argv)
{
  int t, irmin = 1, tstat = 1;
  na_t *na;
  double epsm = 0, Q, eps;

  Q = getchargeQ(tp, &eps);
  debyel = getDebyel(qv, conc, 1, tp);
  if ( argc > 1 ) {
    fninp = argv[1];
  }
  if ( fninp != NULL ) {
    seq = na_getseqpdb(fninp, &irmin);
  }
  na = na_open(seq, model, tp, debyel, uhb0);
  if ( fninp ) {
    if ( na_loadpdb(na, fninp, irmin) != 0 ) {
      fprintf(stderr, "cannot load %s\n", fninp);
      return -1;
    }
    na_writepos(na, na->x, NULL, "init.pos");
    fprintf(stderr, "loaded %s\n", fninp); // getchar();
  }
  printf("%s, seq(%d): %s, Debye-l %g, eps %g, chargeq %g\n",
      fninp, na->nr, na->seq, debyel, eps, Q);

  /* main molecular dynamics loop */
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    na_vv(na, dt);
    if ( tstat ) {
      na->ekin = na_vrescale(na, tp, thdt);
    } else {
      na->ekin = na_ekin(na, na->v);
    }
    if ( t % nstrep == 0 ) {
      fprintf(stderr, "%d, ep %g, ek %g, e %g\n",
          t, na->epot, na->ekin, na->epot + na->ekin);
      na_writepos(na, na->x, na->v, fnpos);
    }
    if ( t <= nequil ) continue;
    epsm += na->epot;
  }

  na_close(na);
  printf("nr %d, tp %g, ep %g\n", na->nr, tp, epsm/nsteps);
  return 0;
}


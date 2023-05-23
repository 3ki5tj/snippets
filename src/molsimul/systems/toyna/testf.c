/* check force */
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
double qv[1] = {1.0};
double conc[1] = {1.0};
const char *fninp = "refs/1RNK.tis.pdb";



int main(int argc, char **argv)
{
  int i, irmin = 1;
  na_t *na;
  double Q, eps;
  double e0, e1, f2, del = 0.001;

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
  }
  printf("%s, seq(%d): %s, Debye-l %g, eps %g, chargeq %g\n",
      fninp, na->nr, na->seq, debyel, eps, Q);

  /* compute the energy and force */
  e0 = na_force(na);
  /* measure the force */
  for ( f2 = 0, i = 0; i < na->na; i++ ) {
    f2 += vsqr( na->f[i] );
  }
  /* move along the force */
  for ( i = 0; i < na->na; i++ ) {
    vsinc(na->x[i], na->f[i], del/f2);
  }
  /* compute the energy again */
  e1 = na_force(na);

  fprintf(stderr, "energy %g -> %g(%g): %g\n",
      e0, e1, na_energy(na), (e0 - e1) / del);

  na_close(na);
  return 0;
}


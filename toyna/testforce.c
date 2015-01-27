/* check force */
#ifndef D
#define D 3
#endif

#include "nacore.h"



const char seq[] = "ACGGUUCAGCU";
double tp = 300.0;
double debyel;
double rc = 1e9;
double dt = 0.002;



int main(void)
{
  int i, t;
  na_t *na;
  double q[1] = {1.0};
  double conc[1] = {1.0};
  double Q, eps;
  double e0, e1, f2, del = 0.001;

  na = na_open(seq, rc, 3);
  debyel = getDebyel(q, conc, 1, tp);
  Q = getchargeQ(tp, &eps);
  printf("Debyel %g, eps %g, chargeq %g\n", debyel, eps, Q);
  for ( t = 1; t <= 1000; t++ ) {
    na_vv(na, dt, tp, debyel);
  }

  /* compute the energy and force */
  e0 = na_force(na, tp, debyel);
  /* measure the force */
  for ( f2 = 0, i = 0; i < na->na; i++ ) {
    f2 += vsqr( na->f[i] );
  }
  /* move along the force */
  for ( i = 0; i < na->na; i++ ) {
    vsinc(na->x[i], na->f[i], del/f2);
  }
  /* compute the energy again */
  e1 = na_force(na, tp, debyel);

  fprintf(stderr, "energy %g -> %g(%g): %g\n",
      e0, e1, na_energy(na, tp, debyel),
      (e0 - e1) / del);

  na_close(na);
  return 0;
}


/* check force */
#ifndef D
#define D 3
#endif

#include "nacore.h"



const char seq[] = "ACGGUUCAGCU";
double tp = 300.0;
double debyel;
double uhb0 = 2.43;
double dt = 0.002;
double qv[1] = {1.0};
double conc[1] = {1.0};



int main(void)
{
  int i, t;
  na_t *na;
  double Q, eps;
  double e0, e1, f2, del = 0.001;

  debyel = getDebyel(qv, conc, 1, tp);
  na = na_open(seq, 3, tp, debyel, uhb0);
  Q = getchargeQ(tp, &eps);
  printf("Debyel %g, eps %g, chargeq %g\n", debyel, eps, Q);
  for ( t = 1; t <= 1000; t++ ) {
    na_vv(na, dt);
  }

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


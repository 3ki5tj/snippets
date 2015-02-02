#define IS2_LB 5
#define L (1 << IS2_LB)
#include "is2.h"



int main(void)
{
  is2_t *is;
  int id, h;
  double beta = 1/2.0;
  unsigned long t, nsteps = 10000000L;
  double sE = 0, eavref, cvref;

  is = is2_open(L);
  IS2_SETPROBA(is, beta);
  for ( t = 0; t < nsteps; t++ ) {
    IS2_PICK(is, id, h);
    if ( h < 0 || mtrand() <= is->uproba[h] ) {
      IS2_FLIP(is, id, h);
    }
    sE += is->E;
  }
  is2_exact(is, beta, &eavref, &cvref);
  printf("E %g vs. %g\n", sE/nsteps, eavref);
  is2_close(is);
  return 0;
}

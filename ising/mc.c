#define IS2_LB 5
#define L (1 << IS2_LB)
#include "is2.h"


static void runmetro(is2_t *is, double tp)
{
  int id, h;
  double beta = 1/tp;
  unsigned long t, nsteps = 100000000L;
  double sE = 0, sEE = 0, eavref, cvref;

  IS2_SETPROBA(is, beta);
  for ( t = 0; t < nsteps; t++ ) {
    IS2_PICK(is, id, h);
    if ( h < 0 || mtrand() <= is->uproba[h] ) {
      IS2_FLIP(is, id, h);
    }
    sE += is->E;
    sEE += is->E * is->E;
  }
  is2_exact(is->l, is->l, beta, &eavref, &cvref);
  h = is->E;
  sE /= nsteps;
  sEE = sEE / nsteps - sE * sE;
  printf("E: average %g vs. %g; final %d vs %d; Cv %g vs. %g\n",
      sE, eavref, h, is2_em(is), sEE*(beta*beta), cvref);
}



static void runwolff(is2_t *is, double tp)
{
  int h;
  double beta = 1/tp, padd;
  unsigned long t, nsteps = 50000L;
  double sE = 0, sEE = 0, eavref, cvref;

  padd = 1 - exp(-2*beta);
  for ( t = 0; t < nsteps; t++ ) {
    is2_wolff(is, padd);
    sE += is->E;
    sEE += 1.0 * is->E * is->E;
  }
  is2_exact(is->l, is->l, beta, &eavref, &cvref);
  h = is->E;
  sE /= nsteps;
  sEE = sEE / nsteps - sE * sE;
  printf("E: average %g vs. %g; final %d vs %d; Cv %g vs. %g\n",
      sE, eavref, h, is2_em(is), sEE*(beta*beta), cvref);
}



int main(int argc, char **argv)
{
  is2_t *is;
  int method = 0;
  double tp = 2.269;

  if ( argc > 1 ) {
    method = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    tp = atof( argv[2] );
  }

  is = is2_open(L);
  if ( method == 0 ) {
    runmetro(is, tp);
  } else {
    runwolff(is, tp);
  }
  is2_close(is);
  return 0;
}


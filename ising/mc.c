#define IS2_LB 5
#define L IS2_L
#include "is2.h"
#include "ave.h"



static void run(is2_t *is, double tp, int method, long nsteps)
{
  int id, h;
  double beta = 1/tp, padd;
  long t;
  double aveE, varE, eavref, cvref;
  av_t avE[1];

  av_clear(avE);
  if ( method == 0 ) {
    is2_setuproba(beta, is->uproba);
  } else {
    padd = 1 - exp(-2*beta);
  }
  for ( t = 0; t < nsteps; t++ ) {
    if ( method == 0 ) { /* Metropolis algorithm */
      //id = is2_pick(is, &h);
      IS2_PICK(is, id, h);
      if ( h <= 0 || mtrand() <= is->uproba[h] ) {
        //is2_flip(is, id, h);
        IS2_FLIP(is, id, h);
      }
    } else { /* Wolff cluster algorithm */
      is2_wolff(is, padd);
    }
    av_add(avE, is->E);
  }
  is2_exact(is->l, is->l, beta, &eavref, &cvref);
  h = is->E;
  aveE = av_getave(avE, &varE);
  printf("E: average %g vs. %g; final %d vs %d; Cv %g vs. %g, stdE %g\n",
      aveE, eavref, h, is2_em(is), varE*(beta*beta), cvref, sqrt(varE));
}



int main(int argc, char **argv)
{
  is2_t *is;
  int method = 0;
  long nsteps = 0;
  double tp = 2.269;

  if ( argc > 1 ) method = atoi( argv[1] );
  if ( argc > 2 ) tp = atof( argv[2] );
  if ( argc > 3 ) nsteps = atol( argv[3] );
  if ( nsteps <= 0 )
    nsteps = ( method == 0 ) ? 100000000L : 50000;

  is = is2_open(L);
  run(is, tp, method, nsteps);
  is2_close(is);
  return 0;
}


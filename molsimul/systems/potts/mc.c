#define POTTS2_LB 5
#define L POTTS2_L
#include "potts2.h"
#include "ave.h"



static void run(potts2_t *pt, double tp, int method, long nsteps)
{
  int id, h, sn = 0;
  double beta = 1/tp, padd;
  long t;
  double aveE, varE;
  av_t avE[1];

  av_clear(avE);
  if ( method == 0 ) {
    potts2_setuproba(beta, pt->uproba);
  } else {
    padd = 1 - exp(-beta);
  }
  for ( t = 0; t < nsteps; t++ ) {
    if ( method == 0 ) { /* Metropolis algorithm */
      //id = potts2_pick(pt, &sn, &h);
      POTTS2_PICK(pt, id, sn, h);
      if ( h <= 0 || mtrand() <= pt->uproba[h] ) {
        //potts2_flip(pt, id, sn, h);
        POTTS2_FLIP(pt, id, sn, h);
      }
      //printf("%d %d %d %d\n", id, h, sn, pt->E);
      //printf("%d\n", potts2_energy(pt));
      //getchar();
    } else { /* Wolff cluster algorithm */
      potts2_wolff(pt, padd);
    }
    av_add(avE, pt->E);
  }
  h = pt->E;
  aveE = av_getave(avE, &varE);
  printf("E: average %g (E/N %g); final %d vs %d; Cv %g, stdE %g\n",
      aveE, aveE/pt->n, h, potts2_energy(pt), varE*(beta*beta), sqrt(varE));
}



int main(int argc, char **argv)
{
  potts2_t *pt;
  int method = 0, q = 10;
  long nsteps = 0;
  double tp = 0.7;

  if ( argc > 1 ) q = atoi( argv[1] );
  if ( argc > 2 ) method = atoi( argv[2] );
  if ( argc > 3 ) tp = atof( argv[3] );
  if ( argc > 4 ) nsteps = atol( argv[4] );
  if ( nsteps <= 0 )
    nsteps = ( method == 0 ) ? 100000000L : 500000;

  pt = potts2_open(L, q);
  run(pt, tp, method, nsteps);
  potts2_close(pt);
  return 0;
}


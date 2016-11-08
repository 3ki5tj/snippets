/* compute the electric energy of a gaussian ion
 * in a uniform negative background by Ewald sum
 * To compile and run
 *    gcc ion.c -lm && ./a.out
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double ewald(double sigma, int km)
{
  int i, j, l, nz;
  double k2, elimit, erecip = 0, etot, mul, del;

  /* reciprocal-space sum */
  for ( i = 0; i <= km; i++ ) {
    for ( j = 0; j <= km; j++ ) {
      for ( l = 0; l <= km; l++ ) {
        k2 = i*i + j*j + l*l;
        nz = (i == 0) + (j == 0) + (l == 0);
        if ( nz == 3 ) {
          continue;
        } else if ( nz == 2 ) {
          mul = 2;
        } else if ( nz == 1 ) {
          mul = 4;
        } else {
          mul = 8;
        }
        erecip += mul * exp(-2*k2*sigma*sigma*M_PI)/k2;
      }
    }
  }
  erecip /= M_PI;

  elimit = -sqrt(2)/sigma;
  etot = erecip + elimit;
  del = (etot + 2.8372974794806)/(sigma*sigma);
  printf("tot %.14f, elimit %.14f, recip %.14f, del %g, sigma %g, kmax %d\n",
      etot, elimit, erecip, del, sigma, km);
  return elimit + erecip;
}


int main(int argc, char **argv)
{
  double sigma = 0.01;
  int kterms = 200;

  if ( argc > 1 ) sigma = atof(argv[1]);
  if ( argc > 2 ) kterms = atof(argv[2]);
  ewald(sigma, kterms);
  return 0;
}

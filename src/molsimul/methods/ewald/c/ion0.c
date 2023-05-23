/* compute the electric energy of a single ion
 * in a uniform negative background by Ewald sum
 * To compile and run
 *    gcc ion.c -lm && ./a.out
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double ewald(double sigma, int xm, int km)
{
  int i, j, l;
  double r, k2, eself = 0, ereal = 0, erecip = 0, ebg = 0;

  /* real-space sum */
  for ( i = -xm; i <= xm; i++ ) {
    for ( j = -xm; j <= xm; j++ ) {
      for ( l = -xm; l <= xm; l++ ) {
        r = sqrt(i*i + j*j + l*l);
        if ( r <= 0 ) continue;
        ereal += erfc(sqrt(M_PI)*r/sigma)/r;
      }
    }
  }

  /* reciprocal-space sum */
  for ( i = -km; i <= km; i++ ) {
    for ( j = -km; j <= km; j++ ) {
      for ( l = -km; l <= km; l++ ) {
        k2 = i*i + j*j + l*l;
        if ( k2 <= 0 ) continue;
        erecip += exp(-k2*sigma*sigma*M_PI)/k2;
      }
    }
  }
  erecip /= M_PI;

  /* self energy */
  eself = -2/sigma;

  /* background energy */
  ebg = -sigma*sigma;
  printf("real %g, recip %g, self %g, background %g\n", ereal, erecip, eself, ebg);
  return eself + ereal + erecip + ebg;
}


int main(int argc, char **argv)
{
  double ene, sigma = 1.2;
  int nterms = 20, kterms = 20;

  if ( argc > 1 ) sigma = atof(argv[1]);
  if ( argc > 2 ) nterms = atoi(argv[2]);
  if ( argc > 3 ) kterms = atoi(argv[3]);
  ene = ewald(sigma, nterms, kterms);
  printf("%.15lf %g %d %d\n", ene, sigma, nterms, kterms);
  return 0;
}

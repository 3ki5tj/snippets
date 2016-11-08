/* compute the electric energy of a uniform-sphere ion
 * in a uniform negative background by Ewald sum
 * To compile and run
 *    gcc ion.c -lm && ./a.out
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int solid = 1;

static double ewald(double r, int km)
{
  int i, j, l, nz;
  double k2, elimit, erecip = 0, etot, mul, del, x, rhok;

  /* reciprocal-space sum */
  for ( i = 0; i <= km; i++ ) {
    for ( j = 0; j <= km; j++ ) {
      for ( l = 0; l <= km; l++ ) {
        k2 = i*i + j*j + l*l;
        if ( k2 > km * km ) continue;
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
        x = sqrt(k2) * 2 * M_PI * r;
        if ( solid ) {
          rhok = 3 * (sin(x) - x * cos(x)) / (x * x * x);
        } else {
          rhok = sin(x) / x;
        }
        erecip += mul * rhok * rhok / k2;
      }
    }
  }
  erecip /= M_PI;

  if ( solid ) {
    elimit = -1.2/r;
  } else {
    elimit = -1.0/r;
  }
  etot = erecip + elimit;
  del = (etot + 2.8372974794806)/(r*r*M_PI);
  printf("%s: tot %.14f, elimit %.14f, recip %.14f, del %g*PI, r %g, kmax %d\n",
      (solid ? "solid" : "hollow"), etot, elimit, erecip, del, r, km);
  return elimit + erecip;
}


int main(int argc, char **argv)
{
  double r = 0.1;
  int kterms = 200;

  if ( argc > 1 ) r = atof(argv[1]);
  if ( argc > 2 ) kterms = atof(argv[2]);
  if ( argc > 3 ) solid = 0;
  ewald(r, kterms);
  return 0;
}

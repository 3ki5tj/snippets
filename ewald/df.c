/* compute the electric energy of a spherical ion
 * in a uniform negative background by Ewald sum
 * To compile and run
 *    gcc ion.c -lm && ./a.out
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cisi.h"


double r = 0.2; /* radius of the cloud charge */
int kterms = 200;

enum { HOLLOW_SPHERE, SOLID_SPHERE, GAUSSIAN, EXPONENTIAL };
const char *type_names[] = { "hollow sphere", "solid sphere", "Gaussian", "exponential" };
int type = SOLID_SPHERE;


/* compute the finite size correction of
 * twice the electrostatic energy of an ion in a unit cubic box */
static double ewald(double rB, int km)
{
  int i, j, l, nz;
  double r, rr, k2, mul, del, x, rhok, xx, ci, si;
  double erecip = 0, elimit, etail = 0, ediff;
  const double PI = 3.1415926535897932385;

  x = km * 2 * PI; /* maximal wave vector */
  /* compute the real radius r and the average of < r^2 > */
  if ( type == HOLLOW_SPHERE ) {
    r = rB;
    rr = r*r;
    x *= r;
    cisi(2*x, &ci, &si);
    etail = 2/(PI*r)*((1-cos(2*x))/(2*x) + PI/2 - si);
  } else if ( type == SOLID_SPHERE ) {
    r = 6./5*rB;
    rr = 3./5*r*r;
    x *= r;
    cisi(2*x, &ci, &si);
    xx = x * x;
    etail = 3/(5*PI*r)*((3+5*xx-(3-xx+2*xx*xx)*cos(2*x)-x*(6+xx)*sin(2*x))/(x*xx*xx) + 2*PI - 4*si);
  } else if ( type == GAUSSIAN ) {
    r = sqrt(2/PI)*rB;
    rr = 1.5*r*r;
    etail = 1/rB*erfc(x*rB/sqrt(PI));
  } else if ( type == EXPONENTIAL ) {
    r = 0.5*rB;
    rr = 6*r*r;
    /* kappa / pi (theta + sin(2*theta)/2)|_{theta_C}^{pi/2} */
    x = atan(x*r);
    etail = 1/rB/PI*(PI-2*x-sin(2*x));
  }

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
        x = sqrt(k2) * 2 * PI * r;
        if ( type == SOLID_SPHERE ) {
          rhok = 3 * (sin(x) - x * cos(x)) / (x * x * x);
          // approximately exp(-x*x/10)
        } else if ( type == HOLLOW_SPHERE ) {
          rhok = sin(x) / x;
          // approximately exp(-x*x/6)
        } else if ( type == GAUSSIAN ) {
          rhok = exp(-x*x/4);
        } else if ( type == EXPONENTIAL ) {
          rhok = 1/(1 + x*x);
        }
        erecip += mul * rhok * rhok / k2;
      }
    }
  }
  /* In Gaussian and natural units, the coefficient is
   * (1/eps) / (2*pi)^2 = (4*pi) / (2*pi)^2 = 1/pi */
  erecip /= PI;

  erecip += etail;

  /* continuous limit */
  elimit = 1./rB;
  ediff = erecip - elimit;
  del = (ediff + 2.8372974794806)/(rr*PI);
  printf("%12s: diff %.14f, elimit %.14f, recip %.14f, del %8.6f pi, tail %.6f*pi, r %8.6f, <r^2> %8.6f, kmax %d\n",
      type_names[type], ediff, elimit, erecip, del, etail/rr/PI, r, rr, km);
  return erecip;
}



int main(int argc, char **argv)
{
  if ( argc > 1 ) type = atoi(argv[1]);
  if ( argc > 2 ) r = atof(argv[2]);
  if ( argc > 3 ) kterms = atof(argv[3]);
  ewald(r, kterms);
  return 0;
}

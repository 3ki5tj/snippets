/* test the Eqs. (5) and (6) of Hub 2010 */
#include "mtrand.h"


#define NB 20  /* number of bins */
double h[NB]; /* histogram */


int n = 1000000;
double dx = 1.0/NB; /* bin size */


int main(void)
{
  double x, y, a, b;
  int i;

  a = 0.6;
  b = sqrt(1 - a * a);

  x = randgaus();
  for ( i = 0; i < n; i++ ) {
    x = a * x + b * randgaus();
    y = (erf(x/sqrt(2)) + 1) / 2;
    h[(int) (y / dx)] += 1;
  }

  printf("histogram:\n");
  for ( i = 0; i < NB; i++ ) {
    printf("%.2f %8.3f\n", i * dx, h[i]/(n*dx));
  }

  return 0;
}

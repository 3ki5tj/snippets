#include "mtrand.h"

double xmin = 1;
double xmax = 2;
double dx = 0.01;
#define N 100
double hist[N] = {0};
const double a = 0.1;

int main(void)
{
  int t = 0, i;
  double x = 1.5, nx, r, nr, xp;

  for ( t = 0; t < 100000000; t++ ) {
    r = randgaus();
    nx = x + x * r * a;
    if ( nx >= xmin && nx < xmax ) {
      nr = x * r / nx;
      xp = (r*r - nr*nr)/2 + log(x/nx);
      if ( xp > 0 || rand01() < exp(xp) ) {
        x = nx;
      }
    }
    i = (int)((x - xmin)/dx);
    //printf("%g %d\n", x, i); getchar();
    hist[i] += 1;
  }
  for ( i = 0; i < N; i++ )
    printf("%g\t%g\n", xmin+(i+0.5)*dx, hist[i]);
  return 0;
}

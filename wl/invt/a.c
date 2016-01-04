#include "mtrand.h"
#include "util.h"


/* the updating magnitude is c/t
 * n in the number bins
 * tau is the correlation time
 * tmax is the number of trials
 * */
static double comperr(double c, double amax,
    int n, double tau, long tmax, double frac)
{
  double *v, a, dv, err;
  int i, j, t;

  xnew(v, n);
  for ( i = 0; i < n; i++ ) {
    v[i] = 0. * amax * randgaus();
  }

  for ( t = 0; t < tmax; t++ ) {
    if ( i < 0 || rand01() * tau < 1 ) {
      j = (int) ( n * rand01() );
      dv = v[j] - v[i];
      if ( dv < 0 || rand01() < exp(-dv) ) {
        i = j;
      }
    }

    a = c * n / (t + 1);
    if ( a > amax ) a = amax;

    if ( frac > 0 ) {
      double rem = 1;

      if ( i > 0 ) {
        v[i-1] += a * frac;
        rem -= frac;
      }
      if ( i < n - 1 ) {
        v[i+1] += a * frac;
        rem -= frac;
      }
      v[i] += a * rem;
        
    } else {
      v[i] += a;
    }
  }

  /* normalize */
  for ( a = 0, i = 0; i < n; i++ ) {
    a += v[i];
  }
  a /= n;
  for ( i = 0; i < n; i++ ) {
    v[i] -= a;
  }

  /* compute the error */
  err = 0;
  for ( i = 0; i < n; i++ ) {
    err += v[i] * v[i];
  }
  return sqrt(err / n);
}



static double invt_test(int n, double c, double frac)
{
  double err, se = 0;
  int i;

  mtscramble( time(NULL) );
  for ( i = 0; i < n; i++ ) {
    err = comperr(c, 0.001, 3, 10.0, 1000000000L, frac);
    se += err;
    printf("%d: err %10.8f, ave %10.8f\n", i, err, se/(i+1));
  }

  err = se / n;
  printf("average error: %g\n", err);
  return err;
}


int main(int argc, char **argv)
{
  int n = 10;
  double c = 1.0, frac = 0.0;

  if ( argc > 1 ) {
    n = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    c = atof( argv[2] );
  }
  if ( argc > 3 ) {
    frac = atof( argv[3] );
  }

  printf("n %d, c %g, frac %g\n", n, c, frac);
  invt_test(n, c, frac);
  return 0;
}

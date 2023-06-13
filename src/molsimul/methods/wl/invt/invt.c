#include "rand/mtrand.h"
#include "util/util.h"



/* return the root-mean-squared error of the inverse time scheme
 *
 * the updating magnitude is `c/t` for t > `nequil`, or `a0` otherwise
 * `n` is the number of bins
 * `tau` is the correlation time
 * `nsteps` is the number of steps
 * */
static double comperr(double c, double alpha0, double nequil,
    int n, double tau, long nsteps,
    double frac,
    int verbose)
{
  double *v, a, dv, err, t0;
  int i, j, t;

  xnew(v, n);
  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  /* constant */
  t0 = c / alpha0;

  for ( t = 0; t < nsteps + nequil; t++ ) {
    /* MCMC sampling */
    if ( i < 0 || rand01() * tau < 1 ) {
      j = (int) ( n * rand01() );
      dv = v[j] - v[i];
      if ( dv < 0 || rand01() < exp(-dv) ) {
        i = j;
      }
    }

    /* compute the updating magnitude
     * the distribution density is p = 1/n
     * this is why we have to multiply by n */
    if ( t >= nequil ) {
      /* the constant t0 makes the transition
       * of alpha at t = t0 smooth */
      a = c * n / (t - nequil + t0);
    } else {
      a = alpha0 * n;
    }

    if ( frac > 0 ) {
      double rem = 1;

      if ( i > 0 ) {
        v[i - 1] += a * frac;
        rem -= frac;
      }
      if ( i < n - 1 ) {
        v[i + 1] += a * frac;
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

  /* print out the error */
  if ( verbose ) {
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", v[i]);
    }
    printf("\n");
  }

  return sqrt(err / n);
}



static double invt_test(int ntrials, int n, double c, double frac)
{
  double err, se = 0;
  int i;

  mtscramble( time(NULL) );
  for ( i = 0; i < ntrials; i++ ) {
    err = comperr(c, 0.01, n*10000,
        n, 1.0, 1000000000L, frac, 1);
    se += err;
    printf("%4d: err %10.8f, ave %10.8f\n", i, err, se/(i+1));
  }

  err = se / ntrials;
  printf("average error: %g\n", err);
  return err;
}



int main(int argc, char **argv)
{
  int ntrials = 10, n = 10;
  double c = 1.0, frac = 0.0;

  if ( argc > 1 ) {
    ntrials = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    n = atoi( argv[2] );
  }
  if ( argc > 3 ) {
    c = atof( argv[3] );
  }
  if ( argc > 4 ) {
    frac = atof( argv[4] );
  }

  printf("ntrials %d, n %d, c %g, frac %g\n", ntrials, n, c, frac);
  invt_test(ntrials, n, c, frac);
  return 0;
}

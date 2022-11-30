#ifndef RNG_H__
#define RNG_H__


#include "rng_engine_manager.h"



typedef struct {
  rng_engine_manager_t* rem;
} rng_t;


__inline static rng_t* rng_open(int type, uint32_t seed)
{
  rng_t* rng;

  if ((rng = (rng_t*) calloc(1, sizeof(*rng))) == NULL) {
    exit(1);
  }

  rng->rem = rng_engine_manager_open(type, seed);

  return rng;
}



__inline static void rng_close(rng_t* rng)
{
  rng_engine_manager_close(rng->rem);
  free(rng);
}



__inline static int rng_save(rng_t* rng, const char* fn)
{
  return rng_engine_manager_save(rng->rem, fn);
}



__inline static int rng_load(rng_t* rng, const char* fn)
{
  return rng_engine_manager_load(rng->rem, fn);
}



__inline static double rng_rand01(rng_t* rng)
{
  return rng_engine_manager_rand01(rng->rem);
}



/* Gaussian distribution with zero mean and unit variance
 * using the ratio method */
__inline static double rng_randgaus(rng_t* rng)
{
  double x, y, u, v, q;

  do {
    u = 1 - rng_rand01(rng);
    v = 1.7156*(rng_rand01(rng) - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));

  return v/u;
}



/* return a random number that satisfies the gamma distribution
 * p(x) = x^(k - 1) exp(-x) / Gamma(k) */
__inline static double rng_randgam(rng_t* rng, double k)
{
  int lt1 = 0;
  double a, b, x, v, u;

  if ( k <= 0 ) return 0;
  if ( k < 1 ) {
    lt1 = 1;
    k += 1;
  }
  a = k - 1./3;
  b = 1./3/sqrt(a);

  for ( ; ; ) {
    do {
      x = rng_randgaus(rng);
      v = 1 + b * x;
    } while ( v <= 0 );
    v *= v * v;
    x *= x;
    u = rng_rand01(rng);
    if ( u <= 1 - 0.331 * x * x ) break;
    u = log(u);
    if ( u <= 0.5 * x + a * (1 - v + log(v)) ) break;
  }

  x = a * v;
  if ( lt1 ) x *= pow(1 - rng_rand01(rng), 1./(k - 1));
  return x;
}



/* return a random number that satisfies the chi-squared distribution,
 * which is the sum of the squares k Gaussian random numbers */
__inline static double rng_randchisqr(rng_t* rng, double k)
{
  return 2*rng_randgam(rng, k*.5);
}



/* a randomly oriented unit vector */
__inline static double* rng_randdir(rng_t* rng, double *v)
{
  double a, b, sq, s;

  do {
    a = 2 * rng_rand01(rng) - 1;
    b = 2 * rng_rand01(rng) - 1;
    sq = a * a + b * b;
  } while ( sq >= 1 );
  s = 2. * sqrt(1 - sq);
  v[0] = a * s;
  v[1] = b * s;
  v[2] = 1 - 2 * sq;
  return v;
}



/* randomly accept by probability exp(-de) */
__inline static int rng_metroacc(rng_t* rng, double de)
{
  double r;
  if ( de <= 0 ) return 1;
  r = rng_rand01(rng);
  return r <= exp(-de);
}


#endif /* RNG_H__ */

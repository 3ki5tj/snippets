#include "mtrand.h"
#include "../hist/hist.h"



/* return a random number that satisfies the gamma distribution
 * p(x) = x^(k - 1) exp(-x) / (k - 1)!
 * old version */
__inline static double randgam_old(int k)
{
  int i;
  double x, k1 = k - 1, r, y, v1, v2, w;

  if ( k <= 0 ) return 0;
  if ( k <= 7 ) {
    /* adding random numbers that satisfy the exponential distribution */
    for ( x = 1.0, i = 0; i < k; i++ )
      x *= 1 - rand01();
    return -log(x);
  }

  w = sqrt(2.*k - 1);
  /* use the rejection method based on the Lorentz distribution */
  for (;;) {
    /* the Lorentz disribution is centered at k1, with width w
     * p(y) = 1/pi/(1 + y^2), x = y*w + k1
     * Int p(y) dy = 1/2 + arctan(y)/pi */
    for (;;) {
      v1 = 2 * rand01() - 1;
      v2 = 2 * rand01() - 1;
      if ( v1 * v1  + v2 * v2 < 1 ) {
        y = v2 / v1;
        x = w * y + k1;
        if (x > 0.) break;
      }
    }
    r = (1 + y*y) * exp(k1 * log(x/k1) - x + k1);
    if ( rand01() <= r ) break;
  }

  return x;
}



static void cmpgam(void)
{
  int i, K = 5;
  hist_t *h;
  double x[2];

  h = hist_open(2, 0, K * 10, 0.01);
  for ( i = 0; i < 10000000; i++ ) {
    x[0] = randgam(K);
    x[1] = randgam_old(K);
    hist_add(h, x, 1, 0);
  }
  hist_save(h, "gam.dat", HIST_ADDAHALF);
  hist_close(h);
}



static double randchisqr_simple(int K)
{
  double x, s = 0;
  int i;

  for ( i = 0; i < K; i++ ) {
    x = gaussrand();
    s += x * x;
  }
  return s;
}



static void cmpchisqr(void)
{
  int i, K = 5;
  hist_t *h;
  double x[2];

  h = hist_open(2, 0, K * 10, 0.01);
  for ( i = 0; i < 10000000; i++ ) {
    x[0] = randchisqr(K);
    x[1] = randchisqr_simple(K);
    hist_add(h, x, 1, 0);
  }
  hist_save(h, "chisqr.dat", HIST_ADDAHALF);
  hist_close(h);
}



int main(void)
{
  cmpgam();
  cmpchisqr();
  return 0;
}

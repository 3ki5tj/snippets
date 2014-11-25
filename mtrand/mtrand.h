#ifndef MTRAND_H__
#define MTRAND_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Mersenne Twister was developed by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000u /* most significant w-r bits */
#define MT_LMASK 0x7fffffffu /* least significant r bits */

int mtonce = 0;
int mtidx = MT_N;
unsigned long mtarr[MT_N] = {5491};

/* scramble the random number state */
__inline static void mtscramble(unsigned seed)
{
  int k;

  mtarr[0] = (seed * 314159265u + 271828183u) & 0xfffffffful;
  for (k = 1; k < MT_N; k++) { /* the final mask is for 64-bit machines */
    mtarr[k] = 1812433253ul * (mtarr[k - 1] ^ (mtarr[k - 1] >> 30)) + k;
    /* mr->arr[k] = (mr->arr[k] + seed) * 22695477ul + 1ul; */
    mtarr[k] = ((mtarr[k] + seed) * 314159265ul + 1ul) & 0xfffffffful;
  }
  mtidx = MT_N; /* request for an update */
  mtonce = 1; /* scrambled */
}

/* return an unsigned random number */
__inline static unsigned mtrand(void)
{
  static const unsigned mag01[2] = {0, 0x9908b0dfu}; /* MATRIX_A */
  unsigned x;
  int k;

  if ( !mtonce ) mtscramble(12345);

  if (mtidx >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+MT_M] ^ (x>>1) ^ mag01[x&1u];
    }
    for (; k < MT_N-1; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1u];
    }
    x = (mtarr[MT_N-1] & MT_UMASK) | (mtarr[0] & MT_LMASK);
    mtarr[MT_N-1] = mtarr[MT_M-1] ^ (x>>1) ^ mag01[x&1u];
    mtidx = 0;
  }
  x = mtarr[ mtidx++ ];
  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680u;
  x ^= (x << 15) & 0xefc60000u;
  x ^= (x >> 18);
  return x;
}


__inline static double rand01(void)
{
  return mtrand() / 4294967296.0;
}



/* Gaussian distribution with zero mean and unit variance
 * using ratio method */
__inline double gaussrand(void)
{
  double x, y, u, v, q;
  do {
    u = 1 - rand01();
    v = 1.7156*(rand01() - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}



/* a randomly oriented unit vector */
__inline static double *randdir(double *v)
{
  double a, b, sq, s;

  do {
    a = 2 * rand01() - 1;
    b = 2 * rand01() - 1;
    sq = a * a + b * b;
  } while ( sq >= 1 );
  s = 2. * sqrt(1 - sq);
  v[0] = a * s;
  v[1] = b * s;
  v[2] = 1 - 2 * sq;
  return v;
}

#endif

#ifndef MTRNG_H__
#define MTRNG_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>



/* Mersenne Twister was developed by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000u /* most significant w-r bits */
#define MT_LMASK 0x7fffffffu /* least significant r bits */



typedef struct {
  int index;
  uint32_t arr[MT_N];
  int inited;
} mtrng_t;



#ifdef MTRNG_USE_RNG_ALIASES__

typedef mtrng_t rng_t;
#define rng_open(seed) mtrng_open(seed)
#define rng_close(rng) mtrng_close(rng)
#define rng_save(rng, fn) mtrng_save(rng, fn)
#define rng_load(rng, fn) mtrng_load(rng, fn)
#define rng_randuint32(rng) mtrng_randuint32(rng)
#define rng_rand01(rng) mtrng_rand01(rng)
#define rng_randgaus(rng) mtrng_randgaus(rng)

#endif /* defined(MTRNG_USE_RNG_ALIASES__) */



/* scramble the random number state from a given seed */
__inline static void mtrng_init_from_seed(mtrng_t* mt, uint32_t seed)
{
  int k;

  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }

  mt->arr[0] = (seed * 314159265u + 271828183u) & 0xfffffffful;
  for (k = 1; k < MT_N; k++) { /* the final mask is for 64-bit machines */
    //mt->arr[k] = seed + k;
    mt->arr[k] = 1812433253ul * (mt->arr[k - 1] ^ (mt->arr[k - 1] >> 30)) + k;
    /* mr->arr[k] = (mr->arr[k] + seed) * 22695477ul + 1ul; */
    mt->arr[k] = ((mt->arr[k] + seed) * 314159265ul + 1ul) & 0xfffffffful;
  }
  mt->index = MT_N; /* request for an update */
  mt->inited = 1; /* scrambled */
}



/* save the current state to file */
__inline static int mtrng_save(mtrng_t* mt, const char *fn)
{
  FILE *fp;
  int k;

  if ( !mt->inited ) return 1; /* never used */
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "MTSEED\n%d\n", mt->index);
  for (k = 0; k < MT_N; k++)
    fprintf(fp, "%lu\n", (long unsigned) mt->arr[k]);
  fclose(fp);
  return 0;
}



/* load state from `fn' */
__inline static int mtrng_load(mtrng_t* mt, const char *fn)
{
  char s[64];
  int k, z, err = 1;
  FILE* fp;

  if ( (fp = fopen(fn, "r")) != NULL ) { /* try to load from file */
    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s is empty\n", fn);
    } else if ( strncmp(s, "MTSEED", 6) != 0 ) { /* to check the first line */
      fprintf(stderr, "%s corrupted\n", fn);
    } else if ( fscanf(fp, "%d", &mt->index) != 1 ) {
      fprintf(stderr, "no index in %s\n", fn);
    } else {
      if ( !mt->inited ) mt->index = MT_N; /* request updating */
      for ( z = 1, k = 0; k < MT_N; k++ ) {
        unsigned long val;
        if ( fscanf(fp, "%lu", &val) != 1 ) {
          break;
        }
        mt->arr[k] = (uint32_t) val;
        if ( mt->arr[k] != 0 ) {
          z = 0; /* a non-zero number */
        }
      }
      if ( k != MT_N ) {
        fprintf(stderr, "%s incomplete %d/%d\n", fn, k, MT_N);
      } else {
        err = z; /* clear error, if array is nonzero */
        mt->inited = 1;
      }
    }
    fclose(fp);
  }

  return err;
}



/* load state from `fn' */
__inline static int mtrng_load_or_init_from_seed(mtrng_t* mt, const char *fn, unsigned long seed)
{
  int err = mtrng_load(mt, fn);

  if (err) {
    mtrng_init_from_seed(mt, seed);
  }

  return !mt->inited;
}



__inline static mtrng_t* mtrng_open(uint32_t seed)
{
  mtrng_t* mt;

  if ((mt = (mtrng_t*) calloc(1, sizeof(*mt))) == NULL) {
    exit(1);
  }

  mtrng_init_from_seed(mt, seed);

  return mt;
}


__inline static void mtrng_close(mtrng_t* mt)
{
  free(mt);  
}



/* return an unsigned random number */
__inline static uint32_t mtrng_randuint32(mtrng_t* mt)
{
  static const uint32_t mag01[2] = {0, 0x9908b0dfu}; /* MATRIX_A */
  uint32_t x;
  int k;

  if (!mt->inited) {
    mtrng_init_from_seed(mt, 0);
  }

  if (mt->index >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mt->arr[k] & MT_UMASK) | (mt->arr[k+1] & MT_LMASK);
      mt->arr[k] = mt->arr[k+MT_M] ^ (x>>1) ^ mag01[x&1u];
    }
    for (; k < MT_N-1; k++) {
      x = (mt->arr[k] & MT_UMASK) | (mt->arr[k+1] & MT_LMASK);
      mt->arr[k] = mt->arr[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1u];
    }
    x = (mt->arr[MT_N-1] & MT_UMASK) | (mt->arr[0] & MT_LMASK);
    mt->arr[MT_N-1] = mt->arr[MT_M-1] ^ (x>>1) ^ mag01[x&1u];
    mt->index = 0;
  }
  x = mt->arr[ mt->index++ ];

  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680u;
  x ^= (x << 15) & 0xefc60000u;
  x ^= (x >> 18);

  return x;
}



/* return an unsigned random number */
__inline static double mtrng_rand01(mtrng_t* mt)
{
  return mtrng_randuint32(mt) * (1.0/4294967296.0);
}



/* Gaussian distribution with zero mean and unit variance
 * using the ratio method */
__inline static double mtrng_randgaus(mtrng_t* mt)
{
  double x, y, u, v, q;

  do {
    u = 1 - mtrng_rand01(mt);
    v = 1.7156*(mtrng_rand01(mt) - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));

  return v/u;
}



#endif /* MTRNG_H__ */

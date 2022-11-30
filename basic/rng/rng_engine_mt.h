#ifndef RNG_ENGINE_MT_H__
#define RNG_ENGINE_MT_H__


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
} rng_engine_mt_t;



/* scramble the random number state from a given seed */
__inline static void rng_engine_mt_init_from_seed(rng_engine_mt_t* mt, uint32_t seed)
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
__inline static int rng_engine_mt_save(rng_engine_mt_t* mt, const char *fn)
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
__inline static int rng_engine_mt_load(rng_engine_mt_t* mt, const char *fn)
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
__inline static int rng_engine_mt_load_or_init_from_seed(rng_engine_mt_t* mt, const char *fn, unsigned long seed)
{
  int err = rng_engine_mt_load(mt, fn);

  if (err) {
    rng_engine_mt_init_from_seed(mt, seed);
  }

  return !mt->inited;
}



__inline static rng_engine_mt_t* rng_engine_mt_open(uint32_t seed)
{
  rng_engine_mt_t* mt;

  if ((mt = (rng_engine_mt_t*) calloc(1, sizeof(*mt))) == NULL) {
    exit(1);
  }

  rng_engine_mt_init_from_seed(mt, seed);

  return mt;
}


__inline static void rng_engine_mt_close(rng_engine_mt_t* mt)
{
  free(mt);  
}



/* return an unsigned random number */
__inline static uint32_t rng_engine_mt_randuint32(rng_engine_mt_t* mt)
{
  static const uint32_t mag01[2] = {0, 0x9908b0dfu}; /* MATRIX_A */
  uint32_t x;
  int k;

  if (!mt->inited) {
    rng_engine_mt_init_from_seed(mt, 0);
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
__inline static double rng_engine_mt_rand01(rng_engine_mt_t* mt)
{
  return rng_engine_mt_randuint32(mt) * (1.0/4294967296.0);
}



#endif /* RNG_ENGINE_MT_H__ */

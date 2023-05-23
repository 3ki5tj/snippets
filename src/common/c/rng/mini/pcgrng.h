#ifndef PCGRNG_H__
#define PCGRNG_H__



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>



typedef struct {
  uint64_t state;
  uint64_t inc;
} pcgrng_t;



#ifdef PCGRNG__USE_RNG_ALIASES__

typedef pcgrng_t rng_t;
#define rng__open(seed) pcgrng__open(seed)
#define rng__close(rng) pcgrng__close(rng)
#define rng__save(rng, fn) pcgrng__save(rng, fn)
#define rng__load(rng, fn) pcgrng__load(rng, fn)
#define rng__randuint32(rng) pcgrng__randuint32(rng)
#define rng__rand01(rng) pcgrng__rand01(rng)
#define rng__randgaus(rng) pcgrng__randgaus(rng)

#endif /* defined(PCGRNG_USE_RNG_ALIASES__) */



__inline static void pcgrng__init(pcgrng_t* pcg, int init_default, int state, int inc)
{
  if (init_default) {
    pcg->state = 0x4d595df4d0f33173ULL;
    pcg->inc = 1442695040888963407ULL;
  } else {
    pcg->state = state;
    pcg->inc = inc;
  }
}


__inline static int pcgrng__save(pcgrng_t* pcg, const char* fn)
{
  FILE* fp;
  int err = 1;

  if ((fp = fopen(fn, "w")) != NULL) {
    fprintf(fp, "# PCG %llu %llu\n", (long long unsigned) pcg->state, (long long unsigned) pcg->inc);
    fclose(fp);
    err = 0;
  }

  return err;
}



__inline static int pcgrng__load(pcgrng_t* pcg, const char* fn)
{
  FILE* fp;
  int err = 1;

  if ((fp = fopen(fn, "r")) != NULL) {
    char buf[1024];
    if (fgets(buf, sizeof buf, fp)) {
      long long unsigned state = 0, inc = 0;
      if (2 == sscanf(buf, "# PCG %llu %llu", &state, &inc)) {
        pcg->state = state;
        pcg->inc = inc;
        err = 0;
      }
    }
    fclose(fp);
  }

  return err;
}


__inline static pcgrng_t* pcgrng__open(uint32_t seed)
{
  pcgrng_t* pcg;

  if ((pcg = (pcgrng_t*) calloc(1, sizeof(pcgrng_t))) == NULL) {
    exit(1);
  }

  pcgrng__init(pcg, 1, 0, 0);

  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }

  pcg->state += seed;

  return pcg;
}



__inline static void pcgrng__close(pcgrng_t* pcg)
{
  free(pcg);
}


/* Permuted congruential generator
 * return a 32-bit unsigned random number */
__inline static uint32_t pcgrng__randuint32(pcgrng_t* pcg)
{

  uint64_t old_state = pcg->state;
  pcg->state = old_state * 6364136223846793005ULL + pcg->inc;
  uint32_t xor_shifted = ((old_state >> 18u) ^ old_state) >> 27u;
  uint32_t rot = old_state >> 59u;
  return (xor_shifted >> rot) | (xor_shifted << ((-rot) & 31));
}



__inline static double pcgrng__rand01(pcgrng_t* pcg)
{
  return pcgrng__randuint32(pcg) * (1.0/4294967296.0);
}


/* Gaussian distribution with zero mean and unit variance
 * using the ratio method */
__inline static double pcgrng__randgaus(pcgrng_t* pcg)
{
  double x, y, u, v, q;

  do {
    u = 1 - pcgrng__rand01(pcg);
    v = 1.7156*(pcgrng__rand01(pcg) - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));

  return v/u;
}



#endif /* PCGRNG_H__ */

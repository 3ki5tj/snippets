#ifndef RNG_ENGINE_PCG_H__
#define RNG_ENGINE_PCG_H__



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>



typedef struct {
  uint64_t state;
  uint64_t inc;
} rng_engine_pcg_t;


__inline static void rng_engine_pcg_init(rng_engine_pcg_t* pcg, int init_default, int state, int inc)
{
  if (init_default) {
    pcg->state = 0x4d595df4d0f33173ULL;
    pcg->inc = 1442695040888963407ULL;
  } else {
    pcg->state = state;
    pcg->inc = inc;
  }
}


__inline static int rng_engine_pcg_save(rng_engine_pcg_t* pcg, const char* fn)
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



__inline static int rng_engine_pcg_load(rng_engine_pcg_t* pcg, const char* fn)
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



__inline static rng_engine_pcg_t* rng_engine_pcg_open(uint32_t seed)
{
  rng_engine_pcg_t* pcg;

  if ((pcg = (rng_engine_pcg_t*) calloc(1, sizeof(rng_engine_pcg_t))) == NULL) {
    exit(1);
  }

  rng_engine_pcg_init(pcg, 1, 0, 0);

  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }

  pcg->state += seed;

  return pcg;
}



__inline static void rng_engine_pcg_close(rng_engine_pcg_t* pcg)
{
  free(pcg);
}


/* Permuted congruential generator
 * return a 32-bit unsigned random number */
__inline static uint32_t rng_engine_pcg_randuint32(rng_engine_pcg_t* pcg)
{

  uint64_t old_state = pcg->state;
  pcg->state = old_state * 6364136223846793005ULL + pcg->inc;
  uint32_t xor_shifted = ((old_state >> 18u) ^ old_state) >> 27u;
  uint32_t rot = old_state >> 59u;
  return (xor_shifted >> rot) | (xor_shifted << ((-rot) & 31));
}



__inline static double rng_engine_pcg_rand01(rng_engine_pcg_t* pcg)
{
  return rng_engine_pcg_randuint32(pcg) * (1.0/4294967296.0);
}



#endif /* RNG_ENGINE_PCG_H__ */


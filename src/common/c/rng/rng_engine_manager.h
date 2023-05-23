#ifndef RNG_ENGINE_MANAGER_H__
#define RNG_ENGINE_MANAGER_H__

#include "rng_engine_mt.h"
#include "rng_engine_pcg.h"


enum {
  RNG_ENGINE_TYPE_MT,
  RNG_ENGINE_TYPE_PCG,
  RNG_ENGINE_TYPE_COUNT
};

typedef struct {
  int type;
  rng_engine_mt_t* mt;
  rng_engine_pcg_t* pcg;
} rng_engine_manager_t;



rng_engine_manager_t* rng_engine_manager_open(int type, uint32_t seed)
{
  rng_engine_manager_t* rem;

  if ((rem = (rng_engine_manager_t*) calloc(1, sizeof(*rem))) == NULL) {
    exit(1);
  }

  rem->type = type;

  if (type == RNG_ENGINE_TYPE_MT) {
    rem->mt = rng_engine_mt_open(seed);
  } else if (type == RNG_ENGINE_TYPE_PCG) {
    rem->pcg = rng_engine_pcg_open(seed);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }

  return rem;
}



void rng_engine_manager_close(rng_engine_manager_t* rem)
{
  if (rem->type == RNG_ENGINE_TYPE_MT) {
    rng_engine_mt_close(rem->mt);
  } else if (rem->type == RNG_ENGINE_TYPE_PCG) {
    rng_engine_pcg_close(rem->pcg);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }
  free(rem);
}



int rng_engine_manager_save(rng_engine_manager_t* rem, const char* fn)
{
  if (rem->type == RNG_ENGINE_TYPE_MT) {
    return rng_engine_mt_save(rem->mt, fn);
  } else if (rem->type == RNG_ENGINE_TYPE_PCG) {
    return rng_engine_pcg_save(rem->pcg, fn);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }
  return -1;
}



int rng_engine_manager_load(rng_engine_manager_t* rem, const char* fn)
{
  if (rem->type == RNG_ENGINE_TYPE_MT) {
    return rng_engine_mt_load(rem->mt, fn);
  } else if (rem->type == RNG_ENGINE_TYPE_PCG) {
    return rng_engine_pcg_load(rem->pcg, fn);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }
  return -1;
}



uint32_t rng_engine_manager_randuint32(rng_engine_manager_t* rem)
{
  if (rem->type == RNG_ENGINE_TYPE_MT) {
    return rng_engine_mt_randuint32(rem->mt);
  } else if (rem->type == RNG_ENGINE_TYPE_PCG) {
    return rng_engine_pcg_randuint32(rem->pcg);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }
  return 0;
}




double rng_engine_manager_rand01(rng_engine_manager_t* rem)
{
  if (rem->type == RNG_ENGINE_TYPE_MT) {
    return rng_engine_mt_rand01(rem->mt);
  } else if (rem->type == RNG_ENGINE_TYPE_PCG) {
    return rng_engine_pcg_rand01(rem->pcg);
  } else {
    fprintf(stderr, "Fatal: unknown RNG engine type\n");
    exit(1);
  }
  return 0;
}



#endif /* RNG_ENGINE_MANAGER_H__ */

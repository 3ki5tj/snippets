#include "rng.h"


int main(void)
{
  rng_t* rng = rng_open(RNG_ENGINE_TYPE_MT, 0);

  printf("%f %f\n", rng_rand01(rng), rng_randgaus(rng));

  return 0;
}

#ifndef PCGRAND_H__
#define PCGRAND_H__

#include <stdint.h>

/* Permuted congruential generator
 * return an unsigned random number */
__inline static uint32_t pcgrand(void)
{
  static uint64_t pcg_state_ = 0x4d595df4d0f33173ULL;
  static uint64_t pcg_inc_ = 1442695040888963407ULL;

  uint64_t old_state = pcg_state_;
  pcg_state_ = old_state * 6364136223846793005ULL + pcg_inc_;
  uint32_t xor_shifted = ((old_state >> 18u) ^ old_state) >> 27u;
  uint32_t rot = old_state >> 59u;
  return (xor_shifted >> rot) | (xor_shifted << ((-rot) & 31));
}

#endif


# Object-oriented random number generators

## Overview

pack.py helps pack all files into a single header [rng_pack.h](rng_pack.h)

Aside from [rng_pack.h](rng_pack.h), there are small minimum-style headers
under the [mini](mini) directory

* [pcgrng.h](pcgrng.h)
* [mtrng.h](mtrng.h)

They only offer `_randuint32()`, `_rand01()`, `_randgaus()`.

## Code structure

```ascii
rng_t
  |
  +- rng_engine_manager_t
      |
      +- rng_engine_mt_t_
      |
      +- rng_engine_pcg_t
```

every engine should have the following methods

* `rng_engine_xxx_open(uint32_t seed)`

* `rng_engine_xxx_randuint32()`

* `rng_engine_xxx_rand01()`

* `rng_engine_xxx_close()`

Two RNG engines are available

* mt: Mersenne twister
* pcg: permuted congruential generator

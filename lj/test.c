#include "lj.h"

int main(void)
{
  lj_t *lj;
  lj = lj_open(108, 0.3, 1e9);
  printf("epot %g\n", lj_energy(lj, lj->x, NULL, NULL, NULL));
  lj_writepos(lj, lj->x, lj->v, "a.pos");
  lj_close(lj);
  return 0;
}

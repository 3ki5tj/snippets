#include "lj.h"



int main(void)
{
  int i, n = 108, t, nsteps = 1000;
  double dt = 0.002, dth = dt/2, tp = 1.0;
  lj_t *lj;

  lj = lj_open(n, 0.3, 1e9);
  for ( t = 0; t < nsteps; t++ ) {
    for (i = 0; i < n; i++) { /* VV part 1 */
      vsinc(lj->v[i], lj->f[i], dth);
      vsinc(lj->x[i], lj->v[i], dt);
    }
    lj->epot = lj_force(lj, lj->x, lj->f, NULL, NULL, NULL);
    for (i = 0; i < n; i++) /* VV part 2 */
      vsinc(lj->v[i], lj->f[i], dth);
    //lj->ekin = lj_ekin(lj->v, n);
    lj->ekin = lj_vrescale(lj->v, lj->n, lj->dof, tp, 20*dt);
    printf("%d, ep %g, ek %g, e %g\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
  }
  lj_writepos(lj, lj->x, lj->v, "a.pos");
  lj_close(lj);
  return 0;
}

#include "cagocore.h"


const char *fnpdb = "pdb/1VII.pdb";
double kb = 200.0;
double ka = 40.0;
double kd1 = 1.0;
double kd3 = 0.5;
double nbe = 1.0;
double nbc = 4.0;
double rc = 5.0;

double mddt = 0.002;
double thdt = 0.1;
double tp = 1.0;
long nequil = 1000;
long nsteps = 10000;
long nstrep = 100;



int main(void)
{
  cago_t *go;
  long t;

  if ( (go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", fnpdb);
    return -1;
  }

  cago_initmd(go, 1, 0.01, tp);

  for ( t = 1; t <= nequil + nsteps; t++ ) {
    cago_vv(go, 1.0, mddt);
    //go->ekin = cago_vrescale(go, go->v, tp, thdt);
    go->ekin = cago_ekin(go, go->v);
    printf("%ld: %g\n", t, go->epot + go->ekin);
    if ( t <= nequil ) {
      continue;
    }
    if ( t % nstrep == 0 ) {
      printf("%ld: %g\n", t, go->epot + go->ekin);
    }
  }

  cago_close(go);
  return 0;
}

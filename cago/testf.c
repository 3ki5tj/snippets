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



int main(void)
{
  cago_t *go;
  long t;
  int i;
  double del = 0.001, f2 = 0.0, ep0, ep1;

  if ( (go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", fnpdb);
    return -1;
  }

  cago_initmd(go, 0, 0.1, tp);

  for ( t = 1; t <= 100; t++ ) {
    cago_vv(go, 1.0, mddt);
  }

  ep0 = cago_force(go, go->x, go->f);
  for ( i = 0; i < go->n; i++ ) {
    f2 += vsqr( go->f[i] );
  }
  for ( i = 0; i < go->n; i++ ) {
    vsinc(go->x[i], go->f[i], del / f2 );
  }
  ep1 = cago_force(go, go->x, go->f);
  printf("ep %g -> %g, %g\n", ep0, ep1, (ep0 - ep1)/del);

  cago_close(go);
  return 0;
}

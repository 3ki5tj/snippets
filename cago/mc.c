#include "cagocore.h"



//const char *fnpdb = "pdb/1VII.pdb";
const char *fnpdb= "pdb/1L2Y.pdb";
double kb = 200.0;
double ka = 40.0;
double kd1 = 1.0;
double kd3 = 0.5;
double nbe = 1.0;
double nbc = 4.0;
double rc = 6.0;

double tp = 1.0;
//double tp = 1.1;
double amp = 0.2;
long nequil = 1000000;
long nsteps = 1000000000;
long nstrep = 1000000;



int main(void)
{
  int acc;
  cago_t *go;
  long t, accsm = 0;

  if ( (go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, tp);
  go->dof = go->n * D;

  for ( t = 1; t <= nequil + nsteps; t++ ) {
    acc = cago_metro(go, amp, 1/tp);
    accsm += acc;
    if ( t <= nequil ) {
      continue;
    }
    if ( t % nstrep == 0 ) {
      double rmsd = cago_rmsd(go, go->x, NULL);
      int nc = cago_ncontacts( go, go->x, 1.2, NULL, NULL);
      printf("%ld: acc %.2f%%, ep %g, rmsd %g, nc %d/%d\n",
          t, 100.0*accsm/t, go->epot, rmsd, nc, go->ncont);
      cago_rmcom(go, go->x, go->v);
    }
  }

  cago_close(go);
  return 0;
}

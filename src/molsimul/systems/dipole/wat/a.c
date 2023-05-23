#include "wat.h"

int nw = 798; // 2394 atoms
const char *fnpdb = "wb.pdb";

//int nw = 1;
//const char *fnpdb = "wb_one.pdb";

//int nw = 2;
//const char *fnpdb = "wb_two.pdb";

double l = 30.0;
double alpha = 0.15;
int km = 15;

int main(void)
{
  wat_t *w;

  w = wat_open(nw, l, alpha, km);
  wat_loadpdb(w, fnpdb);
  w->ene = wat_force(w, w->x, w->f, &w->ene_lj, &w->ene_el);
  printf("energy %g, elec %g, LJ %g\n", w->ene, w->ene_el, w->ene_lj);
  wat_savepos(w, "w.pos");
  wat_close(w);
  return 0;
}

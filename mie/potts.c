#include "potts.h"
#include <time.h>

int q = 6;
int n = 15;
double tp = 1.0; // 1.0e40;
int nsys = 100;
long nsteps = 1000000000L;
int nsttraj = 50;
int nstrep = 20000;

int main(void)
{
  potts_t **p;
  int i;
  long t, t1;
  double ent1, ent2, a1, a2;

  mtscramble(clock());
  xnew(p, nsys);
  for ( i = 0; i < nsys; i++ ) {
    p[i] = potts_open(q, n);
  }
  potts_ent2ref(p[0], 1.0/tp);
  a1 = 0.5 * n * (q - 1);
  a2 = 0.25 * n * (n - 1) * (q - 1) * (q - 1);

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      potts_metro(p[i], 1/tp);
      if ( t % nsttraj == 0 ) potts_reg(p[i]);
    }
    if ( t % nstrep == 0 ) {
      ent1 = ent2 = 0;
      for ( i = 0; i < nsys; i++ ) {
        ent1 += potts_ent1(p[i]);
        ent2 += potts_ent2(p[i]);
      }
      ent1 /= nsys;
      ent2 /= nsys;
      t1 = t / nsttraj;
      printf("%9ld: entropy %8.4f(%8.4f), %8.4f(%8.4f); %8.4f(%6.2f), %8.4f(%6.2f);\n", t,
          ent1, p[0]->ent1r, ent2, p[0]->ent2r,
          (p[0]->ent1r - ent1)*t1, a1,
          (p[0]->ent2r - ent2)*t1, a2);
    }
  }
  for ( i = 0; i < nsys; i++ ) {
    potts_close(p[i]);
  }
  free(p);
  return 0;
}

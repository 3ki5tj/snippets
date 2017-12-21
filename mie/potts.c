#include "potts.h"
#include "argopt.h"

int q = 6;
int n = 15;
double tp = 1.0; // 1.0e40;
int nsys = 100;
long nsteps = 1000000000L;
long nsttraj = 50;
long nstrep = 20000;

static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "Mutual information expansion for 1D Potts model";
  argopt_add(ao, "-q", "%d", &q, "number of states in each spin");
  argopt_add(ao, "-n", "%d", &n, "number of spins");
  argopt_add(ao, "-T", "%lf", &tp, "temperature");
  argopt_add(ao, "-M", "%d", &nsys, "number of copies");
  argopt_add(ao, "-t", "%ld", &nsteps, "number of steps");
  argopt_add(ao, "-j", "%ld", &nsttraj, "number of steps to deposit state");
  argopt_add(ao, "-r", "%ld", &nstrep, "number of steps to report");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

int main(int argc, char **argv)
{
  potts_t **p;
  int i;
  long t, t1;
  double ent1, ent2, a1, a2;

  doargs(argc, argv);
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

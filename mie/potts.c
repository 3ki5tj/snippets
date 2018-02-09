#include "potts.h"
#include "argopt.h"

int q = 6;
int n = 10;
double tp = 0.5; // 1.0e40;
int nsys = 1000;
long nsteps = 100000L;
long nsttraj = 10;
long nstrep = 1000;
int npart = 2; /* number of partitions for the block method */

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
  argopt_add(ao, "-P", "%d", &npart, "number of partitions for the block method");
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
  double ent1, ent1c, ent2, ent2c, a1, a2;

  doargs(argc, argv);
  mtscramble(clock());
  xnew(p, nsys);
  fprintf(stderr, "estimated memory use %gM\n",
      nsteps/nsttraj*nsys*n/(1024.*1024));
  for ( i = 0; i < nsys; i++ ) {
    p[i] = potts_open(q, n, nsteps / nsttraj);
    if ( p[i] == NULL ) return -1;
  }
  potts_ent2ref(p[0], 1.0/tp);
  a1 = 0.5 * n * (q - 1);
  a2 = 0.25 * n * (n - 1) * (q - 1) * (q - 1);

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      potts_metro(p[i], 1/tp);
      if ( t % nsttraj == 0 ) potts_add(p[i]);
    }

    if ( t % nstrep == 0 ) {
      ent1 = ent1c = 0;
      ent2 = ent2c = 0;
      for ( i = 0; i < nsys; i++ ) {
        potts_entropy(p[i], npart);
        ent1 += p[i]->ent1;
        ent2 += p[i]->ent2;
        ent1c += p[i]->ent1c;
        ent2c += p[i]->ent2c;
      }
      ent1 /= nsys;
      ent2 /= nsys;
      ent1c /= nsys;
      ent2c /= nsys;
      t1 = t / nsttraj;
      printf("%9ld: entropy %8.4f,%8.4f(%8.4f), %8.4f,%8.4f(%8.4f),"
             " (%8.4f); %8.4f(%6.2f), %8.4f(%6.2f);\n", t,
          ent1, ent1c, p[0]->ent1r,
          ent2, ent2c, p[0]->ent2r,
          p[0]->entr,
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

#include "potts.h"
#include "argopt.h"
#include "ave.h"

int q = 6;
int n = 10;
double tp = 0.5; // 1.0e40;
int nsys = 1000;
long nsteps = 100000L;
long nsttraj = 10;
long nstrep = 1000;
int npart = 2; /* number of partitions for the block method */
char *fnlog = "potts.log";
int blkave = 0;


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
  argopt_add(ao, "-K", "%b", &blkave, "averaging over blocks for the block method");
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
  double a1, a2;
  av_t avent1[1], avent1c[1], avent2[1], avent2c[1];
  double s1, s1c, s2, s2c;
  double g1, g1c, g2, g2c;
  FILE *fplog;

  doargs(argc, argv);
  mtscramble(time(NULL) + clock());
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
  if ( (fplog = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open the log file %s\n", fnlog);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      potts_metro(p[i], 1/tp);
      if ( t % nsttraj == 0 ) potts_add(p[i]);
    }

    if ( t % nstrep == 0 ) {
      av_clear(avent1);
      av_clear(avent1c);
      av_clear(avent2);
      av_clear(avent2c);
      for ( i = 0; i < nsys; i++ ) {
        potts_entropy(p[i], npart, blkave);
        av_add(avent1,  p[i]->ent1);
        av_add(avent1c, p[i]->ent1c);
        av_add(avent2,  p[i]->ent2);
        av_add(avent2c, p[i]->ent2c);
      }
      s1  = av_getave(avent1,  &g1);
      s1c = av_getave(avent1c, &g1c);
      s2  = av_getave(avent2,  &g2);
      s2c = av_getave(avent2c, &g2c);
      t1 = t / nsttraj;
      printf("%9ld: entropy %8.4f,%8.4f(%8.4f), %8.4f,%8.4f(%8.4f),"
             " (%8.4f); %8.4f(%6.2f), %8.4f(%6.2f);\n", t,
          s1, s1c, p[0]->ent1r,
          s2, s2c, p[0]->ent2r,
          p[0]->entr,
          (p[0]->ent1r - s1)*t1, a1,
          (p[0]->ent2r - s2)*t1, a2);
      fprintf(fplog, "%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t,
          s1, sqrt(g1), s1c, sqrt(g1c), p[0]->ent1r,
          s2, sqrt(g2), s2c, sqrt(g2c), p[0]->ent2r);
    }
  }
  for ( i = 0; i < nsys; i++ ) {
    potts_close(p[i]);
  }
  free(p);
  fclose(fplog);
  return 0;
}

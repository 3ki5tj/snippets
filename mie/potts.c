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
  int i, k;
  long t;
  //double a1, a2;
  av_t avent1[STOT], avent2[STOT], avent3[STOT];
  double s1[STOT], s2[STOT], s3[STOT]; /* average entropies */
  double g1[STOT], g2[STOT], g3[STOT]; /* variances */
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
  potts_entref(p[0], 1.0/tp);
  //a1 = 0.5 * n * (q - 1);
  //a2 = 0.25 * n * (n - 1) * (q - 1) * (q - 1);
  if ( (fplog = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open the log file %s\n", fnlog);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      potts_metro(p[i], 1/tp);
      if ( t % nsttraj == 0 ) potts_add(p[i]);
    }

    if ( t % nstrep == 0 ) {
      for ( k = 0; k < STOT; k++ ) {
        av_clear(&avent1[k]); /* first-order */
        av_clear(&avent2[k]); /* second-order */
        av_clear(&avent3[k]); /* third-order */
      }
      for ( i = 0; i < nsys; i++ ) {
        potts_entropy(p[i], npart);
        for ( k = 0; k < STOT; k++ ) {
          av_add(&avent1[k], p[i]->ent1[k]); /* first-order */
          av_add(&avent2[k], p[i]->ent2[k]); /* second-order */
          av_add(&avent3[k], p[i]->ent3[k]); /* third-order */
        }
      }
      for ( k = 0; k < STOT; k++ ) {
        s1[k] = av_getave(&avent1[k], &g1[k]); /* first-order */
        s2[k] = av_getave(&avent2[k], &g2[k]); /* second-order */
        s3[k] = av_getave(&avent3[k], &g3[k]); /* third-order */
      }
      //t1 = t / nsttraj;
      printf("%9ld: entropy %6.2f,%6.2f,%6.2f,%6.2f(%6.2f), "
             "%6.2f,%6.2f,%6.2f,%6.2f(%6.2f), "
             "%6.2f,%6.2f,%6.2f,%6.2f(%6.2f), "
             "(%8.4f)\n", t,
          s1[0], s1[SBAV], s1[SLIN], s1[SEXP], p[0]->ent1r,
          s2[0], s2[SBAV], s2[SLIN], s2[SEXP], p[0]->ent2r,
          s3[0], s3[SBAV], s3[SLIN], s3[SEXP], p[0]->ent3r,
          p[0]->entr);
      fprintf(fplog, "%ld\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"
                          "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"
                          "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", t,
          s1[0], sqrt(g1[0]), s1[SLIN], sqrt(g1[SLIN]), s1[SEXP], sqrt(g1[SEXP]), p[0]->ent1r,
          s2[0], sqrt(g2[0]), s2[SLIN], sqrt(g2[SLIN]), s2[SEXP], sqrt(g2[SEXP]), p[0]->ent2r,
          s3[0], sqrt(g3[0]), s3[SLIN], sqrt(g3[SLIN]), s3[SEXP], sqrt(g3[SEXP]), p[0]->ent3r);
      fflush(fplog);
    }
  }
  for ( i = 0; i < nsys; i++ ) {
    potts_close(p[i]);
  }
  free(p);
  fclose(fplog);
  return 0;
}

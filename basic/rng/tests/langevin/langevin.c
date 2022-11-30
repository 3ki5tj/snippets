#include "rng.h"


const double move_size = 0.1;


#define HIST_N 100
double xmin = 1;
double xmax = 2;
double dx = 0.01;
double hist[HIST_N] = {0};



void save_hist(int nsteps, int apply_correction)
{
  int i;
  char fn[FILENAME_MAX];
  double norm;
  FILE* fp;

  sprintf(fn, "hist%d.dat", apply_correction);

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "Cannot write file %s\n", fn);
    return;
  }

  norm = nsteps * dx;

  for ( i = 0; i < HIST_N; i++ ) {
    fprintf(fp, "%g\t%g\t%g\n", xmin+(i+0.5)*dx, hist[i], hist[i]/norm);
  }

  fclose(fp);
}



void sample(int nsteps, int apply_correction)
{
  int t = 0;
  double x = 1.5, nx, r, nr, xp;
  rng_t* rng = rng_open(RNG_ENGINE_TYPE_PCG, 0);

  fprintf(stderr, "Sampling %d steps, with correction %d...\n", nsteps, apply_correction);

  for (t = 0; t < nsteps; t++) {
    r = rng_randgaus(rng);
    nx = x + x * r * move_size;
    if ( nx >= xmin && nx < xmax ) {

      if (apply_correction) {
        nr = x * r / nx;
        xp = (r*r - nr*nr)/2 + log(x/nx);
        if ( xp > 0 || rng_rand01(rng) < exp(xp) ) {
          x = nx;
        }
      } else {
        x = nx;
      }

    }

    {
      int i = (int)((x - xmin)/dx);
      hist[i] += 1;
    }
  }

  save_hist(nsteps, apply_correction);

  rng_close(rng);
}



int main(int argc, char** argv)
{
  int nsteps = 100000000;
  int apply_correction = 0;

  if (argc > 1) {
    apply_correction = atoi(argv[1]);
  }

  sample(nsteps, apply_correction);

  return 0;
}

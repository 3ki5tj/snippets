/* molecular dynamics of a two-degrees-of-freedom system
 * measure the distribution of x
 * this model demonstrates the infinite stiffness limit
 * is not uniquely defined */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtrand.h"

/* parameters */
int thermostat = 1; /* whether to use thermostat */
double dt = 0.01; /* MD time step */
double kappa = 1; /* spring stiffness */
double sqrtk = 0; /* to be set to sqrt(kappa) */
double alpha = 2.0; /* x-y coupling constant */
int nsteps = 10000000;
int nequil = 10000;
/* thermostat variables */
double tp = 1.0;
double thermdt = 0.01;

/* MD variables */
#define N 2 /* number of particles */
double x[] = {1.0, 0};
double v[] = {1.0, 1.0};
double f[] = {0, 0};
double ep, ek;

/* histogram variables */
#define HIST_XMAX   10.0
#define HIST_XMIN   (-HIST_XMAX)
#define HIST_DX     0.02
#define HIST_N      (int)(((HIST_XMAX - HIST_XMIN)/HIST_DX) + 1)
double hist[2][HIST_N];
const char *fnhist = "xy.his";

/* compute the force, return the potential energy
 * the potential is
 * V(x, y) = 0.5 x^2 + 0.5 kappa y^2 (1 + alpha x^2)
 * */
static double force(void)
{
  double ep = 0, x2, y2;

  x2 = x[0] * x[0];
  y2 = x[1] * x[1];
  ep = 0.5 * x2 + 0.5 * kappa * y2 * (1 + alpha * x2);
  f[0] = -x[0] * (1 + kappa * alpha * y2);
  f[1] = -kappa * x[1] * (1 + alpha * x2);

  return 0;
}


/* velocity rescaling thermostat */
static double vrescale(double temp, double thdt, int dof)
{
  int i;
  double ek1, ek2, s, c, r, r2;

  /* normal velocity rescaling */
  for ( ek1 = 0, i = 0; i < N; i++ )
    ek1 += 0.5 * v[i] * v[i];

  c = (thdt < 700) ? exp(-thdt) : 0;
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = c * ek1 + (1 - c) * (r2 + r*r) * .5 * temp
      + r * sqrt(c * (1 - c) * 2 * temp * ek1);

  if (ek2 < 1e-30) ek2 = 1e-30;
  s = sqrt(ek2/ek1);
  for (i = 0; i < N; i++) v[i] *= s;

  return ek2;
}



static void runmd(void)
{
  long t;
  int i;
  double ys; /* scaled y */

  for ( t = 1; t <= nsteps + nequil; t++ ) {
    /* basic velocity verlet */
    for ( i = 0; i < N; i++ ) {
      v[i] += f[i] * 0.5 * dt;
      x[i] += v[i] * dt;
    }
    ep = force();
    for ( i = 0; i < N; i++ ) {
      v[i] += f[i] * 0.5 * dt;
    }

    /* apply the thermostat */
    if ( thermostat ) {
      ek = vrescale(tp, thermdt, N);
    } else {
      for ( ek = 0, i = 0; i < N; i++ ) {
        ek += 0.5 * v[i] * v[i];
      }
    }

    if ( t < nequil ) continue;
    /* accumulate histogram */
    if ( x[0] > HIST_XMIN && x[0] < HIST_XMAX ) {
      i = (int) ( (x[0] - HIST_XMIN) / HIST_DX);
      hist[0][i] += 1;
    }
    ys = x[1] * sqrtk;
    if ( ys > HIST_XMIN && ys < HIST_XMAX ) {
      i = (int) ( ( ys - HIST_XMIN) / HIST_DX );
      hist[1][i] += 1;
    }
  }
}

/* plot the x-distribution from histogram */
static int mkplot(const char *fn)
{
  int i, imin = HIST_N, imax = 0;
  double tot[2] = {0, 0};
  FILE *fp;

  /* search the range of x */
  for ( i = 0; i < HIST_N; i++ ) {
    if ( hist[0][i] > 0 || hist[1][i] > 0 ) {
      if ( i > imax ) imax = i;
      if ( i < imin ) imin = i;
      tot[0] += hist[0][i];
      tot[1] += hist[1][i];
    }
  }

  /* write the file */
  if ((fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %g\n", HIST_DX);
  for ( i = imin; i <= imax; i++ ) {
    fprintf(fp, "%g %g %g %g %g\n",
        HIST_XMIN + (i + 0.5) * HIST_DX,
        hist[0][i] / (tot[0] * HIST_DX),
        hist[1][i] / (tot[1] * HIST_DX),
        hist[0][i], hist[1][i]);
  }
  fclose(fp);

  return -1;
}

int main(void)
{
  sqrtk = sqrt(kappa);
  runmd();
  mkplot(fnhist);
  return 0;
}

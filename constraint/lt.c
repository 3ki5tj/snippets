/* molecular dynamics of a two-degrees-of-freedom system
 * elastic string with an angle constraint
 * measure the distribution of l
 * with or without angular constraint */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtrand.h"

/* parameters */
int constraint = 0; /* whether to apply the constraint */
int thermostat = 0; /* whether to use thermostat */
double dt = 0.01; /* MD time step */
double alpha = 16.0; /* stiffness for length */
double kappa = 100; /* spring stiffness for the theta */
double sqrtk = 0; /* to be set to sqrt(kappa) */
int nsteps = 10000000;
int nequil = 10000;
/* thermostat parameters */
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
const char *fnhist[2] = {"l.his", "th.his"};



/* compute the force, return the potential energy
 * the potential is
 *   V(x, y) = 0.5 alpha (l - l0)^2 - kappa cos(theta)
 * with l0 = 1;
 * */
static double force(void)
{
  double ep = 0, l2, l, dl, c;

  f[0] = f[1] = 0;

  /* 1. length force
   * Vl(x, y) = 0.5 alpha (l - l0)^2
   * fx = -dVl/dx = -alpha (l - l0) x / l
   * fy = -dVl/dy = -alpha (l - l0) y / l */
  l2 = x[0] * x[0] + x[1] * x[1];
  l = sqrt(l2);
  dl = l - 1;
  f[0] += -alpha * dl * x[0] / l;
  f[1] += -alpha * dl * x[1] / l;
  ep += 0.5 * alpha * dl * dl;

  /* 2. angular force
   * Va(x, y) = -k cos(theta) = -k x / l = -k x / sqrt(x^2 + y^2)
   * fx = -dVa/dx = k y^2 / l^3
   * fy = -dVa/dx = -k x y / l^3 */
  c = x[0] / l;
  f[0] +=  kappa * x[1] * x[1] / (l * l2);
  f[1] += -kappa * x[0] * x[1] / (l * l2);
  ep += -kappa * c;
  return ep;
}



/* apply the angular constraint */
static void constrain(void)
{
  double l2, l, vl;

  /* compute the length */
  l2 = x[0] * x[0] + x[1] * x[1];
  l = sqrt(l2);

  /* compute the parallel velocity */
  vl = v[0] * x[0] + v[1] * x[1];

  x[0] = l;
  x[1] = 0;
  v[0] = vl / l;
  v[1] = 0;
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
  int i, dof;
  double l, sth; /* scaled, theta */

  /* degrees of freedom for the thermostat */
  dof = ( constraint ? N - 1 : N );

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

    /* apply the constraint */
    if ( constraint ) {
      constrain();
    }

    /* apply the thermostat */
    if ( thermostat ) {
      ek = vrescale(tp, thermdt, dof);
    } else {
      for ( ek = 0, i = 0; i < N; i++ ) {
        ek += 0.5 * v[i] * v[i];
      }
    }

    if ( t < nequil ) continue;
    /* accumulate histogram */
    l = sqrt(x[0] * x[0] + x[1] * x[1]);
    if ( l > HIST_XMIN && l < HIST_XMAX ) {
      i = (int) ( (l - HIST_XMIN) / HIST_DX);
      hist[0][i] += 1;
    }

    if ( !constraint ) {
      sth = sqrtk * atan2(x[1], x[0]);
      if ( sth > HIST_XMIN && sth < HIST_XMAX ) {
        i = (int) ( ( sth - HIST_XMIN) / HIST_DX );
        hist[1][i] += 1;
      }
    }
  }
}

/* plot the x-distribution from histogram */
static int mkplot(const double *hist, const char *fn)
{
  int i, imin = HIST_N, imax = 0;
  double tot = 0;
  FILE *fp;

  /* search the range of x */
  for ( i = 0; i < HIST_N; i++ ) {
    if ( hist[i] > 0 || hist[i] > 0 ) {
      if ( i > imax ) imax = i;
      if ( i < imin ) imin = i;
      tot += hist[i];
    }
  }

  /* write the file */
  if ((fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %g\n", HIST_DX);
  for ( i = imin; i <= imax; i++ ) {
    fprintf(fp, "%g %g %g\n",
        HIST_XMIN + (i + 0.5) * HIST_DX,
        hist[i] / (tot * HIST_DX), hist[i]);
  }
  fclose(fp);

  return -1;
}

int main(void)
{
  sqrtk = sqrt(kappa);
  runmd();
  mkplot(hist[0], fnhist[0]);
  if ( !constraint ) {
    mkplot(hist[1], fnhist[1]);
  }
  return 0;
}

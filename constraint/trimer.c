/* molecular dynamics of the trimer model by Frenkel
 * with or without constraints */
#define D 3
#include "mtrand.h"
#include "vct.h"
#include "mdutil.h"

/* parameters */
int constraint = 0; /* whether to apply the constraint */
int thermostat = 1; /* whether to use thermostat */
int norot = 0; /* annihilate rotation */
double dt = 0.01; /* MD time step */
double kappa = 100; /* spring stiffness for the length */
int nsteps = 100000000;
int nequil = 10000;
/* thermostat parameters */
double tp = 1.0;
double thermdt = 0.01;

/* MD variables */
#define N 3 /* number of particles */
double x[][D] = {{0, 0}, {1, 0}, {0, 1}};
double v[][D] = {{0, 0}, {0, 1}, {1, 0}};
double f[][D] = {{0, 0}, {0, 0}, {0, 0}};
double ep, ek;

/* histogram variables */
#define HIST_XMAX   M_PI
#define HIST_XMIN   (-HIST_XMAX)
#define HIST_DX     (HIST_XMAX/60)
#define HIST_N      ((int)(((HIST_XMAX - HIST_XMIN)/HIST_DX) + 1))
double hist[HIST_N] = {0};
const char *fnhist = "ang.his";



/* remove center of mass motion and rotations */
static void trimer_rmcom(int norotation)
{
  rmcom(x, NULL, N);
  if ( norotation ) {
    rmcom(v, NULL, N);
    shiftang(x, v, NULL, N);
  }
}



/* bond energy 1/2 k (r - r0)^2 */
static double potbond(double *a, double *b,
    double r0, double kr, double *fa, double *fb)
{
  double dx[D], r, dr, amp;

  r = vnorm( vdiff(dx, a, b) );
  dr = r - r0;
  amp = kr * dr / r;
  vsinc(fa, dx, -amp);
  vsinc(fb, dx,  amp);
  return 0.5 * kr * dr * dr;
}



/* compute the force, return the potential energy
 * the potential is the two bond lengths */
static double force(void)
{
  double ep = 0;
  int i;

  for ( i = 0; i < N; i++ ) {
    vzero( f[i] );
  }

  if ( constraint ) return ep;

  /* bonds 0-1 and 0-2 */
  for ( i = 1; i < N; i++ ) {
    ep += potbond(x[i], x[0], 1.0, kappa, f[i], f[0]);
  }

  return ep;
}



/* apply the constraints with SHAKE and RATTTLE */
static void constrain(int itmax, double tol)
{
  int i, it;
  double vl, dev, maxdev, dot;
  double dx0[N][D], dx[D], dv[D];

  for ( i = 1; i < N; i++ ) {
    vdiff(dx0[i], x[i], x[0]);
  }

  /* SHAKE */
  for ( it = 0; it < itmax; it++ ) {
    maxdev = 0;
    for ( i = 1; i < N; i++ ) {
      /* constraint the distance between atoms 0 and i */
      vdiff(dx, x[i], x[0]);
      dev = vsqr(dx) - 1;
      if ( dev > maxdev ) maxdev = dev;
      dot = vdot(dx, dx0[i]);
      if ( dot < 0.2 ) { /* really bad */
        dot = 0.2;
      }
      dev = dev / (2 * dot);
      vsinc(x[i], dx0[i], -dev * 0.5);
      vsinc(x[0], dx0[i],  dev * 0.5);
    }
    if ( maxdev < tol ) break;
  }

  /* RATTLE */
  for ( it = 0; it < itmax; it++ ) {
    maxdev = 0;
    for ( i = 1; i < N; i++ ) {
      /* annihilate the parallel velocity */
      vdiff(dx, x[i], x[0]);
      vdiff(dv, v[i], v[0]);
      /* the relative velocity should be perpendicular
       * to dr, if not, subtract the nonperpendicular one */
      vl = vdot(dv, dx);
      /* because the reference length is 1.0
       * we do not have to divide it */
      vsinc(v[i], dx, -vl * 0.5);
      vsinc(v[0], dx, +vl * 0.5);
      dev = fabs(vl);
      if ( dev > maxdev ) {
        maxdev = dev;
      }
    }
    if ( maxdev < tol ) break;
  }
}



static void runmd(void)
{
  long t;
  int i, j, dof;
  double ang;

  /* degrees of freedom for the thermostat */
  if ( norot ) {
    dof = ( constraint ? 1 : 3);
  } else {
    dof = constraint ? D * N - 2 : D * N;
  }

  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < D; j++ ) {
      v[i][j] = randgaus();
    }
    #if D == 3
    /* annihilate the last dimension velocities */
    if ( norot ) {
      v[i][2] = 0;
    }
    #endif
  }
  trimer_rmcom(norot);

  for ( t = 1; t <= nsteps + nequil; t++ ) {
    /* basic velocity verlet */
    for ( i = 0; i < N; i++ ) {
      vsinc(v[i], f[i], 0.5 * dt);
      vsinc(x[i], v[i], dt);
    }
    if ( !constraint ) {
      ep = force();
    }
    for ( i = 0; i < N; i++ ) {
      vsinc(v[i], f[i], 0.5 * dt);
    }

    ang = vang(x[1], x[0], x[2], NULL, NULL, NULL);
    if ( t % 1000000 == 0 ) {
      printf("%g %g %g %g %g %g %g %g %g, ang %g, ep %g\n", x[0][0], x[0][1], x[0][2], x[1][0], x[1][1], x[1][2], x[2][0], x[2][1], x[2][2], ang, ep);
      // getchar();
    }
    /* apply the constraint */
    if ( constraint ) {
      constrain(1000, 1e-6);
    }

    /* apply the thermostat */
    if ( thermostat ) {
      if ( !norot && t % 20 == 0 ) {
        /* Andersen thermostat */
        i = (int) (rand01() * N);
        for ( j = 0; j < D; j++ ) {
          v[i][j] = randgaus() * sqrt(tp);
        }
        if ( constraint ) {
          constrain(1000, 1e-6);
        }
        /* constraint? */
      } else {
        ek = md_vrescale(v, NULL, N, dof, tp, thermdt);
      }
    } else {
      for ( ek = 0, i = 0; i < N; i++ ) {
        ek += 0.5 * vsqr( v[i] );
      }
    }

    if ( t % 10000 == 0 ) {
      trimer_rmcom(norot);
    }

    if ( t < nequil ) continue;
    /* accumulate histogram of the angle */
    if ( ang > HIST_XMIN && ang < HIST_XMAX ) {
      i = (int) ( (ang - HIST_XMIN) / HIST_DX);
      hist[i] += 1;
    }
  }
}

/* plot the x-distribution from histogram */
static int mkplot(void)
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
  if ((fp = fopen(fnhist, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fnhist);
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
  runmd();
  mkplot();
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mtrand.h" /* random number generator */
#include "eig.h" /* eigenvalues */
#include "lu.h" /* LU decomposition */



/* To change the number of oscillators, define N in compiling
 *    icc -DN=4 ho.c
 * */
#ifndef N
#define N 3
#endif

#ifndef THERMOSTAT
#define THERMOSTAT LANGEVIN
#endif

#ifndef TRIGMASS
#define TRIGMASS 0.0
#endif



enum { ANDERSEN, LANGEVIN, VRESCALE, RANDOM_COLLISION_MC, RANDOM_COLLISION_LANGEVIN };



double hess[N];
double mass[N];
double x[N];
double v[N];
double f[N];
double trigmass = TRIGMASS; /* triangular mass distribution */
double temp = 1; /* temperature */
double dt = 0.005; /* time step for molecular dynamics */
double tstatdt = 0.01; /* thermostat time step */
int thermostat = THERMOSTAT;

double quartic_a = 0.0;
double quartic_b = 1.0;

double fsmallworld = 0.0; /* frequency of adding small-world springs */
double ksmallworld = 1.0; /* stiffness of small-world springs */
int smallworld[N][N];

long long nsteps = 100000000; /* number of molecular dynamics steps */
long long nstblah = 1000000; /* frequency of reporting */
long long nstequiv = 10000; /* number of steps for equilibration */

int nstsamp = 10; /* frequency of sampling */
long long nsamps; /* number of samples collected so far */
double qsum[N];
double qqsum[N][N]; /* correlations */



/* initialize the force field */
static void initff(void)
{
  int i, j;

  for ( i = 0; i < N; i++ ) {
    hess[i] = 1;
  }

  /* add small-world springs */
  for ( i = 0; i < N; i++ )
    for ( j = i + 2; j < N; j++ ) {
      smallworld[i][j] = ( rand01() < fsmallworld );
      if ( smallworld[i][j] )
        printf("adding small-world springs %d - %d\n", i, j);
    }

  for ( i = 0; i < N; i++ ) {
    /* triangular mass distribution,
     * trigmass > 0, increase from 1 to 1 + trigmass, then decrease back to 1
     * trigmass < 0, decrease from 1 + trigmass to 1, then increase back to trigmass */
    if ( fabs(trigmass) > DBL_MIN && N > 1 ) {
      mass[i] = -fabs(2.*i/(N-1) - 1)*trigmass + 1;
      if ( trigmass > 0 ) mass[i] += trigmass;
    } else {
      mass[i] = 1;
    }
    printf("mass %d: %g\n", i, mass[i]);
  }
}



static double getpairpot(double x, double k, double a, double b, double *df)
{
  double dx2 = x*x - b*b;
  *df = k * x + a * dx2 * x;
  return 0.5 * k * x * x + 0.25 * a * dx2 * dx2;
}



/* compute the force, return the energy */
static double force(void)
{
  int i, j;
  double df, U = 0;

  for ( i = 0; i < N; i++ ) f[i] = 0;

  if ( N == 1 ) {
    U = getpairpot(-x[0], hess[0], quartic_a, quartic_b, &f[0]);
  } else {
    for ( i = 0; i < N - 1; i++ ) { /* loop over springs */
      U += getpairpot(x[i+1] - x[i], hess[i], quartic_a, quartic_b, &df);
      f[i]   += df;
      f[i+1] -= df;
    }

    // force from small-world springs
    if ( fsmallworld > 0 ) {
      for ( i = 0; i < N; i++ ) {
        for ( j = i + 2; j < N; j++ ) {
          if ( smallworld[i][j] ) {
            U += getpairpot(x[j] - x[i], ksmallworld, quartic_a, quartic_b, &df);
            f[i] += df;
            f[j] -= df;
          }
        }
      }
    }
  }
  return U;
}



/* remove the center of mass motion */
static void rmcom(double *arr)
{
  int i;
  double p = 0, mt = 0;

  if ( N <= 1 ) return;
  for ( i = 0; i < N; i++ ) {
    p += arr[i] * mass[i];
    mt += mass[i];
  }
  p /= mt;
  for ( i = 0; i < N; i++ )
    arr[i] -= p;
}



static double getEk(double *vel)
{
  int i;
  double ek = 0;

  for ( i = 0; i < N; i++ )
    ek += .5 * mass[i] * vel[i] * vel[i];
  return ek;
}



/* Langevin-like thermostat */
static double langevin(double thdt)
{
  int i;

  for ( i = 0; i < N; i++ )
    v[i] += (-thdt * v[i] + sqrt(2*temp*thdt) * gaussrand()) / mass[i];
  rmcom(v);
  return getEk(v);
}



/* Andersen thermostat
 * Note, DOF should be N in this case */
static double andersen(void)
{
  int i;
  double ek, r;

  if ( rand01() < 0.1 ) {
    i = (int) (N * rand01());
    r = gaussrand();
    v[i] = sqrt(temp/mass[i]) * r;
  }
  for (ek = 0, i = 0; i < N; i++)
    ek += .5 * mass[i] * v[i] * v[i];
  return ek;
}



/* return the sum of the squares of n Gaussian random numbers  */
__inline static double randgausssum(int n)
{
  double x = 0., r;
  if (n > 0) {
    x = 2.0 * randgam(n/2);
    if (n % 2) { r = gaussrand(); x += r*r; }
  }
  return x;
}



/* velocity rescaling thermostat
 * do not use, not ergodic! */
static double vrescale(double thdt, int dof)
{
  int i;
  double ekav = .5f*temp*dof, ek1, ek2, s, c, r, r2;

  /* normal velocity rescaling */
  ek1 = getEk(v);

  /* approximate algorithm of computing ek2
   * only valid for a small thdt <= 0.001 */
/*
  c = sqrt(2*ek1*temp*thdt);
  ek2 = ek1 + (ekav - ek1)*thdt + c*gaussrand();
*/
  c = (thdt < 700) ? exp(-thdt) : 0;
  r = gaussrand();
  r2 = randgausssum(dof - 1);
  ek2 = c * ek1 + (1 - c) * (r2 + r*r) * .5 * temp
      + r * sqrt(c * (1 - c) * 2 * temp * ek1);

  if (ek2 < 1e-30) ek2 = 1e-30;
  s = sqrt(ek2/ek1);
  for (i = 0; i < N; i++) v[i] *= s;

  return ek2;
}



/* randomly collide two particles */
static double random_collision_mc(void)
{
  int i, j;
  double m, frac, vi, vj, dv, r, amp;

  i = (int) (N * rand01());
  j = ((int) ((N - 1) * rand01()) + i + 1) % N;

  frac = mass[i] / (mass[i] + mass[j]);
  m = mass[j] * frac;

  amp = sqrt(2 * temp / m);
  dv = (rand01()*2 - 1) * amp;
  /* distribute dv to i and j such that
   * mass[i] * v[i] + mass[j] * v[j] is conserved */
  vi = v[i] + dv * (1 - frac);
  vj = v[j] - dv * frac;
  r = .5 * mass[i] * (vi*vi - v[i]*v[i])
    + .5 * mass[j] * (vj*vj - v[j]*v[j]);
  if ( r < 0 || rand01() < exp(-r/temp) ) {
    v[i] = vi;
    v[j] = vj;
  }
  return getEk(v);
}




/* randomly collide two particles */
static double random_collision_langevin(double thdt)
{
  int i, j;
  double m, frac, vij, dv;

  i = (int) (N * rand01());
  j = ((int) ((N - 1) * rand01()) + i + 1) % N;

  frac = mass[i] / (mass[i] + mass[j]);
  m = mass[j] * frac;
  vij = v[i] - v[j];

  /* do a step of Langevin equation */
  dv = (-thdt * vij + sqrt(2*temp*thdt) * gaussrand()) / m;
  /* distribute dv to i and j such that
   * mass[i] * v[i] + mass[j] * v[j] is conserved */
  v[i] += dv * (1 - frac);
  v[j] -= dv * frac;
  return getEk(v);
}



/* collect position correlations */
static void sample()
{
  int i, j;

  for ( i = 0; i < N; i++ ) {
    qsum[i] += x[i];
    for ( j = i; j < N; j++ )
      qqsum[i][j] += x[i] * x[j];
  }
  nsamps += 1;
}



static void domd(void)
{
  long long t;
  double Ep, Ek, smT = 0, smU = 0, r;
  int i = 0, dof;

  /* initialize the random velocities */
  for ( Ek = 0, i = 0; i < N; i++ ) {
    x[i] = 0;
    r = gaussrand();
    v[i] = sqrt(temp/mass[i]) * r;
    Ek += .5 * mass[i] * v[i] * v[i];
  }
  rmcom(x);
  rmcom(v);

  dof = ( thermostat == ANDERSEN || N == 1 ) ? N : N - 1;
  Ep = force();
  for ( t = 1; t <= nsteps; t += 1 ) {
    /* integrate Newton's equation of motion */
    for ( i = 0; i < N; i++ ) {
      v[i] += f[i]/mass[i] * .5 * dt;
      x[i] += v[i] * dt;
    }

    Ep = force();

    for ( i = 0; i < N; i++ ) {
      v[i] += f[i]/mass[i] * .5 * dt;
    }

    /* apply the thermostat */
    if ( thermostat == ANDERSEN || N == 1 ) {
      Ek = andersen();
    } else if ( thermostat == LANGEVIN ) {
      Ek = langevin(tstatdt);
    } else if ( thermostat == VRESCALE ) {
      Ek = vrescale(tstatdt, dof);
    } else if ( thermostat == RANDOM_COLLISION_MC ) {
      Ek = random_collision_mc();
    } else if ( thermostat == RANDOM_COLLISION_LANGEVIN ) {
      Ek = random_collision_langevin(tstatdt);
    } else {
      fprintf(stderr, "unknown thermostat %d\n", thermostat);
      exit(1);
    }

    smT += 2*Ek/dof;
    smU += Ep/N;

    if ( t > nstequiv && t % nstsamp == 0 ) {
      if ( N > 1 && thermostat != ANDERSEN ) rmcom(v);
      if ( N > 1 ) rmcom(x);
      sample();
    }

    if ( t % nstblah == 0 )
      printf("t %g: Ek %g, Ep %g, Ek + Ep %g, T %g, U %g\n",
          1.*t, Ek, Ep, Ek + Ep, smT/t, smU/t);
  }
}



#define prvec(a, nm) printvec(a, N, nm)

static void printvec(double *a, int n, const char *name)
{
  int i;

  printf("%s:\n", name);
  for ( i = 0; i < n; i++ )
    printf("%12g ", a[i]);
  printf("\n");
}



#define prmat(m, nm) printmat((double *) m, N, nm)

static void printmat(double *m, int n, const char *name)
{
  int i, j;

  printf("%s:\n", name);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      printf("%12g ", m[i*n+j]);
    printf("\n");
  }
  printf("\n");
}



const double eval_min = 1e-8;

/* add the zero modes */
static void add_zero_modes(double c[N][N], double eval[N], double evec[N][N])
{
  int k, i, j;
  const double big = N*100;
  double v[N];

  for ( k = 0; k < N; k++ ) {
    if ( eval[k] > eval_min ) continue;
    for ( i = 0; i < N; i++ )
      v[i] = evec[i][k];
    for ( i = 0; i < N; i++ )
      for ( j = 0; j < N; j++ )
        c[i][j] += v[i] * v[j] * big;
  }
}



/* principle component analysis */
static void pca(void)
{
  int i, j;
  double c[N][N], invc[N][N], eval[N], evec[N][N], avq[N], y;

  for ( i = 0; i < N; i++ )
    avq[i] = qsum[i] / nsamps;

  for ( i = 0; i < N; i++ )
    for ( j = i; j < N; j++ ) {
      y = qqsum[i][j]/nsamps - avq[i]*avq[j];
      c[j][i] = c[i][j] = y*temp*sqrt(mass[i]*mass[j]);
    }
  prmat(c, "c");

  eigsym((double *) c, eval, (double *) evec, N);
  prvec(eval, "eigenvalues (omega^(-2))");
  prmat(evec, "eigenvectors");

  /* add zero-frequency modes */
  add_zero_modes(c, eval, evec);
  luinv((double *) c, (double *) invc, N, DBL_MIN);
  prmat(invc, "c^{-1}");
}



int main(void)
{
  initff();
  domd();
  pca();
  return 0;
}

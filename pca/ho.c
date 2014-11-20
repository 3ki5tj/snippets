#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mtrand.h" /* random number generator */
#include "eig.h" /* eigenvalues */
#include "lu.h" /* LU decomposition */


#define N 13


double hess[N];
double mass[N];
double x[N];
double v[N];
double f[N];
double trigmass = 0; /* triangular mass distribution */
double temp = 1; /* temperature */
double dt = 0.005; /* time step for molecular dynamics */
double tstatdt = 0.01; /* thermostat time step */
int useandersen = 0; /* use the Andersen thermostat */

double fsmallworld = 0.0; /* frequency of adding small-world springs */
double ksmallworld = 1.0; /* stiffness of small-world springs */
int smallworld[N][N];

long long nsteps = 1000000; /* number of molecular dynamics steps */
long long nstblah = 100000; /* frequency of reporting */
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



/* compute the force, return the energy */
static double force(void)
{
  int i, j;
  double dx, k, U = 0;

  for ( i = 0; i < N; i++ ) f[i] = 0;

  if ( N == 1 ) {
    dx = -x[0];
    k = hess[0];
    f[0] = k * dx;
    U = .5 * k * dx * dx;
  } else {
    for ( i = 0; i < N - 1; i++ ) { /* loop over springs */
      dx = x[i+1] - x[i];
      k = hess[i];
      f[i]   += k * dx;
      f[i+1] -= k * dx;
      U += .5 * k * dx * dx;
    }

    // force from small-world springs
    if ( fsmallworld > 0 ) {
      for ( i = 0; i < N; i++ ) {
        for ( j = i + 2; j < N; j++ ) {
          if ( smallworld[i][j] ) {
            dx = x[j] - x[i];
            k = ksmallworld;
            f[i] += k * dx;
            f[j] -= k * dx;
            U += .5 * k * dx * dx;
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



/* velocity rescaling thermostat */
static double vrescale(double thdt, int dof)
{
  int i, j;
  double ekav = .5f*temp*dof, ek1, ek2, s, amp;

  /* hack: randomly swap two velocities to increase the randomness
   * for a complex system, we shouldn't need this */
  if ( rand01() < 0.1 ) {
    double tmp;
    i = (int) (N * rand01());
    j = ((int) ((N - 1) * rand01()) + i + 1) % N;
    tmp = mass[i] * v[i];
    v[i] = v[j] * mass[j] / mass[i];
    v[j] = tmp / mass[j];
  }

  /* normal velocity rescaling */
  for ( ek1 = 0, i = 0; i < N; i++ )
    ek1 += .5 * mass[i] * v[i] * v[i];
  amp = 2*sqrt(ek1*ekav*thdt/dof);
  ek2 = ek1 + (ekav - ek1)*thdt + (amp*gaussrand());
  if (ek2 < 1e-6) ek2 = 1e-6;
  s = sqrt(ek2/ek1);
  for (i = 0; i < N; i++)
    v[i] *= s;

  return ek2;
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

  dof = ( useandersen || N == 1 ) ? N : N - 1;
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
    if ( useandersen || N == 1 ) {
      Ek = andersen();
    } else {
      Ek = vrescale(tstatdt, dof);
    }

    smT += 2*Ek/dof;
    smU += Ep/N;

    if ( t > nstequiv && t % nstsamp == 0 ) {
      if ( N > 1 && !useandersen ) rmcom(v);
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

  if ( useandersen ) {
    luinv((double *) c, (double *) invc, N, DBL_MIN);
    prmat(invc, "c^{-1}");
  }
}



int main(void)
{
  initff();
  domd();
  pca();
  return 0;
}

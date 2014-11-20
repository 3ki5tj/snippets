/* compute the Madelung constant by Ewald sum */
#include <stdio.h>
#include <math.h>

#define N 8
double unitcell[3] = {2.0, 2.0, 2.0};
double x[N][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
double q[N] = {1,-1,-1,1,-1,1,1,-1};

static double dot(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double dist(double a[3], double b[3], double ucell[3], int ishift[3])
{
  double c[3];
  int i;
  for (i = 0; i < 3; i++ )
    c[i] = a[i] - b[i] + ucell[i]*ishift[i];
  return sqrt(dot(c, c));
}

static double ewald(int n, double (*x)[3], double *q,
    double *ucell, double alpha, int km)
{
  int i, j, l, m, ishift[3];
  double r, eself = 0, ereal = 0, erecip = 0;
  double sqrta = sqrt(alpha), k[3], k2, re, im, phase;

  /* real-space sum */
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      /* loop over possible periodic images */
      for ( ishift[0] = -1; ishift[0] <= 1; ishift[0]++ ) {
        for ( ishift[1] = -1; ishift[1] <= 1; ishift[1]++ ) {
          for ( ishift[2] = -1; ishift[2] <= 1; ishift[2]++ ) {
            r = dist(x[i], x[j], ucell, ishift);
            ereal += q[i]*q[j]*erfc(sqrta*r)/r;
          }
        }
      }
    }
  }

  /* reciprocal space sum */
  for ( i = -km; i <= km; i++ ) {
    k[0] = i*2*M_PI/ucell[0];
    for ( j = -km; j <= km; j++ ) {
      k[1] = j*2*M_PI/ucell[1];
      for ( l = -km; l <= km; l++ ) {
        if ( i == 0 && j == 0 && l == 0) continue;
        k[2] = l*2*M_PI/ucell[2];
        k2 = dot(k, k);
        re = im = 0;
        for ( m = 0; m < n; m++ ) {
          phase = dot(k, x[m]);
          re += q[m] * cos(phase);
          im += q[m] * sin(phase);
        }
        erecip += (re*re + im*im) * exp(-k2/4/alpha)/k2;
      }
    }
  }
  erecip *= 2*M_PI/(ucell[0]*ucell[1]*ucell[2]);

  /* self energy */
  for ( i = 0; i < n; i++ )
    eself += q[i]*q[i];
  eself *= -sqrt(alpha/M_PI);
  return eself + ereal + erecip;
}


int main(void)
{
  double ene = ewald(N, x, q, unitcell, 4.0, 5);
  printf("energy %g, Madelung constant %g\n", ene, 2*ene/N);
  return 0;
}

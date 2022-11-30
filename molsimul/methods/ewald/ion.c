/* compute the electric energy of a single ion
 * in a uniform negative background by Ewald sum
 * To compile and run
 *    gcc ion.c -lm && ./a.out
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double echarge = 1.6021766208e-19;
double eps0 = 8.854187817e-12;
double NA = 6.022140857e23;
double cal = 4.184;
double atm = 101325;



/* vector cross product */
static double *vcross(double *c, double *a, double *b)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

static void vscale(double *v, double s)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

static double getvol(double mat[3][3])
{
  return fabs(mat[0][0] * mat[1][1] * mat[2][2]);
}

static double getmindim(double mat[3][3])
{
  double x = fabs( mat[0][0] );
  if ( fabs(mat[1][1]) < x ) x = fabs(mat[1][1]);
  if ( fabs(mat[2][2]) < x ) x = fabs(mat[2][2]);
  return x;
}

static double inverse(double invmat[3][3], double mat[3][3])
{
  double invvol = 1. / getvol(mat);
  vscale( vcross(invmat[0], mat[1], mat[2]), invvol );
  vscale( vcross(invmat[1], mat[2], mat[0]), invvol );
  vscale( vcross(invmat[2], mat[0], mat[1]), invvol );
}

static double getdist2(double x, double y, double z, double mat[3][3])
{
  double v[3] = {x, y, z}, u[3] = {0, 0, 0}, dis2 = 0;
  int i, j;
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )
      u[i] += v[j] * mat[j][i];
    dis2 += u[i] * u[i];
  }
  return dis2;
}

static double ewald(double sigma, double tol, double mat[3][3])
{
  int i, j, l, xm[3], km[3];
  double r, k2, eself = 0, ereal = 0, erecip = 0, ebg = 0;
  double sqrta, inva;
  double vol = getvol(mat), invmat[3][3], dis2, dim, invdim;

  dim = getmindim(mat);
  inverse(invmat, mat);
  invdim = getmindim(invmat);

  sqrta = sigma > 0 ? sigma : sqrt(M_PI*invdim/dim);
  inva = M_PI * M_PI / (sqrta * sqrta);

  /* estimate the number of neighboring cells */
  for ( i = 0; i < 3; i++ ) {
    for ( xm[i] = 1; xm[i] <= 1000; xm[i]++ ) {
      r = fabs(mat[i][i]) * xm[i];
      if ( erfc(sqrta*r)/r < tol ) break;
    }
    printf("dim %d, xm %d, sqrta %g, r %g, %g\n", i, xm[i], sqrta, r, erfc(sqrta*r)/r);
  }

  /* estimate the number of wave vectors */
  for ( i = 0; i < 3; i++ ) {
    for ( km[i] = 1; km[i] <= 1000; km[i]++ ) {
      r = fabs(invmat[i][i]) * km[i];
      k2 = r * r;
      if ( exp(-k2*inva)/k2/M_PI/vol < tol ) break;
    }
    printf("dim %d, km %d, inva %g, k %g, %g\n", i, km[i], inva, r, exp(-k2*inva)/k2/M_PI/vol);
  }

  /* real-space sum */
  for ( i = -xm[0]; i <= xm[0]; i++ ) {
    for ( j = -xm[1]; j <= xm[1]; j++ ) {
      for ( l = -xm[2]; l <= xm[2]; l++ ) {
        dis2 = getdist2(i, j, l, mat);
        if ( dis2 <= 0 ) continue;
        r = sqrt(dis2);
        ereal += erfc(sqrta*r)/r;
      }
    }
  }

  /* reciprocal-space sum */
  for ( i = -km[0]; i <= km[0]; i++ ) {
    for ( j = -km[1]; j <= km[1]; j++ ) {
      for ( l = -km[2]; l <= km[2]; l++ ) {
        k2 = getdist2(i, j, l, invmat);
        if ( k2 <= 0 ) continue;
        erecip += exp(-k2*inva)/k2;
      }
    }
  }
  erecip /= M_PI * vol;

  /* self energy */
  eself = -2 * sqrta / sqrt(M_PI);

  /* background energy */
  ebg = -M_PI / (sqrta * sqrta * vol);
  printf("%g %g %g %g\n", ereal, erecip, eself, ebg);
  return eself + ereal + erecip + ebg;
}



int main(int argc, char **argv)
{
  double ene, epot, sigma = 0.0, tol = 1e-14;
  double mat[3][3] = {{10, 0, 0}, {0, 10, 0}, {0, 0, 10}};
  double K0, P0, vol = getvol(mat), pres;

  if ( argc > 1 ) sigma = atof(argv[1]);
  K0 = echarge * echarge / (4 * M_PI * eps0) * NA * 1e10 / (cal*1e3);
  // 332.0716 -- CHARMM
  ene = ewald(sigma, tol, mat);
  epot = 0.5 * K0 * ene;
  P0 = K0 / (NA*1e7/cal) * 1e40 / atm;
  pres = 0.5 * P0 * ene / (3 * vol);
  printf("ene = %.15f, u = %.14f, p = %.14f, K0 %.8f\n", ene, epot, pres, K0);
  return 0;
}

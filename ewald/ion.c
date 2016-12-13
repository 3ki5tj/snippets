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
  return mat[0][0] * mat[1][1] * mat[2][2];
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
      u[i] += mat[i][j]*v[j];
    dis2 += u[i] * u[i];
  }
  return dis2;
}

static double ewald(double sqrta, int xm, int km,
    double mat[3][3])
{
  int i, j, l;
  double r, k2, eself = 0, ereal = 0, erecip = 0, ebg = 0;
  double inva = M_PI * M_PI / (sqrta * sqrta);
  double vol = getvol(mat), invmat[3][3], dis2;

  /* real-space sum */
  for ( i = -xm; i <= xm; i++ ) {
    for ( j = -xm; j <= xm; j++ ) {
      for ( l = -xm; l <= xm; l++ ) {
        dis2 = getdist2(i, j, l, mat);
        if ( r <= 0 ) continue;
        r = sqrt(dis2);
        ereal += erfc(sqrta*r)/r;
      }
    }
  }

  inverse(invmat, mat);

  /* reciprocal-space sum */
  for ( i = -km; i <= km; i++ ) {
    for ( j = -km; j <= km; j++ ) {
      for ( l = -km; l <= km; l++ ) {
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
  return eself + ereal + erecip + ebg;
}


int main(int argc, char **argv)
{
  double ene, epot, sigma = 0.0;
  double mat[3][3] = {{10, 0, 0}, {0, 10, 0}, {0, 0, 10}};
  int nterms = 100, kterms = 100;
  double K0 = 332.0716, P0, vol = getvol(mat), pres;

  if ( argc > 1 ) sigma = atof(argv[1]);
  if ( sigma == 0 ) sigma = 5.6 / mat[0][0];
  if ( argc > 2 ) nterms = atoi(argv[2]);
  if ( argc > 3 ) kterms = atoi(argv[3]);
  //K0 = echarge * echarge / (4 * M_PI * eps0) * NA * 1e10 / (cal*1e3);
  ene = ewald(sigma, nterms, kterms, mat);
  epot = 0.5 * K0 * ene;
  P0 = K0 / (NA*1e7/cal) * 1e40 / atm;
  pres = 0.5 * P0 * ene / (3 * vol);
  printf("%.15f, u = %.14f, p = %.14f, K0 %.10f\n", ene, epot, pres, K0);
  return 0;
}

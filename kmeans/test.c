#include "kmeans.h"
#include "mtrand.h"



int dim = 2;
int n = 4000;
int K = 4;
double **dat;



static void mkdata(void)
{
  int i, k, n1;
  double x[2], c, s;
  double av[4][2] = {{0, 0}, {1.5, 0}, {-2, 2}, {-2, -2}};
  double sig[4][2] = {{1, 1}, {0.2, 1.5}, {1., 0.1}, {1., 0.4}};
  double theta[4] = {0, 0, M_PI/4, -M_PI/4};
  FILE *fp;

  xnew(dat, n);
  n1 = n / K;
  for ( i = 0; i < n; i++ )
    xnew(dat[i], dim);
  for ( k = 0; k < K; k++ ) {
    for ( i = k*n1; i < (k+1)*n1; i++ ) {
      x[0] = sig[k][0] * gaussrand();
      x[1] = sig[k][1] * gaussrand();
      c = cos(theta[k]);
      s = sin(theta[k]);
      dat[i][0] = av[k][0] + x[0] * c + x[1] * s;
      dat[i][1] = av[k][1] + x[0] * (-s) + x[1] * c;
    }
  }

  fp = fopen("x.dat", "w");
  for ( i = 0; i < n; i++ )
    fprintf(fp, "%g %g\n", dat[i][0], dat[i][1]);
  fclose(fp);
}



int main(void)
{
  kmeans_t *km;
  int i;

  mkdata();
  km = kmeans_open(dim, n, K, dat);
  for ( i = 0; i < 1000; i++ ) {
    kmeans_estep(km);
    kmeans_mstep(km);
  }
  kmeans_print(km);
  kmeans_close(km);
  return 0;
}

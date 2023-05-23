#include "kmeans.h"
#include "mtrand.h"



int dim = 2;
int n = 4000;
int K = 4;
double **dat;


typedef struct {
  double xc, yc;
  double a, b, theta;
  int npt;
} ellipse_t;



static void mkdata(void)
{
  int i, k, n1;
  double x[2], c, s;
  ellipse_t elp[4] = {
    { 0.0,  0.0, 1.0, 1.0, 0, 1000},
    { 1.5,  0.0, 0.2, 1.5, 0, 1000},
    {-2.0,  2.0, 1.0, 0.1, -M_PI/4, 1000},
    {-2.0, -2.0, 1.0, 0.4, M_PI/4, 1000},
  };
  FILE *fp;

  xnew(dat, n);
  n1 = n / K;
  for ( i = 0; i < n; i++ )
    xnew(dat[i], dim);
  for ( k = 0; k < K; k++ ) {
    for ( i = k*n1; i < (k+1)*n1; i++ ) {
      x[0] = elp[k].a * gaussrand();
      x[1] = elp[k].b * gaussrand();
      c = cos(elp[k].theta);
      s = sin(elp[k].theta);
      dat[i][0] = elp[k].xc + x[0] * c - x[1] * s;
      dat[i][1] = elp[k].yc + x[0] * s + x[1] * c;
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

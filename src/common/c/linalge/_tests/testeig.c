#include "eig.h"


#define N 4

int main(void)
{
  double mat[N][N] = {
    {3, 1, 1, 0},
    {1, 4, 1, 1},
    {1, 1, 5, 1},
    {0, 1, 1, 6}
  }, eval[N], evec[N][N];
  int i, j;

  printf("original matrix:\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%11g ", mat[i][j]);
    printf("\n");
  }
  printf("\n");

  eigsym((double *) mat, eval, (double *) evec, N);

  printf("eigenvalues:\n");
  for ( i = 0; i < N; i++ ) {
    printf("%11g ", eval[i]);
  }
  printf("\n");

  printf("eigenvectors:\n");
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ )
      printf("%11g ", evec[i][j]);
    printf("\n");
  }

  return 0;
}



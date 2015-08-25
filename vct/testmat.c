#include "mat.h"



int main(void)
{
  double m1[D][D] = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  double invm1[D][D] = {{1, -3, 2}, {-3, 3, -1}, {2, -1, 0}};
  double m2[D][D] = {{1, 2, 3}, {2, 4, 6}, {3.00000000000001, 6, 9}};
  double m3[D][D] = {{1, 2, 3}, {4, 5, 5}, {5, 7, 8}};
  double m4[D][D] = {{1, -1, 0}, {-1, 2, -1}, {0, -1, 1}};
  double t1[D][D], t2[D][D], t3[D][D], v[D];
  double tol = 1e-7;
  int i, j, n;

  /* test the determinant and inverse */
  printf("det %g (-1), %g (0)\n", mdet(m1), mdet(m2));
  minv(t1, m1);
  printf("Inverse matrix:\n");
  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      printf("%8.3f(%8.3f)\t", t1[i][j], invm1[i][j]);
    }
    printf("\n");
  }

  printf("solvezero I:\n");
  n = msolvezero(m2, t2, tol);
  for ( i = 0; i < n; i++ ) {
    printf("solution %d/%d:\t", i, n);
    for ( j = 0; j < D; j++ ) {
      printf("%8.3f\t", t2[i][j]);
    }
    printf("\n");
  }

  printf("solvezero II:\n");
  n = msolvezero(m3, t3, tol);
  for ( i = 0; i < n; i++ ) {
    printf("solution %d/%d:\t", i, n);
    for ( j = 0; j < D; j++ ) {
      printf("%8.3f\t", t3[i][j]);
    }
    printf("\n");
  }

  meigsys(v, t1, m4, 1);
  for ( i = 0; i < D; i++ ) {
    printf("eigval %8.3f: eigvec ", v[i]);
    for ( j = 0; j < D; j++ ) {
      printf("%8.3f ", t1[i][j]);
    }
    printf("\n");
  }
  return 0;
}

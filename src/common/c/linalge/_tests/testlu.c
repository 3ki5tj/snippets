#include <stdio.h>
#include "lu.h"

int main(void)
{
  double a[3][3] = {{1,2,3}, {2,3,5}, {3,2,1}}, b[3][3];
  int i, j, n = 3;

  luinv((double *)a, (double *)b, n, 1e-10);

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      printf("%8.3f", b[i][j]);
    }
    printf("\n");
  }

  return 0;
}

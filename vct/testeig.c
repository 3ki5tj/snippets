#include "mat.h"

int main(void)
{
  double a[3][3] = {{6.67838e6, 4.19371e6, -157685}, {4.19371e6, 2.70432e6, -109931}, {-157685, -109931, 17187.6}};
  double s[3], v[3][3];
  int ret;
  ret = meigsys(s, v, a, 1);
  printf("ret %d\n", ret);
  printf("eigval: %g, %g, %g\n", s[0], s[1], s[2]);
  return 0;
}

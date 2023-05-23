#define D 3
#include "mat.h"

#if D == 2

double a1[D][D] = {{4, 2}, {2, 1}};
double a2[D][D] = {{2, 1}, {-1, 0}};

static void testeig(double a[2][2])
{
  double s[2], v[2][2];
  int ret;

  ret = meigsys(s, v, a, 1);
  printf("ret %d: eigval: %14.8f, %14.8f, tol %g\n", ret, s[0], s[1], msolvezero_lasttol);
  printf("vec 1: %14.8f %14.8f\n", v[0][0], v[0][1]);
  printf("vec 2: %14.8f %14.8f\n\n", v[1][0], v[1][1]);
}

#else // D == 3

/* for the following matrix, the tolerance must be
 * increased significantly for a solution */
double a1[3][3] = {{ 247146.37794751967886,  46767.34439483784081,-234432.36461535561830},{  46767.34439483784081,   9741.70930413716087, -44517.71041657723981},{-234432.36461535561830, -44517.71041657723981, 223992.19184470266919}};

/* the two largest eigenvectors of the following matrix are degenerate */
double a2[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 0.5}};

static void testeig(double a[3][3])
{
  double s[3], v[3][3];
  int ret;

  ret = meigsys(s, v, a, 1);
  printf("ret %d: eigval: %14.8f, %14.8f, %14.8f, tol %g\n", ret, s[0], s[1], s[2], msolvezero_lasttol);
  printf("vec1: %14.8f %14.8f %14.8f\n", v[0][0], v[0][1], v[0][2]);
  printf("vec2: %14.8f %14.8f %14.8f\n", v[1][0], v[1][1], v[1][2]);
  printf("vec3: %14.8f %14.8f %14.8f\n\n", v[2][0], v[2][1], v[2][2]);
}

#endif


int main(void)
{
  testeig(a1);
  testeig(a2);

  return 0;
}

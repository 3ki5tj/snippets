#include "fsamp.h"
#include "mtrand.h"

double p1[] = {0.50, 0.50};
double p2[] = {0.33, 0.67};

static void test1(int n, double *pa, double *pb,
    int ntrial)
{
  int t, i, id = 0;
  double *cnt, *his, *p;

  xnew(cnt, n);
  xnew(his, n);

  id = 0;
  p = pa;
  for ( t = 0; t < ntrial; t++ ) {
    i = fsamp_select(n, cnt, p);
    cnt[i] += 1;
    his[i] += 1;

    if ( t % 200 == 100 ) {
      fsamp_truncate(n, cnt, p);
      id = !id;
      p = ( id == 0 ) ? pa : pb;
    }
  }

  for ( i = 0; i < n; i++ ) {
    printf("%4d: %.0f\n", i, his[i]);
  }
  free(cnt);
  free(his);
}

int main(void)
{
  test1(2, p1, p2, 1000000);
  return 0;
}

/* basic tests for tconv.h */
#include "tconv.h"


double p1[] = {0.15, 0.24, 0.26, 0.35};
double p2[] = {0.15, 0.15, 0.25, 0.45};
double p3[] = {0.05, 0.1, 0.15, 0.7};

static void test_printmat(void)
{
  tc_printmat(1, 4, p1, "Type 1, smooth 1");
  tc_printmat(1, 4, p2, "Type 1, smooth 2");
  tc_printmat(1, 4, p3, "Type 1, peak");

  tc_printmat(2, 4, p1, "Type 2, smooth 1");
  tc_printmat(2, 4, p2, "Type 2, smooth 2");
  tc_printmat(2, 4, p3, "Type 2, peak");
}


/* conduct sampling */
__inline static void tc_select_wrap(int type, int n, const double *p,
    const char *name, int nsamp)
{
  int r, c, t;
  double *cp, *mat;

  xnew(cp, n + 1);
  xnew(mat, n * n);
  for ( c = 0; c < n; c++ ) {
    for ( t = 0; t < nsamp; t++ ) {
      if ( type == 1 ) {
        r = tc1_select(n, c, p, cp);
      } else {
        r = tc2_select(n, c, p, cp);
      }
      mat[r*n + c] += 1;
    }
    for ( r = 0; r < n; r++ )
      mat[r*n + c] /= nsamp;
  }

  tc_pmat(n, mat, p, name, 0);
  tc_printmat(type, n, p, name);

  free(cp);
  free(mat);
}



static void test_select(void)
{
  int nsamp = 1000000;

  tc_select_wrap(1, 4, p1, "Type 1, smooth 1", nsamp);
  tc_select_wrap(1, 4, p2, "Type 1, smooth 2", nsamp);
  tc_select_wrap(1, 4, p3, "Type 1, dominant", nsamp);

  tc_select_wrap(2, 4, p1, "Type 2, smooth 1", nsamp);
  tc_select_wrap(2, 4, p2, "Type 2, smooth 2", nsamp);
  tc_select_wrap(2, 4, p3, "Type 2, dominant", nsamp);
}


double p1_us[] = {0.26, 0.24, 0.15, 0.35};
double p2_us[] = {0.25, 0.15, 0.15, 0.45};
double p3_us[] = {0.7, 0.1, 0.15, 0.05};

/* conduct sampling (unsorted distribution) */
__inline static void tc_select_us_wrap(int type, int n, const double *p,
    const char *name, int nsamp)
{
  int r, c, t, *id;
  double *ps, *cp, *mat;

  xnew(id, n);
  xnew(ps, n);
  xnew(cp, n + 1);
  xnew(mat, n * n);
  for ( c = 0; c < n; c++ ) {
    for ( t = 0; t < nsamp; t++ ) {
      r = tc_select(type, n, c, p, id, ps, cp);
      mat[r*n + c] += 1;
    }
    for ( r = 0; r < n; r++ )
      mat[r*n + c] /= nsamp;
  }

  tc_pmat(n, mat, p, name, 0);

  free(id);
  free(ps);
  free(cp);
  free(mat);
}



static void test_select_us(void)
{
  int nsamp = 1000000;

  tc_select_us_wrap(1, 4, p1_us, "Type 1, smooth 1", nsamp);
  tc_select_us_wrap(1, 4, p2_us, "Type 1, smooth 2", nsamp);
  tc_select_us_wrap(1, 4, p3_us, "Type 1, dominant", nsamp);

  tc_select_us_wrap(2, 4, p1_us, "Type 2, smooth 1", nsamp);
  tc_select_us_wrap(2, 4, p2_us, "Type 2, smooth 2", nsamp);
  tc_select_us_wrap(2, 4, p3_us, "Type 2, dominant", nsamp);
}


int main(int argc, char **argv)
{
  if ( argc > 1 ) {
    if ( argv[1][0] == 'p' ) {
      test_printmat();
    } else if ( argv[1][0] == 's' ) {
      test_select();
    } else {
      test_select_us();
    }
  }
  return 0;
}

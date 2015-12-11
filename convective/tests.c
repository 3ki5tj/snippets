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


static void test_select(void)
{
  int nsamp = 1000000;

  tc_sample(1, p1, "Type 1, smooth 1", nsamp);
  tc_sample(1, p2, "Type 1, smooth 2", nsamp);
  tc_sample(1, p3, "Type 1, dominant", nsamp);
  
  tc_sample(2, p1, "Type 2, smooth 1", nsamp);
  tc_sample(2, p2, "Type 2, smooth 2", nsamp);
  tc_sample(2, p3, "Type 2, dominant", nsamp);
}


int main(void)
{
  test_printmat();
  test_select();
  return 0;
}

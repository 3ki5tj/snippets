#ifndef FSAMP_H__
#define FSAMP_H__



/* determininistic sampling based on frequency */


#include "util.h"



/* select an item with the minimal height */
__inline static
int fsamp_select(int n, double *cnt, const double *p)
{
  int i, im = 0;
  double y, ym = DBL_MAX;

  for ( i = 0; i < n; i++ ) {
    if ( p[i] <= 0 ) continue;

    y = cnt[i] / p[i];
    if ( y < ym ) {
      im = i;
      ym = y;
    }
  }

  return im;
}



/* truncate the counts based on the given probability p */
__inline static
void fsamp_truncate(int n, double *cnt, const double *p)
{
  int i;
  double ptot = 0, ctot = 0, have;

  for ( i = 0; i < n; i++ ) {
    ptot += p[i];
    ctot += cnt[i];
  }

  have = ctot / ptot;

  for ( i = 0; i < n; i++ ) {
    cnt[i] -= have * p[i];
  }
}



#endif /* FSAMP_H__ */

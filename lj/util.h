#ifndef UTIL_H__
#define UTIL_H__



#include "mdutil.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



/* return the larger of a and b */
__inline double dblmax(double a, double b)
{
  return a > b ? a : b;
}



/* return the smaller of a and b */
__inline double dblmin(double a, double b)
{
  return a < b ? a : b;
}



#endif


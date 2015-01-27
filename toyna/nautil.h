#ifndef UTIL_H__
#define UTIL_H__



#include "vct.h"



#ifndef PI
#define PI 3.141592653589793
#endif



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



#define D2R(x)  ((x)*PI/180.0)
#define R2D(x)  ((x)*180.0/PI)



#endif


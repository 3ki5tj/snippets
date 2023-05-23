#ifndef CAGOUTIL_H__
#define CAGOUTIL_H__



#include "mdutil.h"
#include "mtrand.h"



#ifndef xnew
#define xnew(x, n) { \
  if ( (x = calloc((n), sizeof(*(x)))) == NULL ) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



#ifndef xrenew
#define xrenew(x, n) { \
  if ( (x = realloc(x, sizeof(*(x)) * (n))) == NULL ) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



/* return the larger of x and y */
__inline static int intmax(int x, int y)
{
  return x > y ? x : y;
}



/* return the smaller of x and y */
__inline static int intmin(int x, int y)
{
  return x < y ? x : y;
}



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



/* update the MC move size according to the acceptance ratio */
__inline static void update_mcamp(double *amp, double acc,
    double *nacc, double *ntot)
{
  double x = sqrt( *nacc / *ntot / acc );
  if ( x > 2 ) x = 2;
  else if ( x < 0.5 ) x = 0.5;
  *amp *= x;
  fprintf(stderr, "acc %g%%, amp %g\n", 100*(*nacc)/(*ntot), *amp);
  *nacc = 0;
  *ntot = DBL_MIN;
}



const char *aanames[] = {
  "GLY", "ALA", "VAL", "LEU", "ILE",
  "PRO", "SER", "THR", "CYS", "MET",
  "ASN", "GLN", "ASP", "GLU",
  "LYS", "ARG", "HIS",
  "PHE", "TYR", "TRP"};

const char aaletters[] = "GAVLIPSTCMNQDEKRHFYW";

/* residue name to integer of amino acid */
__inline static int res2iaa(const char *res)
{
  int i;

  for ( i = 0; i < 20; i++ ) {
    if ( strcmp(res, aanames[i]) == 0 ) {
      return i;
    }
  }
  return -1;
}



#endif /* CAGOUTIL_H__ */


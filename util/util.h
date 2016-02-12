#ifndef UTIL_H__
#define UTIL_H__





#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



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




/* String routines */



/* remove leading and trailing spaces */
__inline static char *strstrip(char *s)
{
  char *p, *q;

  /* remove trailing spaces */
  for ( p = s + strlen(s) - 1; p >= s && isspace(*p); p-- )
    *p = '\0';

  /* remove leading spaces */
  for ( p = s; *p && isspace(*p); p++ )  ;
  if ( p != s )
    for ( q = s; (*q++ = *p++) != '\0'; ) ;
  return s;
}



#define strcmpfuzzy(s, t) strncmpfuzzy(s, t, INT_MAX)

/* comparison, ignoring cases, spaces and punctuations */
__inline static int strncmpfuzzy(const char *s, const char *t, int n)
{
  int is, it, i;
  const char cset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789()[]{}";

  for ( i = 0; i < n; s++, t++, i++ ) {
    while ( *s != '\0' && strchr(cset, *s) == NULL ) s++;
    while ( *t != '\0' && strchr(cset, *t) == NULL ) t++;
    is = tolower(*s);
    it = tolower(*t);
    if ( is != it ) return is - it;
    if ( *s == '\0' ) return 0;
  }
  return 0;
}



/* check if `s' starts with `t', using fuzzy comparison */
__inline static int strstartswith(const char *s, const char *t)
{
  return strncmpfuzzy(s, t, strlen(t)) == 0;
}



/* check if a string is a nonnegative integer */
__inline static int striscnum(const char *s)
{
  for ( ; *s; s++ ) if ( !isdigit(*s) ) return 0;
  return 1;
}



#endif /* UTIL_H__ */

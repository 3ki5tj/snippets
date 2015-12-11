#ifndef TCONV_H__
#define TCONV_H__



#include "mtrand.h"


/* determine the case ID for
 * the type I convective transition matrix */
__inline static
int tc1_getcase(int n, const double *p)
{
  int i;
  double cp = p[n-2];

  for ( i = 1; i <= n - 2; i++ ) {
    cp += p[n-2-i];
    if ( p[n-1] <= cp ) return i;
  }
  return n - 1;
}


/* determine the case ID for
 * the type II convective transition matrix */
__inline static
int tc2_getcase(int n, const double *p)
{
  int i;
  double cp = p[0];

  for ( i = 1; i <= n - 2; i++ ) {
    cp += p[i];
    if ( p[n-1] <= cp ) return i;
  }
  return n - 1;
}


/* return the transition probability
 * for convective transition matrix
 * for the dominant case */
__inline static
double tc_getprob0(int n, int r, int c, const double *p)
{
  int i;
  double cp;

  if ( c < n - 1 ) {
    return r == n - 1;
  } else if ( r < n - 1 ) {
    return p[r]/p[n-1];
  } else {
    for ( cp = 0, i = 0; i < n-1; i++ )
      cp += p[i];
    return 1 - cp/p[n-1];
  }
}


/* return the transition probability
 * for type I convective transition matrix */
__inline static
double tc1_getprob(int n, int r, int c, const double *p)
{
  int i, j;
  double cp;

  i = tc1_getcase(n, p);

  if ( i < n - 1 ) {
    if ( c == 0 ) {

      return r == n - 1;

    } else if ( c < n - i - 1 ) {

      cp = p[c-1] / p[c];
      return ( r == c - 1 ) ? cp : ( r == n - 1 ) ? 1 - cp : 0;

    } else if ( c < n - 1 ) {

      for ( cp = 0, j = n - i - 1; j < n - 1; j++ )
        cp += p[j];
      cp = ( p[n - 1] - p[n - i - 2] ) / cp;

      return ( r == n - 1 ) ? cp : ( r == n - i - 2 ) ? 1 - cp : 0;

    } else {

      if ( r == n - i - 2 ) {
        for ( cp = 0, j = n - i - 1; j < n - 1; j++ )
          cp += p[j];
        return 1 - cp / p[n-1];
      } else if ( n - i - 1 <= r && r < n - 1 ) {
        return p[r] / p[n-1];
      } else {
        return 0;
      }

    }
  } else {
    return tc_getprob0(n, r, c, p);
  }
}



/* return the transition probability
 * for type II convective transition matrix */
__inline static
double tc2_getprob(int n, int r, int c, const double *p)
{
  int i, j;
  double cp;

  i = tc2_getcase(n, p);

  if ( i < n - 1 ) {
    if ( c < i ) {

      if ( r <  i ) {
        return 0;
      } else {
        for ( cp = 0, j = 0; j < i; j++ )
          cp += p[j];
        if ( r == i ) {
          return 1 - (p[n-1] - p[i])/cp;
        } else {
          return (p[r] - p[r-1])/cp;
        }
      }

    } else if ( c < n - 1 ) {

      return r == c + 1;

    } else { /* c == n - 1 */

      if ( r < i ) {
        return p[r]/p[n-1];
      } else if ( r == i ) {
        for ( cp = 0, j = 0; j < i; j++ )
          cp += p[j];
        return 1 - cp/p[n-1];
      } else {
        return 0;
      }

    }
  } else {
    return tc_getprob0(n, r, c, p);
  }
}



/* select an item of an array of n items according to
 * probabilites p[0..n-1],
 * cp[0..n] are space for the cumulative probabilities
 * return n on failure */
__inline static
int heatbath_select(int n, const double *p, double *cp)
{
  int i;
  double r;

  /* compute the cumulative probabilities */
  cp[0] = 0;
  for ( i = 0; i < n; i++ )
    cp[i+1] = cp[i] + p[i];

  r = cp[n] * rand01();

  /* find the bracket containing r
   * the following algorithm is good for a distribution
   * dominated by p[n-1] */
  for ( i = n - 1; i >= 0; i-- ) {
    if ( r >= cp[i] ) {
      break;
    }
  }
  return i;
}



/* select a state according to
 * the type I convective transition matrix */
__inline static
int tc1_select(int n, int c, const double *p, double *cp)
{
  int i, j;
  double rn;

  /* 1. compute the case number 1..n-1
   *    also compute the cumulative probabilities */
  i = n - 1;
  cp[i] = 0;
  for ( j = n - 2; j >= 0; j-- ) {
    cp[j] = cp[j + 1] + p[j];
    if ( p[n - 1] <= cp[j] ) {
      i = n - j - 2;
      break;
    }
  }

  if ( c < n - 1 ) {
    if ( i == n - 1 || c == 0 ) {
      return n - 1;
    } else if ( c < n - i - 1 ) {
      return ( rand01() < p[c - 1] / p[c] )  ? c - 1 :  n - 1;
    } else {
      return ( rand01() < ( p[n - 1] - p[n - i - 2] ) / cp[n - i - 1] ) ? n - 1 : n - i - 2;
    }
  }

  rn = p[n - 1] * rand01();

  /* find the bracket containing rn */
  for ( j = n - 2; j >= n - i - 1; j-- )
    if ( rn <= cp[j] )
      break;

  return j < 0 ? n - 1 : j;
}



/* return the transition probability
 * for type II convective transition matrix */
__inline static
double tc2_select(int n, int c, const double *p, double *cp)
{
  int i, j;
  double rn;

  /* 1. compute the type */
  cp[0] = 0;
  for ( i = 0; i < n - 1; i++ ) {
    cp[i + 1] = cp[i] + p[i];
    if ( p[n - 1] <= cp[i + 1] ) break;
  }
  //printf("case %d, %g, %g vs %g\n", i, p[i], cp[i], p[n-1]); getchar();

  if ( c < n - 1 ) { 
    if ( i < n - 1 ) {
      if ( c < i ) {

        /* select j with probability (p[j]-p[j-1])/(p[0]+...+p[i-1]);
         * Note: (p[i+1]-p[i]) + ... + (p[n-1] - p[n-2]) = p[n-1] - p[i] */
        rn = rand01();
        for ( j = n - 1; j > i; j-- ) {
          if ( rn >= (p[j - 1] - p[n - 1])/cp[i] + 1 )
            break;
        }
        return j;

      } else {

        return c + 1;

      }
    }

    return n - 1;
  }

  /* select j with probability p[j]/p[n-1]; */
  rn = rand01() * p[n-1];
  for ( j = i; j >= 0; j-- )
    if ( rn >= cp[j] )
      break;
  return j < 0 ? n - 1 : j;
}



/* print a transition matrix */
__inline static
void tc_pmat(int n, const double *tmat, const double *p,
    const char *name, int casenum)
{
  int r, c;
  double sum;

  if ( name != NULL ) {
    printf("%s (case %d):\n", name, casenum);
  }

  printf("p:  ");
  for ( r = 0; r < n; r++ ) {
    printf("%7.5f ", p[r]);
  }
  printf("\n---------------------------------------------\n");
  for ( r = 0; r < n; r++ ) {
    printf("%d:  ", r);
    sum = 0;
    for ( c = 0; c < n; c++ ) {
      sum += p[c] * tmat[r*n + c];
      printf("%7.5f ", tmat[r*n + c]);
    }
    printf("| %7.5f\n", sum);
  }
  printf("\n");
}


#define NMAX 10
#define N2MAX (NMAX * NMAX)

__inline static
void tc_printmat(int type, int n, const double *p, const char *name)
{
  int r, c, casenum = 0;
  double tmat[N2MAX];

  if ( type == 1 ) {
    casenum = tc1_getcase(n, p);
  } else {
    casenum = tc2_getcase(n, p);
  }

  for ( r = 0; r < n; r++ ) {
    for ( c = 0; c < n; c++ ) {
      if ( type == 1 ) {
        tmat[r*n+c] = tc1_getprob(n, r, c, p);
      } else {
        tmat[r*n+c] = tc2_getprob(n, r, c, p);
      }
    }
  }

  tc_pmat(n, tmat, p, name, casenum);
}



/* conduct sampling */
static void tc_sample(int type, const double *p,
    const char *name, int nsamp)
{
  int r, c, t;
  int n = 4;
  double cp[NMAX+1], mat[NMAX*NMAX] = {0};

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
}


#endif /* TCONV_H__ */

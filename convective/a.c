#include <stdio.h>


/* determine the case ID for
 * the type I convective transition matrix */
static int tc1_getcase(int n, const double *p)
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
static int tc2_getcase(int n, const double *p)
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
__inline double tc_getprob0(int n, int r, int c, const double *p)
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
static double tc1_getprob(int n, int r, int c, const double *p)
{
  int i, j;
  double cp;

  i = tc1_getcase(n, p);

  if ( i < n - 1 ) {
    if ( c == 0 ) {

      return r == n - 1;

    } else if ( c + 1 < n - i ) {

      return (r == c - 1) ? p[c-1]/p[c] : (r == n - 1) ? 1 - p[c-1]/p[c] : 0;

    } else if ( c < n - 1 ) {

      if ( r + 1 == n - i - 1 ) {
        for ( cp = 0, j = n - i - 1; j < n - 1; j++ )
          cp += p[j];
        return 1 - (p[n-1] - p[n-i-2]) / cp;
      } else if ( r == n - 1 ) {
        for ( cp = 0, j = n - i - 1; j < n - 1; j++ )
          cp += p[j];
        return (p[n-1] - p[n-i-2]) / cp;
      } else {
        return 0;
      }

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
static double tc2_getprob(int n, int r, int c, const double *p)
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



/* print a transition matrix */
static void tmat_print(int n, const double *tmat, const double *p,
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


static void tc1_print(int n, const double *p, const char *name)
{
  int r, c, casenum = 0;
  double tmat[100];

  casenum = tc1_getcase(n, p);

  for ( r = 0; r < n; r++ ) {
    for ( c = 0; c < n; c++ ) {
      tmat[r*n+c] = tc1_getprob(n, r, c, p);
    }
  }

  tmat_print(n, tmat, p, name, casenum);
}



static void tc2_print(int n, const double *p, const char *name)
{
  int r, c, casenum = 0;
  double tmat[100];

  casenum = tc2_getcase(n, p);

  for ( r = 0; r < n; r++ ) {
    for ( c = 0; c < n; c++ ) {
      tmat[r*n+c] = tc2_getprob(n, r, c, p);
    }
  }

  tmat_print(n, tmat, p, name, casenum);
}



int main(void)
{
  double p1[] = {0.15, 0.25, 0.25, 0.35};
  double p2[] = {0.15, 0.15, 0.25, 0.45};
  double p3[] = {0.05, 0.1, 0.15, 0.7};

  tc1_print(4, p1, "Type 1, smooth 1");
  tc1_print(4, p2, "Type 1, smooth 2");
  tc1_print(4, p3, "Type 1, peak");

  tc2_print(4, p1, "Type 2, smooth 1");
  tc2_print(4, p2, "Type 2, smooth 2");
  tc2_print(4, p3, "Type 2, peak");
  return 0;
}

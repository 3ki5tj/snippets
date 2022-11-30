#ifndef TCONV_H__
#define TCONV_H__



#include "mtrand.h"
#include "util.h"
#include "fsamp.h"



#define TC_HEATBATH     0x00000000
#define TC_MATRIX1      0x00000001
#define TC_MATRIX2      0x00000002
#define TC_METROPOLIS   0x00000003
#define TC_MATRIXMASK   0x000000ff

/* sampling using frequency instead of probability */
#define TC_FREQ         0x00010000



/* determine the case ID for
 * the type I convective transition matrix */
__inline static
int tc1_getcase(int n, const double *p)
{
  int i;
  double cp = p[n - 2];

  for ( i = 1; i <= n - 2; i++ ) {
    cp += p[n - 2 - i];
    if ( p[n - 1] <= cp ) return i;
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
    if ( p[n - 1] <= cp ) return i;
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
    return p[r] / p[n - 1];
  } else {
    for ( cp = 0, i = 0; i < n-1; i++ )
      cp += p[i];
    return 1 - cp/p[n-1];
  }
}



/* return the transition probability of all destinations
 * for convective transition matrix
 * for the dominant case */
__inline static
void tc_getprobarr0(int n, int c, const double *p, double *tp)
{
  int i;
  double cp;

  if ( c < n - 1 ) {
    for ( i = 0; i < n - 1; i++ ) {
      tp[i] = 0;
    }
    tp[n - 1] = 1;
  } else {
    cp = 0;
    for ( i = 0; i < n - 1; i++ ) {
      tp[i] = p[i] / p[n-1];
      cp += p[i];
    }
    tp[n - 1] = cp / p[n-1];
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



/* return the transition probability of all destinations
 * for type I convective transition matrix */
__inline static
void tc1_getprobarr(int n, int c, const double *p, double *tp)
{
  int i, j;
  double cp;

  i = tc1_getcase(n, p);

  if ( i < n - 1 ) {
    for ( j = 0; j < n; j++ ) {
      tp[j] = 0;
    }

    if ( c == 0 ) {

      tp[n - 1] = 1;

    } else if ( c < n - i - 1 ) {

      cp = p[c - 1] / p[c];
      tp[c - 1] = cp;
      tp[n - 1] = 1 - cp;

    } else if ( c < n - 1 ) {

      cp = 0;
      for ( j = n - i - 1; j < n - 1; j++ )
        cp += p[j];
      cp = ( p[n - 1] - p[n - i - 2] ) / cp;

      tp[n - 1] = cp;
      tp[n - i - 2] = 1 - cp;

    } else {

      cp = 0;
      for ( j = n - i - 1; j < n - 1; j++ ) {
        tp[j] = p[j] / p[n - 1];
        cp += tp[j];
      }
      tp[n - i - 2] = 1 - cp;

    }
  } else {
    tc_getprobarr0(n, c, p, tp);
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



/* return the transition probability of all destinations
 * for type II convective transition matrix */
__inline static
void tc2_getprobarr(int n, int c, const double *p, double *tp)
{
  int i, j;
  double cp;

  i = tc2_getcase(n, p);

  if ( i < n - 1 ) {

    for ( j = 0; j < n; j++ ) {
      tp[j] = 0;
    }

    if ( c < i ) {

      for ( cp = 0, j = 0; j < i; j++ )
        cp += p[j];

      tp[i] = 1 - (p[n - 1] - p[i]) / cp;
      for ( j = i + 1; j < n; j++ ) {
        tp[j] = (p[j] - p[j - 1]) / cp;
      }

    } else if ( c < n - 1 ) {

      tp[c + 1] = 1;

    } else { /* c == n - 1 */

      cp = 0;
      for ( j = 0; j < i; j++ ) {
        tp[j] = p[j] / p[n - 1];
        cp += tp[j];
      }
      tp[i] = 1 - cp;

    }
  } else {
    tc_getprobarr0(n, c, p, tp);
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



/* select using the metropolis way */
__inline static
int metropolis_select(int n, int c, const double *p)
{
  int r, acc = 0;

  r = ( c + 1 + (int) ( (n - 1) * rand01() ) ) % n;
  if ( p[r] >= p[c] ) {
    acc = 1;
  } else if ( rand01() < p[r] / p[c] ) {
    acc = 1;
  }

  return acc ? r : c;
}



/* select a state according to
 * the type I convective transition matrix
 * p has been sorted in ascending order */
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



/* choose an item according to the
 * type II convective transition matrix
 * p has been sorted in ascending order */
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



/* sort the probabilities `p` in ascending order
 * save the sorted array in `ps`
 * if `cnt` exists, the sorted version is saved in `scnt`
 * return the index of `c` */
__inline static
int tc_sort(int n, int c, const double *p,
    int *idmap, double *ps,
    double *cnt, double *scnt)
{
  int i, j, jm, ic = c;
  double x, y;

  for ( i = 0; i < n; i++ ) idmap[i] = i;

  /* bubble sort */
  for ( i = 0; i < n; i++ ) {
    /* find the ith smallest item */
    jm = i;
    x = p[ idmap[jm] ];
    for ( j = i + 1; j < n; j++ ) {
      y = p[ idmap[j] ];
      if ( y < x ) {
        x = y;
        jm = j;
      }
    }
    /* handling index
     * idmap[i] is the index of the ith smallest item */
    if ( jm != i ) {
      /* swap idmap[i] and idmap[jm] */
      j = idmap[i];
      idmap[i] = idmap[jm];
      idmap[jm] = j;
    }
    ps[i] = p[ idmap[i] ];
    if ( cnt != NULL ) {
      scnt[i] = cnt[ idmap[i] ];
    }
    if ( idmap[i] == c ) {
      ic = i;
    }
  }

  return ic;
}



/* randomly select an item from an unsorted array */
__inline static
int tc_select(int type, int n, int c, const double *p,
    int *idmap, double *ps, double *cp)
{
  int i, ic = c;

  type &= TC_MATRIXMASK;

  if ( type == TC_HEATBATH ) {

    return heatbath_select(n, p, cp);

  } else if ( type == TC_METROPOLIS ) {

    return metropolis_select(n, c, p);

  }

  /* 1. sort the probabilities p in ascending order */
  ic = tc_sort(n, c, p, idmap, ps, NULL, NULL);

  //for ( i = 0; i < n; i++ ) printf("%d: %8.5f  %d %8.5f\n", i, p[i], idmap[i], ps[i]); getchar();

  /* 2. select from the sorted probability */
  if ( type == TC_MATRIX1 ) {
    i = tc1_select(n, ic, ps, cp);
  } else if ( type == TC_MATRIX2 ) {
    i = tc2_select(n, ic, ps, cp);
  }

  /* 3. map back to the original index */
  return idmap[i];
}



/* deterministically select an item from an unsorted array */
__inline static
int tc_next(int type, int n, int c, const double *p,
    int *idmap, double *ps, double *tp,
    double *cnt, double *scnt)
{
  int i, ii, j, ic;

  type &= TC_MATRIXMASK;

  if ( type == TC_HEATBATH ) {
    ii = fsamp_select(n, p, cnt);
    cnt[ii] += 1;
    fsamp_truncate(n, p, cnt);
    return ii;
  }

  /* 1. sort the probabilities p in ascending order */
  ic = tc_sort(n, c, p, idmap, ps, cnt, scnt);

  /* 2. get the array of transition probabilities */
  if ( type == TC_MATRIX1 ) {
    tc1_getprobarr(n, ic, ps, tp);
  } else if ( type == TC_MATRIX2 ) {
    tc2_getprobarr(n, ic, ps, tp);
  } else if ( type == TC_METROPOLIS ) {
    for ( j = 0; j < n; j++ ) {
      tp[j] = ps[j];
    }
  }

  /* 3. select deterministically */
  i = fsamp_select(n, tp, scnt);
  ii = idmap[i];
  fsamp_truncate(n, tp, scnt);

  /* 4. convert the counts back to the index */
  for ( j = 0; j < n; j++ ) {
    cnt[ idmap[j] ] = scnt[j];
  }

  cnt[ ii ] += 1;

  return ii;
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



__inline static
void tc_printmat(int type, int n, const double *p, const char *name)
{
  int r, c, casenum = 0;
  double *tmat, *tcol;

  xnew(tmat, n * n);
  xnew(tcol, n);

  if ( type == TC_MATRIX1 ) {
    casenum = tc1_getcase(n, p);
  } else if ( type == TC_MATRIX2 ) {
    casenum = tc2_getcase(n, p);
  }

  for ( c = 0; c < n; c++ ) {
    if ( type == TC_MATRIX1 ) {
      tc1_getprobarr(n, c, p, tcol);
    } else if ( type == TC_MATRIX2 ) {
      tc2_getprobarr(n, c, p, tcol);
    }
    for ( r = 0; r < n; r++ ) {
      tmat[r*n + c] = tcol[r];
    }
  }
/*
  for ( r = 0; r < n; r++ ) {
    for ( c = 0; c < n; c++ ) {
      if ( type == TC_MATRIX1 ) {
        tmat[r*n+c] = tc1_getprob(n, r, c, p);
      } else if ( type == TC_MATRIX2 ) {
        tmat[r*n+c] = tc2_getprob(n, r, c, p);
      }
    }
  }
*/

  tc_pmat(n, tmat, p, name, casenum);
  free(tmat);
  free(tcol);
}



#endif /* TCONV_H__ */


#include <math.h>
#include <float.h>

/* Numerical recipes
 * http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-9.pdf
 * https://searchcode.com/codesearch/view/13427461/ */
__inline static void cisi(double x, double *ci, double *si)
{
  int i, k, odd;
  double a, err, fact, sign, sum, sumc, sums, t, term;
  double h[2], b[2], c[2], d[2], del[2], y[2], yy;
  const double BIG = DBL_MAX*DBL_EPSILON;
  const double EULER = 0.577215664901532860606512090;

  t = fabs(x);
  if ( t <= 0 ) {
    *si = 0;
    *ci = -BIG;
    return;
  }

  if ( t > 2 ) {
    /* use the continued fraction */
    b[0] = 1;
    b[1] = t;
    c[0] = BIG;
    c[1] = 0;
    yy = 1 + t*t;
    d[0] = h[0] =  1/yy;
    d[1] = h[1] = -t/yy;
    for ( i = 1; i < 100; i++ ) {
      a = -i*i;
      b[0] += 2;
      /* d = 1/(a*d+b) */
      y[0] = a*d[0] + b[0];
      y[1] = a*d[1] + b[1];
      yy = y[0]*y[0] + y[1]*y[1];
      d[0] =  y[0]/yy;
      d[1] = -y[1]/yy;
      /* c = b + a/c
       * for i = 1, c is infinity, so need to be careful */
      yy = c[0]*c[0] + c[1]*c[1];
      c[0] = b[0] + a*c[0]/yy;
      c[1] = b[1] - a*c[1]/yy;
      /* del = c * d */
      del[0] = c[0] * d[0] - c[1] * d[1];
      del[1] = c[0] * d[1] + c[1] * d[0];
      /* h *= del */
      y[0] = h[0] * del[0] - h[1] * del[1];
      h[1] = h[0] * del[1] + h[1] * del[0];
      h[0] = y[0];
      if ( fabs(del[0] - 1) + fabs(del[1]) < DBL_EPSILON )
        break;
      //printf("yy %g, c %g + %gI, del %g + %gI\n", yy, c[0], c[1], del[0], del[1]);
    }
    /* h *= e^{-it} */
    c[0] = cos(t);
    c[1] = -sin(t);
    y[0] = h[0] * c[0] - h[1] * c[1];
    h[1] = h[1] * c[0] + h[0] * c[1];
    h[0] = y[0];
    *ci = -h[0];
    *si = M_PI/2 + h[1];
  } else {
    /* power series */
    if ( t < 2 * sqrt(DBL_MIN) ) {
      sumc = 0;
      sums = t;
    } else {
      sum = sums = sumc = 0.0;
      sign = fact = 1.0;
      odd = 1;
      for ( k = 1; k <= 100; k++ ) {
        fact *= t / k;
        term = fact / k;
        sum += sign * term;
        err = term / fabs(sum);
        if ( odd ) {
          sign = -sign;
          sums = sum;
          sum = sumc;
        } else {
          sumc = sum;
          sum = sums;
        }
        if ( err < DBL_EPSILON ) break;
        odd = !odd;
      }
    }
    *si = sums;
    *ci = sumc + log(t) + EULER;
  }
  if ( x < 0 ) *si = -(*si);
}


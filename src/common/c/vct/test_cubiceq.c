#include <stdio.h>
#include <math.h>

/* compute the roots of a cubic equation
 *
 * x^3 - 3*p*x - 2*q = 0
 *
 * The roots are arranged in descending order
 *
 * The code in this function is partly in meigval() of mat.h
 * */
__inline static double *solve_cubic_eq(double p, double q, double v[3])
{
  const double m = 0;
  double discr = q*q - p*p*p;

  if (discr <= 0 || p <= 0) { /* the p <= 0 case actually means p = q = 0 */
    /* try to solve the equation with
     * x = 2 sqrt(p) cos(t)
     *
     * 4 cos^3(t) - 3 cos(t) = q/p^(3/2) = cos(3 t)
     *
     * This solution holds if |q/p^(3/2)| <= 1, or if |p|^3 >= q^2 */
    double pr = sqrt(p);
    double cos3t = q / (pr * p);
    double t;

    /* currently, t holds the value of cos(3 t)
     * Its absolute value should not exceed 1.0,
     * but we allow some margin of error that is typical
     * in degenerate cases
     *
     * For example:
     *    x^3 - 3 x + 2 = (x + 2) (x - 1)^2 = 0
     * or
     *    x^3 - 3 x - 2 = (x - 2) (x + 1)^2 = 0
     * */
    if (cos3t > 1) cos3t = 1;
    else if (cos3t < -1) cos3t = -1;

    t = acos(cos3t) / 3; /* 0 < t < pi/3 */
    v[0] = m + 2.0 * pr * cos(t);  /* largest */
    v[1] = m + 2.0 * pr * cos(t - 2*M_PI/3); /* second largest */
    v[2] = m + 2.0 * pr * cos(t + 2*M_PI/3); /* smallest */

  } else {

    /* Cardano's formula */
    double del = sqrt(discr);

    /* use cbrt(z) instead of pow(z, 1.0/3) as
     * the latter would fail for a negative z */
    v[2] = m + cbrt(q + del);
    v[2] += cbrt(q - del);
    v[0] = v[1] = v[2];

  }

  return v;
}


void test_cubic_eq(double p, double q, double x0, double x1, double x2)
{
  double x[3];
  const double tol = 1e-5;

  solve_cubic_eq(p, q, x);

  if (fabs(x[0] - x0) > tol
   || fabs(x[1] - x1) > tol
   || fabs(x[2] - x2) > tol) {
    fprintf(stderr, "Error: failed to solve x^3 %+g x %+g = 0\n",
      -3*p, -2*q);
  } else {
    fprintf(stderr, "Equation x^3 %+g x %+g = 0 solved successfully\n",
      -3*p, -2*q);
  }
}

int main(void)
{
  test_cubic_eq(2, 3, 2.84732, 2.84732, 2.84732);
  test_cubic_eq(2, 1, 2.60168, -0.339877, -2.2618);
  test_cubic_eq(1, 1, 2, -1, -1);
  test_cubic_eq(-2, 1, 0.32748, 0.32748, 0.32748);
  test_cubic_eq(-1, -1, -0.596072, -0.596072, -0.596072);
  test_cubic_eq(4, -8, 2, 2, -4);

  return 0;
}


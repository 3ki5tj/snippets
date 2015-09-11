#include "mtrand.h" /* Mersenne-Twister random number generator */



int N = 200000; /* sample size */
int M = 100; /* number of bootstrap samples */
double tcorr = 9.5; /* autocorrelation time */
int verbose = 0;



/* generate a random sample with correlation time tau */
static void gensample(double *x, int n, double tau)
{
  int i;
  double gam = tau > 0 ? exp(-1/tau) : 0;
  double kappa = sqrt(1 - gam * gam);

  x[0] = randgaus(); /* generate a normal random number */

  for ( i = 1; i < n; i++ ) {
    x[i] = x[i-1] * gam + kappa * randgaus();
  }
}



/* compute the average of `f(x)` */
static double getave(const double *x, int n,
    double (*f)(double))
{
  int i;
  double ave = 0;

  for ( i = 0; i < n; i++ )
    ave += f(x[i]);
  return ave /= n;
}



/* estimate the autocorrelation time of `x`
 * truncate the autocorrelation function
 * when it falls under `cutoff`  */
static double gettcorr(const double *x, int n, double cutoff)
{
  int i, j;
  double tau = 0.5, x1, x2, s, ave, sig;

  /* 1. compute the average and variance */
  ave = sig = 0;
  for ( i = 0; i < n; i++ ) {
    ave += x[i];
    sig += x[i] * x[i];
  }
  ave /= n;
  sig = sqrt(sig / n - ave * ave);

  /* 2. compute the integral correlation time */
  for ( i = 1; i < n - 1; i++ ) {
    /* compute the autocorrelation function at separation i */
    s = 0;
    for ( j = 0; j < n - i; j++ ) {
      x1 = (x[j] - ave) / sig;
      x2 = (x[j + i] - ave) / sig;
      s += x1 * x2;
    }
    s /= (n - i);
    if ( s < cutoff ) break;
    /* integrate the autocorrelation function */
    tau += s;
  }

  return tau;
}



/* generate a bootstrap sample and compute the average of f(x) */
static double bootstrap(const double *x, int n,
    double tau, double (*f)(double))
{
  int i, j = 0;
  double xi, s = 0;
  double gam = tau > 0 ? exp(-1/tau) : 0;

  for ( i = 0; i < n; i++ ) {
    /* replace the current trajectory frame with probability `gam` */
    if ( rand01() >= gam || i == 0 ) {
      j = (int) (rand01() * n);
    }
    xi = x[j];
    s += (*f)(xi);
  }
  return s / n;
}



/* generate a bootstrap sample and compute the average of f(x)
 * using the Hub 2010 method
 * `arr` should be sorted */
static double bootstrap_hub(const double *arr, int n,
    double tau, double (*f)(double))
{
  int i, j = 0;
  double xi, s = 0;
  double gam = tau > 0 ? exp(-1/tau) : 0;
  double mag = sqrt(1 - gam * gam);
  double x, y;

  x = randgaus();
  for ( i = 0; i < n; i++ ) {
    /* update `x` with autocorrelation time `tau` */
    x = x * gam + mag * randgaus();
    /* map the normally distributed function to (0, 1) */
    y = (erf(x/sqrt(2)) + 1) / 2;
    /* compute the frame corresponding to y */
    j = (int) (y * n);

    xi = arr[j];
    s += (*f)(xi);
  }
  return s / n;
}



static int dblcmp(const void *pa, const void *pb)
{
  double a = *((double *) pa);
  double b = *((double *) pb);
  if ( a > b ) {
    return 1;
  } else if ( a < b ) {
    return -1;
  } else {
    return 0;
  }
}



/* use bootstrapping to estimate the error of f(x) */
static double getbserr(const double *x, int n, int m,
    double tau, double (*f)(double), int hub)
{
  int i;
  double y, ave = 0, sig = 0;
  double *xs; /* sorted array to facilitate binary search */

  xs = calloc(n, sizeof(*xs));
  for ( i = 0; i < n; i++ ) xs[i] = x[i];
  qsort(xs, n, sizeof(*xs), &dblcmp);

  /* generate M samples */
  for ( i = 0; i < m; i++ ) {
    if ( hub ) {
      y = bootstrap_hub(xs, n, tau, f);
    } else {
      y = bootstrap(x, n, tau, f);
    }
    if ( verbose ) printf("bootstrap %d: %g\n", i, y);
    ave += y;
    sig += y * y;
  }
  ave /= m;
  sig = sqrt((sig - m * ave * ave) / (m - 1));
  free(xs);
  return sig;
}



static double fx(double x)
{
  return 2 * x;
}



static double fexp(double x)
{
  return exp(2 * x);
}



int main(int argc, char **argv)
{
  double *x, ave, err, errhub, tau;

  if ( argc >= 2 ) N = atoi(argv[1]);
  if ( argc >= 3 ) M = atoi(argv[2]);
  if ( argc >= 4 ) tcorr = atof(argv[3]);
  printf("N %d, M %d, tcorr %g\n", N, M, tcorr);


  /* generate a test sample */
  x = calloc(N, sizeof(double));
  gensample(x, N, tcorr);


  /* estimate the integral correlation time */
  tau = gettcorr(x, N, 0.01);
  printf("integral correlation time %g vs %g\n\n", tau, tcorr);


  /* compute the average of f(x) */
  ave = getave(x, N, fx);
  /* use bootstrapping to estimate the error of f(x) */
  err = getbserr(x, N, M, tau, fx, 0);
  errhub = getbserr(x, N, M, tau, fx, 1);
  printf("1. <2*x> %g (should be 0), err %g (simple) / %g (Hub) (should be %g)\n\n",
      ave, err, errhub, 2*sqrt((1+tcorr*2)/N));


  /* try the same for the exponetial function */
  ave = getave(x, N, fexp);
  err = getbserr(x, N, M, tau, fexp, 0);
  errhub = getbserr(x, N, M, tau, fexp, 1);
  printf("2. <exp(2*x)> %g (should be %g), err %g (simple) / %g (Hub)\n",
      ave, exp(2), err, errhub);
  printf("  log<exp(2*x)> %g (should be 2), err %g (simple) / %g (Hub), cf. case 1\n\n",
      log(ave), err/ave, errhub/ave);


  free(x);
  return 0;
}


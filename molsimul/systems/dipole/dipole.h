#ifndef DIPOLE_H__
#define DIPOLE_H__

#include "vct.h"
#include "util.h"
#include "mtrand.h"

const double Kelec = 1;

typedef struct {
  int n;
  double l;
  double disp;
  double *q;
  double (*x)[3];
  double (*x1)[3];
  double (*v)[3];
  double (*f)[3];
  double (*cs)[2];
  double sig;
  double eps;
  double alpha;
  int km;
  double ene_lj;
  double ene_el;
  double ene;
} dipole_t;

static dipole_t *dipole_open(int n, double l, double disp,
    double sig, double eps, double alpha, int km)
{
  dipole_t *dp;
  int i, ix, iy, iz, nx;

  xnew(dp, 1);
  n = (n/2) * 2;
  dp->n = n;
  dp->l = l;
  dp->disp = disp;
  dp->sig = sig;
  dp->eps = eps;
  dp->alpha = alpha;
  dp->km = km;
  xnew(dp->q, n);
  xnew(dp->x, n);
  xnew(dp->x1, n);
  xnew(dp->v, n);
  xnew(dp->f, n);
  xnew(dp->cs, n);
  for ( i = 0; i < n; i++ ) {
    dp->q[i] = (i % 2 == 0 ? 1 : -1);
    vzero(dp->f[i]);
    vzero(dp->v[i]);
  }
  nx = (int) (pow(n*0.5, 1./3) + 0.999999);
  i = 0;
  for ( ix = 0; ix < nx && i < n; ix++ ) {
    for ( iy = 0; iy < nx && i < n; iy++ ) {
      for ( iz = 0; iz < nx; iz++ ) {
        dp->x[i][0] = ix * l / nx;
        dp->x[i][1] = iy * l / nx;
        dp->x[i][2] = iz * l / nx;
        dp->x[i+1][0] = dp->x[i][0];
        dp->x[i+1][1] = dp->x[i][1];
        dp->x[i+1][2] = dp->x[i][2] + disp;
        i += 2;
        if ( i > n ) break;
      }
    }
  }
  dp->ene_lj = 0;
  dp->ene_el = 0;
  dp->ene = 0;
  return dp;
}

static int dipole_savepos(dipole_t *dp, const char *fn)
{
  int i, n = dp->n;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d\n", n);
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%g %g %g\n", dp->x[i][0], dp->x[i][1], dp->x[i][2]);
    if ( i % 2 == 1 ) fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

static void dipole_close(dipole_t *dp)
{
  free(dp->q);
  free(dp->x);
  free(dp->x1);
  free(dp->v);
  free(dp->f);
  free(dp->cs);
  free(dp);
}

static double dist2_pbc(double *a, double *b, double *dx, double l)
{
  double x;
  int i;
  for (i = 0; i < 3; i++ ) {
    x = (a[i] - b[i]) / l + 100;
    dx[i] = (x - (int) (x + 0.5)) * l;
  }
  return vsqr(dx);
}

/* NOTE: zero f before calling */
static double ewald(int n, double (*x)[3], double (*f)[3],
    const double *q, double (*cs)[2], double K,
    double boxl, double alpha, int km)
{
  int i, j, l, m;
  double r, r2, eself, ereal = 0, erecip = 0;
  double dx[3], ene, fs, qij, ivol, xp;
  double sqrta = sqrt(alpha), k[3], k2, re, im, phase;
  double sqrtpi = sqrt(M_PI);

  /* real-space sum */
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      /* loop over possible periodic images */
      r2 = dist2_pbc(x[i], x[j], dx, boxl);
      r = sqrt(r2);
      qij = q[i] * q[j];
      ene = qij*erfc(sqrta*r)/r;
      ereal += K*ene;
      if ( f != NULL ) {
        fs = K*(qij*2*sqrta/sqrtpi*exp(-alpha*r*r) + ene)/r2;
        vsinc(f[i], dx,  fs);
        vsinc(f[j], dx, -fs);
      }
    }
  }

  ivol = 2*M_PI/(boxl*boxl*boxl);
  /* reciprocal space sum */
  for ( i = -km; i <= km; i++ ) {
    k[0] = i*2*M_PI/boxl;
    for ( j = -km; j <= km; j++ ) {
      k[1] = j*2*M_PI/boxl;
      for ( l = -km; l <= km; l++ ) {
        if ( i == 0 && j == 0 && l == 0) continue;
        k[2] = l*2*M_PI/boxl;
        k2 = vdot(k, k);
        re = im = 0;
        for ( m = 0; m < n; m++ ) {
          phase = vdot(k, x[m]);
          cs[m][0] = q[m] * cos(phase);
          cs[m][1] = q[m] * sin(phase);
          re += cs[m][0];
          im += cs[m][1];
        }
        xp = K*ivol*exp(-k2/4/alpha)/k2;
        erecip += (re*re + im*im) * xp;
        if ( f != NULL ) {
          for ( m = 0; m < n; m++ ) {
            fs = xp*2*(re*cs[m][1]-im*cs[m][0]);
            vsinc(f[m], k, fs);
          }
        }
      }
    }
  }

  /* self energy */
  for ( eself = 0, i = 0; i < n; i++ )
    eself += q[i] * q[i];
  eself *= -K*sqrt(alpha/M_PI);
  return eself + ereal + erecip;
}

/* compute the energy and force */
static double dipole_force(dipole_t *dp,
    double (*x)[3], double (*f)[3],
    double *e_lj, double *e_el)
{
  int i, j, n = dp->n;
  double dx[3], r2, ir2, ir6, rc2, fs, sig2;
  double ene_lj = 0, ene_el = 0;

  if ( f != NULL ) {
    for ( i = 0; i < n; i++ )
      vzero(f[i]);
  }

  rc2 = dp->l * 0.5;
  rc2 *= rc2;
  sig2 = dp->sig * dp->sig;
  for ( i = 0; i < n; i++ ) {
    /* skip atoms in the same dipole */
    for ( j = i + 2; j < n; j++ ) {
      r2 = dist2_pbc(x[i], x[j], dx, dp->l);
      if ( r2 > rc2 ) continue;
      ir2 = sig2 / r2;
      ir6 = ir2 * ir2 * ir2;
      ene_lj += 4 * ir6 * (ir6 - 1);
      if ( f != NULL ) {
        fs = ir6 * (48 * ir6 - 24) / r2;
        vsinc(f[i], dx,  fs);
        vsinc(f[j], dx, -fs);
      }
    }
  }

  ene_el = ewald(dp->n, x, f, dp->q, dp->cs,
      Kelec, dp->l, dp->alpha, dp->km);
  if ( e_lj != NULL ) *e_lj = ene_lj;
  if ( e_el != NULL ) *e_el = ene_el;
  return ene_lj + ene_el;
}

/* test if the force matches the energy */
static double dipole_testforce(dipole_t *dp)
{
  double ene1, ene2, s, del = 0.0001;
  int i, n = dp->n;

  ene1 = dipole_force(dp, dp->x, dp->f, NULL, NULL);
  for ( s = 0, i = 0; i < n; i++ )
    s += vsqr(dp->f[i]);
  s = del / s;
  for ( i = 0; i < n; i++ ) {
    vsadd(dp->v[i], dp->x[i], dp->f[i], s);
  }
  ene2 = dipole_force(dp, dp->v, NULL, NULL, NULL);
  fprintf(stderr, "ene %g -> %g, de %g\n", ene1, ene2, (ene1-ene2)/del);
}

/* translational move */
static int dipole_trans(dipole_t *dp, double amp)
{
  int i, j, n = dp->n;
  double dx;

  i = (int) (rand01() * (n/2)) * 2;
  for ( j = 0; j < n; j++ )
    vcopy(dp->x1[j], dp->x[j]);
  /* move i and i+1 */
  for ( j = 0; j < D; j++ ) {
    double dx = amp * (2 * rand01() - 1);
    dp->x1[i][j] += dx;
    dp->x1[i+1][j] += dx;
  }
  return i;
}

/* Metropolis move */
static int dipole_metro(dipole_t *dp, double amp, double beta)
{
  double r, ene0, ene1, ene_lj1, ene_el1;
  int acc, i;

  i = dipole_trans(dp, amp);
  ene0 = dp->ene;
  ene1 = dipole_force(dp, dp->x1, NULL, &ene_lj1, &ene_el1);
  acc = 0;
  if ( ene1 < ene0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp(-beta*(ene1 - ene0)) );
  }
  if ( acc ) {
    vcopy(dp->x[i], dp->x1[i]);
    vcopy(dp->x[i+1], dp->x1[i+1]);
    /* copy force */
    dp->ene_lj = ene_lj1;
    dp->ene_el = ene_el1;
    dp->ene = ene1;
  }
  return acc;
}

#endif /* DIPOLE_H__ */

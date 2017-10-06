#ifndef WAT_H__
#define WAT_H__

#include "vct.h"
#include "util.h"
#include "mtrand.h"

#define Q_O   -0.834
#define Q_H    0.417
#define R_OH   0.9572
#define A_HOH  (104.52/180*M_PI)
/* LJ parameters from CHARMM */
#define O_EPS  0.1521
#define O_SIG  1.7682
#define H_EPS  0.046
#define H_SIG  0.2245


/* LJ parameter for OO */
//#define LJ_E12 582.0e3
//#define LJ_E6  595.0

const double Kelec = 332.0636; /* A kcal/mol, same as COULOMB in common.h */


typedef struct {
  int nw;
  int n;
  double l;
  double *q;
  double (*x)[3];
  double (*x1)[3];
  double (*v)[3];
  double (*f)[3];
  double (*cs)[2];
  double eps[3][3], sig[3][3];
  double alpha;
  int km;
  double ene_lj;
  double ene_el;
  double ene;
} wat_t;

static wat_t *wat_open(int nw, double l,
    double alpha, int km)
{
  wat_t *w;
  int i, j, ix, iy, iz, nx, n;
  double dx, dy;

  xnew(w, 1);
  n = nw * 3;
  w->nw = nw;
  w->n = n;
  w->l = l;
  w->alpha = alpha;
  w->km = km;
  xnew(w->q, n);
  xnew(w->x, n);
  xnew(w->x1, n);
  xnew(w->v, n);
  xnew(w->f, n);
  xnew(w->cs, n);
  for ( i = 0; i < n; i++ ) {
    w->q[i] = (i % 3 == 0 ? Q_O : Q_H);
    vzero(w->f[i]);
    vzero(w->v[i]);
  }
  nx = (int) (pow(n/3, 1./3) + 0.999999);
  i = 0;
  dx = R_OH * cos(A_HOH/2);
  dy = R_OH * sin(A_HOH/2);
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      if ( i == 0 && j == 0 ) {
        w->eps[i][j] = O_EPS;
        w->sig[i][j] = O_SIG * 2;
      } else if ( i > 0 && j > 0 ) {
        w->eps[i][j] = H_EPS;
        w->sig[i][j] = H_SIG * 2;
      } else {
        w->eps[i][j] = sqrt(O_EPS * H_EPS);
        w->sig[i][j] = O_SIG + H_SIG;
      }
    }
  }
  for ( ix = 0; ix < nx && i/3 < nw; ix++ ) {
    for ( iy = 0; iy < nx && i/3 < nw; iy++ ) {
      for ( iz = 0; iz < nx; iz++ ) {
        w->x[i][0] = ix * l / nx;
        w->x[i][1] = iy * l / nx;
        w->x[i][2] = iz * l / nx;
        w->x[i+1][0] = w->x[i][0] + dx;
        w->x[i+1][1] = w->x[i][1] + dy;
        w->x[i+1][2] = w->x[i][2];
        w->x[i+2][0] = w->x[i][0] - dx;
        w->x[i+2][1] = w->x[i][1] + dy;
        w->x[i+2][2] = w->x[i][2];
        i += 3;
        if ( i >= n ) break;
      }
    }
  }
  w->ene_lj = 0;
  w->ene_el = 0;
  w->ene = 0;
  return w;
}

static int wat_loadpdb(wat_t *w, const char *fn)
{
  FILE *fp;
  char buf[1024] = "";
  int i = 0;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  while ( fgets(buf, sizeof buf, fp) ) {
    if ( strncmp(buf, "ATOM  ", 6) != 0 ) continue;
    sscanf(buf+30, "%lf%lf%lf", &w->x[i][0], &w->x[i][1], &w->x[i][2]);
    //printf("%d %g %g %g\n", i, w->x[i][0], w->x[i][1], w->x[i][2]);
    i += 1;
  }
  fclose(fp);
  return 0;
}

static int wat_savepos(wat_t *w, const char *fn)
{
  int i, n = w->n;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d\n", n);
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%g %g %g\n", w->x[i][0], w->x[i][1], w->x[i][2]);
    if ( i % 3 == 2 ) fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

static void wat_close(wat_t *w)
{
  free(w->q);
  free(w->x);
  free(w->x1);
  free(w->v);
  free(w->f);
  free(w->cs);
  free(w);
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
static double wat_force(wat_t *w,
    double (*x)[3], double (*f)[3],
    double *e_lj, double *e_el)
{
  int i, j, n = w->n;
  double dx[3], r2, ir2, ir6, rc2, fs, sig, eps;
  double ene_lj = 0, ene_el = 0, ene_el0 = 0;

  if ( f != NULL ) {
    for ( i = 0; i < n; i++ )
      vzero(f[i]);
  }

  rc2 = w->l * 0.5;
  rc2 *= rc2;
  for ( i = 0; i < n - 3; i ++ ) {
    /* skip atoms in the same dipole */
    for ( j = (i/3 + 1) * 3; j < n; j ++ ) {
      r2 = dist2_pbc(x[i], x[j], dx, w->l);
      if ( r2 > rc2 ) continue;
      sig = w->sig[i%3][j%3];
      eps = w->eps[i%3][j%3];
      //printf("i %d, j %d, sig %g, eps %g\n", i, j, sig, eps);
      ir2 = sig * sig / r2;
      ir6 = ir2 * ir2 * ir2;
      //if ( LJ_E12*ir6*ir6 > 10 ) {
      //  printf("i %d, j %d, r %g, LJ12 %g, LJ6 %g\n", i, j, sqrt(r2), LJ_E12*ir6*ir6, LJ_E6*ir6); getchar();
      //}
      ene_lj += eps * (ir6 - 2) * ir6;
      if ( f != NULL ) {
        fs = eps * ir6 * 12 * (ir6 - 1) / r2;
        vsinc(f[i], dx,  fs);
        vsinc(f[j], dx, -fs);
      }
    }
  }

  ene_el0 = 0;
  for ( i = 0; i < 2; i++ ) {
    for ( j = i + 1; j < 3; j++ ) {
      r2 = dist2_pbc(x[i], x[j], dx, w->l);
      ene_el0 += Kelec*w->q[i]*w->q[j]/sqrt(r2);
      //printf("%d %d %g\n", i, j, ene_el0);
    }
  }
  ene_el0 *= w->nw;

  ene_el = ewald(w->n, x, f, w->q, w->cs,
      Kelec, w->l, w->alpha, w->km);
  ene_el -= ene_el0;

  if ( e_lj != NULL ) *e_lj = ene_lj;
  if ( e_el != NULL ) *e_el = ene_el;
  return ene_lj + ene_el;
}

/* test if the force matches the energy */
static double wat_testforce(wat_t *w)
{
  double ene1, ene2, s, del = 0.0001;
  int i, n = w->n;

  ene1 = wat_force(w, w->x, w->f, NULL, NULL);
  for ( s = 0, i = 0; i < n; i++ )
    s += vsqr(w->f[i]);
  s = del / s;
  for ( i = 0; i < n; i++ ) {
    vsadd(w->v[i], w->x[i], w->f[i], s);
  }
  ene2 = wat_force(w, w->v, NULL, NULL, NULL);
  fprintf(stderr, "ene %g -> %g, de %g\n", ene1, ene2, (ene1-ene2)/del);
}

/* translational move */
static int wat_trans(wat_t *w, double amp)
{
  int i, j, n = w->n;
  double dx;

  i = (int) (rand01() * (n/2)) * 2;
  for ( j = 0; j < n; j++ )
    vcopy(w->x1[j], w->x[j]);
  /* move i and i+1 */
  for ( j = 0; j < D; j++ ) {
    double dx = amp * (2 * rand01() - 1);
    w->x1[i][j] += dx;
    w->x1[i+1][j] += dx;
  }
  return i;
}

/* Metropolis move */
static int wat_metro(wat_t *w, double amp, double beta)
{
  double r, ene0, ene1, ene_lj1, ene_el1;
  int acc, i;

  i = wat_trans(w, amp);
  ene0 = w->ene;
  ene1 = wat_force(w, w->x1, NULL, &ene_lj1, &ene_el1);
  acc = 0;
  if ( ene1 < ene0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp(-beta*(ene1 - ene0)) );
  }
  if ( acc ) {
    vcopy(w->x[i], w->x1[i]);
    vcopy(w->x[i+1], w->x1[i+1]);
    /* copy force */
    w->ene_lj = ene_lj1;
    w->ene_el = ene_el1;
    w->ene = ene1;
  }
  return acc;
}

#endif /* WAT_H__ */

#ifndef LJCORE_H__
#define LJCORE_H__



/* this file collect routines common to all dimensions
 * define the dimension D before including this file
 * Note: coordinates are not reduced */



#include "mtrand.h"
#include "ljutil.h"



typedef struct {
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  double rho;
  double l, vol;
  double rc2, rc;
  double rcdef; /* preferred cutoff */

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double epot;
  double ekin;
  double epot_shift;
  double epot_tail;
  double p_tail;
} lj_t;



/* the following functions are dimension D dependent */
static void lj_initfcc(lj_t *lj);
static double lj_gettail(lj_t *lj, double rho, int n, double *ptail);
static void lj_shiftang(double (*x)[D], double (*v)[D], int n);



/* set density and compute tail corrections */
static void lj_setrho(lj_t *lj, double rho)
{
  double irc;
  lj->rho = rho;
  lj->vol = lj->n/rho;
  lj->l = pow(lj->vol, 1./D);
  if ((lj->rc = lj->rcdef) > lj->l/2) lj->rc = lj->l/2;
  lj->rc2 = lj->rc * lj->rc;
  irc = 1/lj->rc;
  irc *= irc * irc;
  irc *= irc;
  lj->epot_shift = 4*irc*(irc - 1);
  lj->epot_tail = lj_gettail(lj, rho, lj->n, &lj->p_tail);
}



/* remove the center of mass motion */
static void lj_rmcom(double (*x)[D], int n)
{
  int i;
  double rc[D] = {0};

  for ( i = 0; i < n; i++ )
    vinc(rc, x[i]);
  vsmul(rc, 1./n);
  for ( i = 0; i < n; i++ )
    vdec(x[i], rc);
}



/* open an LJ system */
static lj_t *lj_open(int n, double rho, double rcdef)
{
  lj_t *lj;
  int i, d;

  xnew(lj, 1);
  lj->n = n;
  lj->dof = n * D - D * (D+1)/2;
  lj->rcdef = rcdef;

  xnew(lj->x, n);
  xnew(lj->v, n);
  xnew(lj->f, n);

  lj_setrho(lj, rho);

  lj_initfcc(lj);

  /* init. random velocities */
  for (i = 0; i < n; i++)
    for ( d = 0; d < D; d++ )
      lj->v[i][d] = randgaus();

  lj_rmcom(lj->v, lj->n);
  lj_shiftang(lj->x, lj->v, lj->n);

  return lj;
}



static void lj_close(lj_t *lj)
{
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj);
}



/* write positions (and possibly velocities) */
static int lj_writepos(lj_t *lj, double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d, n = lj->n;

  if (fn == NULL) fn = "lj.pos";
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %.14e\n", D, n, (v != NULL), lj->l);
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ )
      fprintf(fp, "%.14e ", x[i][d]);
    if ( v != NULL )
      for ( d = 0; d < D; d++ )
        fprintf(fp, "%.14e ", v[i][d]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



#define LJ_PBC(x, l, invl) { (x) -= ((int)((x)*invl + 1000.5) - 1000.)*l; }



static double *lj_vpbc(double *v, double l, double invl)
{
  int d;
  for ( d = 0; d < D; d++ )
    LJ_PBC(v[d], l, invl);
  return v;
}



static double lj_pbcdist2(double *dx, const double *a, const double *b,
    double l, double invl)
{
  return vsqr( lj_vpbc(vdiff(dx, a, b), l, invl) );
}



/* compute force and virial, return energy */
__inline static double lj_energy(lj_t *lj, double (*x)[D],
    double *virial, double *ep0, double *eps)
{
  double dx[D], dr2, dr6, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if (dr2 > rc2) continue;
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      vir += dr6 * (48 * dr6 - 24); /* f.r */
      ep += 4 * dr6 * (dr6 - 1);
      npr++;
    }
  }
  if (virial) *virial = vir;
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  return ep + lj->epot_tail; /* unshifted energy */
}



/* compute force and virial, return energy */
__inline static double lj_force(lj_t *lj, double (*x)[D], double (*f)[D],
    double *virial, double *ep0, double *eps)
{
  double dx[D], fi[D], dr2, dr6, fs, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if (dr2 > rc2) continue;
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48 * dr6 - 24); /* f.r */
      vir += fs; /* f.r */
      fs *= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += 4 * dr6 * (dr6 - 1);
      npr++;
    }
    vinc(f[i], fi);
  }
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  if (virial) *virial = vir;
  return ep + lj->epot_tail; /* unshifted energy */
}



/* velocity-verlet */
static void lj_vv(lj_t *lj, double dt)
{
  int i, n = lj->n;
  double dth = dt * .5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(lj->v[i], lj->f[i], dth);
    vsinc(lj->x[i], lj->v[i], dt);
  }
  lj->epot = lj_force(lj, lj->x, lj->f, NULL, NULL, NULL);
  for (i = 0; i < n; i++) /* VV part 2 */
    vsinc(lj->v[i], lj->f[i], dth);
}



/* compute the kinetic energy */
static double lj_ekin(double (*v)[D], int n)
{
  int i;
  double ek = 0;
  for ( i = 0; i < n; i++ ) ek += vsqr( v[i] );
  return ek/2;
}



#define lj_vrescale(lj, tp, dt) \
  lj_vrescale_low(lj->v, lj->n, lj->dof, tp, dt)

/* exact velocity rescaling thermostat */
static double lj_vrescale_low(double (*v)[D], int n, int dof, double tp, double dt)
{
  int i;
  double ek1, ek2, s, c, r, r2;

  c = (dt < 700) ? exp(-dt) : 0;
  ek1 = lj_ekin(v, n);
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp / 2);
  if (ek2 < 0) ek2 = 0;
  s = sqrt(ek2/ek1);
  for (i = 0; i < n; i++) vsmul(v[i], s);
  return ek2;
}



#endif /* LJCORE_H__ */

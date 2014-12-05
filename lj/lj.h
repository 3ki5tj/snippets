#ifndef LJ_H__
#define LJ_H__



/* This module is a simplification of the zcom module of the same name.
 * It is more self-contained.
 * Several differences are listed below:
 *
 *  o Dimension is 3
 *  o Coordinates are not reduced
 *
 * */



#define D 3
#include "util.h"
#include "mtrand.h"



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



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  double a;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l/n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        lj->x[id][0] = (i + .5) * a;
        lj->x[id][1] = (j + .5) * a;
        lj->x[id][2] = (k + .5) * a;
        id++;
      }
}



/* get the tail correction */
static double lj_gettail(lj_t *lj, double rho, int n, double *ptail)
{
  double irc, irc3, irc6, utail;

  irc = 1/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = ( 8*M_PI*rho*n/9*(irc6 - 3)*irc3 );
  if (ptail != NULL)
    *ptail = ( 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3 );
  return utail;
}



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
  double rc[D] = {0, 0, 0};

  for ( i = 0; i < n; i++ )
    vinc(rc, x[i]);
  vsmul(rc, 1./n);
  for ( i = 0; i < n; i++ )
    vdec(x[i], rc);
}



/* annihilate the angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
static void lj_shiftang(double (*x)[D], double (*v)[D], int n)
{
  int i;
  double xc[D] = {0, 0, 0}, xi[D], ang[D], am[D] = {0, 0, 0}, dv[D], mat[D][D], inv[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) vinc(xc, x[i]);
  vsmul(xc, 1.f/n);
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(ang, xi, v[i]);
    vinc(am, ang);
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  rm3_inv(inv, mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  /* ang is the solution of M^(-1) * I */
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(dv, ang, xi);
    vinc(v[i], dv);
  }
}



/* compute the kinetic energy */
static double lj_ekin(double (*v)[D], int n)
{
  int i;
  double ek = 0;
  for ( i = 0; i < n; i++ ) ek += vsqr( v[i] );
  return ek/2;
}



/* exact velocity rescaling thermostat */
static double lj_vrescale(double (*v)[D], int n, int dof, double tp, double dt)
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



static double lj_pbc(double x, double l)
{
  return x - ((int)(x/l + 1000.5) - 1000)*l;
}



static double *lj_vpbc(double *v, double l)
{
  int d;
  for ( d = 0; d < D; d++ )
    v[d] = lj_pbc(v[d], l);
  return v;
}



static double lj_pbcdist2(double *dx,
    const double *a, const double *b, double l)
{
  return vsqr( lj_vpbc(vdiff(dx, a, b), l) );
}



/* compute force and virial, return energy */
__inline static double lj_energy(lj_t *lj, double (*x)[D],
    double *virial, double *ep0, double *eps)
{
  double dx[D], dr2, dr6, ep, vir, l = lj->l, rc2 = lj->rc2;
  int i, j, npr = 0, n = lj->n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l);
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
  double dx[D], fi[D], dr2, dr6, fs, ep, vir, l = lj->l, rc2 = lj->rc2;
  int i, j, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l);
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



#endif /* LJ_H__ */

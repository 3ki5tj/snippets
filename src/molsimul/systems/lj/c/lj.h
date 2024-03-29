#ifndef LJ_H__
#define LJ_H__



/* define the dimension D before including this file
 * Note: coordinates are not reduced */



#include "rand/mtrand.h"
#include "vct/vct.h"
#include "vct/mat.h"
#include "mdutil.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



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
  double *r2ij; /* pair distances */
  double *r2i; /* pair distances from i */
  double epot, ep0, ep6, ep12, eps;
  double vir;
  double ekin;
  double epot_shift;
  double epot_tail;
  double p_tail;
} lj_t;



/* functions that are dimension D dependent */
#if D == 2



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++) {
      if ((i+j) % 2 != 0) continue;
      /* add some noise to prevent two atoms happened to
       * be separated by precisely some special cutoff distance,
       * which might be half of the box */
      lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
      lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
      id++;
    }
}



/* get the tail correction */
static double lj_gettail(double rc, double rho, int n, double *ptail)
{
  double irc, irc3, irc6, utail;

  irc = 1 / rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = M_PI*rho*n*(.4*irc6 - 1)*irc3*irc;
  if (ptail != NULL)
    *ptail = M_PI*rho*rho*(2.4*irc6 - 3)*irc3*irc;
  return utail;
}



#else /* D == 3 */



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        /* add some noise to prevent two atoms happened to
         * be separated by precisely some special cutoff distance,
         * which might be half of the box */
        lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][2] = (k + .5) * a + noise * (2*rand01() - 1);
        id++;
      }
}



/* get the tail correction */
static double lj_gettail(double rc, double rho, int n, double *ptail)
{
  double irc, irc3, irc6, utail;

  irc = 1 / rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = 8*M_PI*rho*n/9*(irc6 - 3)*irc3;
  if (ptail != NULL)
    *ptail = 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3;
  return utail;
}



#endif



/* set density and compute tail corrections */
static void lj_setrho(lj_t *lj, double rho)
{
  double irc;

  lj->rho = rho;
  lj->vol = lj->n/rho;
  lj->l = pow(lj->vol, 1.0 / D);
  lj->rc = lj->l * 0.5;
  if ( lj->rc > lj->rcdef ) {
    lj->rc = lj->rcdef;
  }
  lj->rc2 = lj->rc * lj->rc;
  irc = 1 / lj->rc;
  irc *= irc * irc;
  irc *= irc;
  lj->epot_shift = 4 * irc * (irc - 1);
  lj->epot_tail = lj_gettail(lj->rc, rho, lj->n, &lj->p_tail);
}



/* open an LJ system */
static lj_t *lj_open(int n, double rho, double rcdef, int dopr)
{
  lj_t *lj;
  int i, d;

  xnew(lj, 1);
  lj->n = n;
  lj->dof = n * D - D; // - D * (D+1)/2;
  lj->rcdef = rcdef;

  xnew(lj->x, n);
  xnew(lj->v, n);
  xnew(lj->f, n);

  if ( dopr ) {
    xnew(lj->r2ij, n * n);
    xnew(lj->r2i, n);
  } else {
    lj->r2ij = NULL;
    lj->r2i = NULL;
  }

  lj_setrho(lj, rho);

  lj_initfcc(lj);

  /* initialize random velocities */
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ ) {
      lj->v[i][d] = randgaus();
    }
  }

  rmcom(lj->v, NULL, n);
  shiftang(lj->x, lj->v, NULL, n);

  return lj;
}



/* close the lj object */
static void lj_close(lj_t *lj)
{
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj->r2ij);
  free(lj->r2i);
  free(lj);
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
  lj_vpbc(vdiff(dx, a, b), l, invl);
  return vsqr( dx );
}



#define lj_energy(lj) \
  lj->epot = lj_energy_low(lj, lj->x, lj->r2ij, \
      &lj->vir, &lj->ep0, &lj->ep6, &lj->ep12, &lj->eps)

/* compute force and virial, return energy */
__inline static double lj_energy_low(lj_t *lj, double (*x)[D],
    double *r2ij, double *virial, double *ep0,
    double *pep6, double *pep12, double *eps)
{
  double dx[D], dr2, ir6, ep, ep6, ep12, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  ep6 = ep12 = 0;
  for ( i = 0; i < n - 1; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) continue;
      dr2 = 1 / dr2;
      ir6 = dr2 * dr2 * dr2;
      ep12 += ir6 * ir6;
      ep6 += ir6;
      npr++;
    }
  }

  ep6 *= 4;
  ep12 *= 4;
  ep = ep12 - ep6;
  if ( ep0 ) *ep0 = ep;
  if ( eps ) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  if ( pep6 ) *pep6 = ep6;
  if ( pep12 ) *pep12 = ep12;
  if ( virial ) *virial = 12 * ep12 - 6 * ep6;

  return ep + lj->epot_tail; /* unshifted energy */
}



#define lj_force(lj) \
  lj->epot = lj_force_low(lj, lj->x, lj->f, \
      lj->r2ij, &lj->vir, &lj->ep0, \
      &lj->ep6, &lj->ep12, &lj->eps)

/* compute force and virial, return energy
 * the pair distances are recomputed */
__inline static double lj_force_low(lj_t *lj, double (*x)[D], double (*f)[D],
    double *r2ij, double *virial, double *ep0,
    double *pep6, double *pep12, double *eps)
{
  double dx[D], fi[D], dr2, ir6, fs;
  double ep, ep6, ep12, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }
  ep6 = ep12 = 0;
  for ( i = 0; i < n - 1; i++ ) {
    vzero(fi);
    for ( j = i + 1; j < n; j++ ) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) continue;
      dr2 = 1 / dr2;
      ir6 = dr2 * dr2 * dr2;
      fs = ir6 * (48 * ir6 - 24); /* f.r */
      fs *= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep6 += ir6;
      ep12 += ir6 * ir6;
      npr++;
    }
    vinc(f[i], fi);
  }

  ep6 *= 4;
  ep12 *= 4;
  ep = ep12 - ep6;
  if ( ep0 ) *ep0 = ep;
  if ( eps ) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  if ( pep6 ) *pep6 = ep6;
  if ( pep12 ) *pep12 = ep12;
  if ( virial ) *virial = 12 * ep12 - 6 * ep6;

  return ep + lj->epot_tail; /* unshifted energy */
}



/* compute pressure */
static double lj_calcp(lj_t *lj, double tp)
{
  return (lj->dof * tp + lj->vir) / (D * lj->vol) + lj->p_tail;
}



/* velocity-verlet */
__inline static void lj_vv(lj_t *lj, double dt)
{
  int i, n = lj->n;
  double dth = dt * .5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(lj->v[i], lj->f[i], dth);
    vsinc(lj->x[i], lj->v[i], dt);
  }
  lj_force(lj);
  for (i = 0; i < n; i++) /* VV part 2 */
    vsinc(lj->v[i], lj->f[i], dth);
}



/* compute the kinetic energy */
#define lj_ekin(v, n) md_ekin(v, NULL, n)



/* exact velocity-rescaling thermostat */
#define lj_vrescale(lj, tp, dt) \
  md_vrescale(lj->v, NULL, lj->n, lj->dof, tp, dt)



/* Nose-Hoover chain thermostat */
#define lj_nhchain(lj, tp, dt, nhclen, zeta, zmass) \
  md_nhchain(lj->v, NULL, lj->n, lj->dof, tp, dt, nhclen, zeta, zmass)



__inline static double lj_langevin(lj_t *lj, double tp, double dt)
{
  int n = lj->n;

  lj->dof = n * D;
  md_langevin(lj->v, NULL, n, tp, dt);
  return md_ekin(lj->v, NULL, n);
}




/* position Langevin barostat, with coordinates only
 * set cutoff to half of the box */
__inline static void lj_langp0(lj_t *lj, double dt,
    double tp, double pext, int ensx)
{
  double pint, amp, s, dlnv;
  int i;

  pint = lj_calcp(lj, tp);
  amp = sqrt(2 * dt);
  dlnv = ((pint - pext) * lj->vol / tp + 1 - ensx) * dt + amp * randgaus();
  s = exp( dlnv / D );
  lj->vol *= exp( dlnv );
  lj_setrho(lj, lj->n / lj->vol);
  for ( i = 0; i < lj->n; i++ ) {
    vsmul(lj->x[i], s);
  }
  lj_force(lj);
}



/* displace a random particle i, return i */
static int lj_randmv(lj_t *lj, double *xi, double amp)
{
  int i, d;

  i = (int) (rand01() * lj->n);
  for ( d = 0; d < D; d++ )
    xi[d] = lj->x[i][d] + (rand01() * 2 - 1) * amp;
  return i;
}



/* compute pair energy */
__inline static int lj_pair(double dr2,
    double rc2, double *u6, double *u12, double *vir)
{
  if ( dr2 < rc2 ) {
    double invr2 = 1 / dr2;
    double invr6 = invr2 * invr2 * invr2;
    *vir = invr6 * (48 * invr6 - 24); /* f.r */
    *u6  = 4 * invr6;
    *u12 = 4 * invr6 * invr6;
    return 1;
  } else {
    *vir = 0;
    *u6  = 0;
    *u12 = 0;
    return 0;
  }
}



/* return the energy change from displacing x[i] to xi */
__inline static double lj_depot(lj_t *lj, int i, double *xi,
    double *u6, double *u12, double *vir)
{
  int j, n = lj->n;
  double l = lj->l, invl = 1/l, rc2 = lj->rc2;
  double du6, du12, dvir;
  double dx[D], r2;

  *u6 = *u12 = *vir = 0.0;
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) {
      if ( lj->r2i != NULL )
        lj->r2i[j] = 0;
      continue;
    }
    if ( lj->r2ij != NULL ) {
      r2 = lj->r2ij[i*n + j];
    } else {
      r2 = lj_pbcdist2(dx, lj->x[i], lj->x[j], l, invl);
    }
    if ( lj_pair(r2, rc2, &du6, &du12, &dvir) ) {
      *u6 -= du6;
      *u12 -= du12;
      *vir -= dvir;
    }
    r2 = lj_pbcdist2(dx, xi, lj->x[j], l, invl);
    if ( lj_pair(r2, rc2, &du6, &du12, &dvir) ) {
      *u6 += du6;
      *u12 += du12;
      *vir += dvir;
    }
    if ( lj->r2i != NULL )
      lj->r2i[j] = r2;
  }
  return *u12 - *u6;
}



/* commit a particle displacement */
__inline static void lj_commit(lj_t *lj, int i, const double *xi,
    double du6, double du12, double dvir)
{
  int j, n = lj->n;
  double du = du12 - du6;

  vwrap( vcopy(lj->x[i], xi), lj->l );
  lj->ep6 += du6;
  lj->ep12 += du12;
  lj->ep0 += du;
  lj->epot += du;
  lj->vir += dvir;
  if ( lj->r2ij != NULL ) {
    for ( j = 0; j < n; j++ ) {
      lj->r2ij[i*n + j] = lj->r2i[j];
      lj->r2ij[j*n + i] = lj->r2i[j];
    }
  }
}



/* Metropolis algorithm */
__inline static int lj_metro(lj_t *lj, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0, du6 = 0, du12 = 0, dvir = 0;

  i = lj_randmv(lj, xi, amp);
  du = lj_depot(lj, i, xi, &du6, &du12, &dvir);
  if ( du < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -bet * du ) );
  }
  if ( acc ) {
    lj_commit(lj, i, xi, du6, du12, dvir);
    return 1;
  }
  return 0;
}



/* Monte Carlo volume move
 * the scaling is r = r*s, p = p/s;
 * set cutoff to half of the box
 * the ensemble distribution is
 *   distr. ~ V^(dof/D - ensx) exp(-beta * ep)
 * */
__inline static int lj_mcvmov(lj_t *lj, double lnvamp, double tp, double pext,
    int ensx)
{
  int acc = 0, i;
  double lnlo, lnln, lo, ln, vo, vn, s, ss, is6, dex;
  double epo, epn, ep6n, ep12n, eptailn, bet = 1/tp;

  vo = lj->vol;
  lo = lj->l;
  lnlo = log(lo);
  lnln = lnlo + lnvamp/D * (rand01() * 2 - 1);
  ln = exp( lnln );
  vn = exp( D * lnln );

  epo = lj->epot;
  s = ln / lo;
  ss = s * s;
  is6 = 1 / (ss * ss * ss);
  ep6n = lj->ep6 * is6;
  ep12n = lj->ep12 * is6 * is6;
  /* note: assuming half-box cutoff here */
  eptailn = lj_gettail(ln/2, lj->n / vn, lj->n, NULL);
  epn = ep12n - ep6n + eptailn;
  dex = bet * (epn - epo + pext * (vn - vo))
      + lj->dof * (lnlo - lnln)
      + D * (lnlo - lnln) * (1 - ensx);
  //printf("l %g -> %g, vol %g -> %g, s %g\n", lo, ln, vo, vn, s);
  //printf("ep %g, %g, %g -> %g, %g, %g\n", lj->ep6, lj->ep12, lj->epot, ep6n, ep12n, epn);

  acc = 1;
  if ( dex > 0 ) {
    double r = rand01();
    acc = ( r < exp(-dex) );
  }
  if ( acc ) { /* scale the velocities */
    lj_setrho(lj, lj->n / vn);
    lj->ep12 = ep12n;
    lj->ep6 = ep6n;
    lj->ep0 = ep12n - ep6n;
    lj->vir = ep12n * 12 - ep6n * 6;
    lj->epot = lj->ep0 + lj->epot_tail;
    /* scale the coordinates */
    for ( i = 0; i < lj->n; i++ ) {
      vsmul(lj->x[i], s);
    }
    if ( lj->r2ij != NULL ) {
      for ( i = 0; i < lj->n * lj->n; i++ ) {
        lj->r2ij[i] *= ss;
      }
    }
  }
  return acc;
}



/* wrap coordinates such that particles stay in the box */
__inline static int lj_wrapbox(lj_t *lj,
    double (*xin)[D], double (*xout)[D])
{
  int i, n = lj->n;
  double l = lj->l;

  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }
  return 0;
}



/* write positions (and possibly velocities) */
__inline static int lj_writepos(lj_t *lj,
    double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d, n = lj->n;

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



#endif /* LJ_H__ */

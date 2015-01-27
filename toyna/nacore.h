#ifndef NACORE_H__
#define NACORE_H__



#ifndef D
#define D 3
#endif



#include "mtrand.h"
#include "nautil.h"
#include "potential.h"
#include <ctype.h>



#define BOLTZK 0.0019872041 /* Boltzmann constant in kcal/mol/K */

#define TIS_IP 0
#define TIS_IS 1
#define TIS_IB 2

/* bond length parameters */
const double R0_PS = 4.247;
const double K0_PS = 64.0;
const double R0_SP = 3.878;
const double K0_SP = 23.0;
const double R0_SB[4] = {4.684, 4.083, 4.830, 4.086};
const double K0_SB = 10.0;

/* bond angle parameters */
const double A0_PSB[4] = {D2R(100.03), D2R(93.19),  D2R(104.57), D2R(92.99)};
const double A0_BSP[4] = {D2R(105.10), D2R(107.58), D2R(103.97), D2R(107.63)};
const double A0_PSP = D2R(88.01);
const double KA_PSB = 5.0;
const double KA_BSP = 5.0;
const double KA_PSP = 20.0;

/* WCA parameters */
#define WCA_SIG 3.2
const double WCA_SIG2 = WCA_SIG * WCA_SIG;
const double WCA_EPS = 1.0;



typedef struct {
  int nr; /* number of residues */
  char *seq; /* sequence */
  int *iseq; /* A: 0, C: 1, G: 2, U: 3 */

  int apr; /* number of atoms per residue */
  int na; /* number of atoms */
  int dof; /* degrees of freedom */
  double rc2, rc;

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double ekin;
  double epot;
} na_t;



/* initialize an RNA chain */
static void na_initchain2(na_t *na)
{
  int i, ib, is, nr = na->nr;
  double rb = 10.0, rs = 7.0, ang = 2*PI/10, dh = 3.4, th, c, s;

  for ( i = 0; i < nr; i++ ) {
    th = ang * i;
    c = cos(th);
    s = sin(th);
    ib = i*2;
    is = i*2 + 1;
    na->x[ib][0] = rb * c;
    na->x[ib][1] = rb * s;
    na->x[ib][2] = dh * i;
    na->x[is][0] = rs * c;
    na->x[is][1] = rs * s;
    na->x[is][2] = dh * i;
  }
}



/* initialize an RNA chain */
static void na_initchain3(na_t *na)
{
  int i, ic, nr = na->nr;
  double ang = D2R(32.7), dh = 2.81;
  const double rp = 8.710, thp = -D2R(70.502), dhp = 3.750;
  const double rs = 9.064, ths = -D2R(43.651), dhs = 2.806;
  const double rbarr[4] = {5.458, 5.643, 5.631, 5.633};
  const double thbarr[4] = {-D2R(25.976), -D2R(33.975), -D2R(22.124), -D2R(34.102)};
  const double dhbarr[4] = {0.742, 0.934, 0.704, 0.932};
  double rb, thb, dhb, th;

  for ( i = 0; i < nr; i++ ) {
    th = ang * i;

    na->x[i*3 + TIS_IP][0] = rp * cos(th + thp);
    na->x[i*3 + TIS_IP][1] = rp * sin(th + thp);
    na->x[i*3 + TIS_IP][2] = dh * i + dhp;

    na->x[i*3 + TIS_IS][0] = rs * cos(th + ths);
    na->x[i*3 + TIS_IS][1] = rs * sin(th + ths);
    na->x[i*3 + TIS_IS][2] = dh * i + dhs;

    ic = na->iseq[ i ];
    rb = rbarr[ ic ];
    thb = thbarr[ ic ];
    dhb = dhbarr[ ic ];
    na->x[i*3 + TIS_IB][0] = rb * cos(th + thb);
    na->x[i*3 + TIS_IB][1] = rb * sin(th + thb);
    na->x[i*3 + TIS_IB][2] = dh * i + dhb;
  }
}



/* remove the center of mass motion */
static void na_rmcom(double (*x)[D], int n)
{
  int i;
  double rc[D] = {0};

  for ( i = 0; i < n; i++ )
    vinc(rc, x[i]);
  vsmul(rc, 1./n);
  for ( i = 0; i < n; i++ )
    vdec(x[i], rc);
}



/* annihilate the total angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
static void na_shiftang(double (*x)[D], double (*v)[D], int n)
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



/* open an LJ system */
static na_t *na_open(const char *seq, double rc, int model)
{
  na_t *na;
  int i, d, n, nr;

  xnew(na, 1);
  na->nr = nr = strlen(seq);
  xnew(na->seq, nr);
  xnew(na->iseq, nr);
  for ( i = 0; i < nr; i++ ) {
    na->seq[i] = toupper( seq[i] );
    if ( na->seq[i] == 'A' ) {
      na->iseq[i] = 0;
    } else if ( na->seq[i] == 'C' ) {
      na->iseq[i] = 1;
    } else if ( na->seq[i] == 'G' ) {
      na->iseq[i] = 2;
    } else if ( na->seq[i] == 'U' || na->seq[i] == 'T' ) {
      na->iseq[i] = 3;
    } else {
      fprintf(stderr, "unknown residue type %c at %d\n", na->seq[i], i);
      free(na->seq);
      free(na->iseq);
      free(na);
      return NULL;
    }
  }

  na->apr = model;
  na->na = n = nr * na->apr;
  na->dof = n * D - D * (D+1)/2;
  na->rc = rc;
  na->rc2 = rc * rc;

  xnew(na->x, n);
  xnew(na->v, n);
  xnew(na->f, n);

  if ( model == 2 ) {
    na_initchain2(na);
  } else {
    na_initchain3(na);
  }

  /* initalize random velocities */
  for (i = 0; i < n; i++) {
    for ( d = 0; d < D; d++ ) {
      na->v[i][d] = randgaus();
    }
  }

  na_rmcom(na->v, n);
  na_shiftang(na->x, na->v, n);

  return na;
}



static void na_close(na_t *na)
{
  free(na->x);
  free(na->v);
  free(na->f);
  free(na);
}



/* write positions (and possibly velocities) */
__inline static int na_writepos(na_t *na,
    double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d, n = na->na;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %d\n", D, na->nr, na->na, (v != NULL));
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ )
      fprintf(fp, "%.14e ", x[i][d]);
    fprintf(fp, "%d ", i % 2);
    if ( v != NULL )
      for ( d = 0; d < D; d++ )
        fprintf(fp, "%.14e ", v[i][d]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



#define na_energy(na) \
  na->epot = na_energyTIS_low(na, na->x)

/* compute force, return energy */
__inline static double na_energyTIS_low(na_t *na, double (*x)[D])
{
  int i, j, ic, nr = na->nr, n = na->na;
  double ep = 0;

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS],
                   NULL, NULL);
    ep += ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   NULL, NULL);
    if ( i < nr - 1 ) {
      ep += ebondlen(R0_SP, K0_SP,
                     x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL);
    }
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondang(A0_PSB[ic], KA_PSB,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   NULL, NULL, NULL);
    if ( i < nr - 1 ) {
      ep += ebondang(A0_BSP[ic], KA_BSP,
                     x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL, NULL);
      ep += ebondang(A0_PSP, KA_PSP,
                     x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL, NULL);
    }
  }

  // excluded volume
  for ( i = 0; i < n - 1; i++ ) {
    for (j = i + 1; j < n; j++) {
      ep += ewca(WCA_SIG2, WCA_EPS, x[i], x[j], NULL, NULL);
    }
  }
  return ep;
}



#define na_force(na) \
  na->epot = na_forceTIS_low(na, na->x, na->f)

/* compute force, return energy */
__inline static double na_forceTIS_low(na_t *na, double (*x)[D], double (*f)[D])
{
  int i, j, ic, nr = na->nr, n = na->na;
  double ep = 0;

  for (i = 0; i < n; i++) vzero(f[i]);

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS],
                   f[i*3 + TIS_IP], f[i*3 + TIS_IS]);
    ep += ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondlen(R0_SP, K0_SP,
                     x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    }
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondang(A0_PSB[ic], KA_PSB,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondang(A0_BSP[ic], KA_BSP,
                     x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IB], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
      ep += ebondang(A0_PSP, KA_PSP,
                     x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    }
  }

  // excluded volume
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      ep += ewca(WCA_SIG2, WCA_EPS, x[i], x[j], f[i], f[j]);
    }
  }
  return ep;
}



/* velocity-verlet */
__inline static void na_vv(na_t *na, double dt)
{
  int i, n = na->na;
  double dth = dt * 0.5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(na->v[i], na->f[i], dth);
    vsinc(na->x[i], na->v[i], dt);
  }
  na_force(na);
  for (i = 0; i < n; i++) { /* VV part 2 */
    vsinc(na->v[i], na->f[i], dth);
  }
}



/* compute the kinetic energy */
static double na_ekin(double (*v)[D], int n)
{
  int i;
  double ek = 0;
  for ( i = 0; i < n; i++ ) {
    ek += vsqr( v[i] );
  }
  return ek/2;
}



#define na_vrescale(na, tp, dt) \
  na_vrescale_low(na->v, na->na, na->dof, tp, dt)

/* exact velocity rescaling thermostat */
__inline static double na_vrescale_low(double (*v)[D], int n,
    int dof, double tp, double dt)
{
  int i;
  double ek1, ek2, s, c, r, r2;

  tp *= BOLTZK;
  c = (dt < 700) ? exp(-dt) : 0;
  ek1 = na_ekin(v, n);
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp / 2);
  if (ek2 < 0) ek2 = 0;
  s = sqrt(ek2/ek1);
  for (i = 0; i < n; i++) vsmul(v[i], s);
  return ek2;
}



#if 0



/* displace a random particle i, return i */
static int na_randmv(na_t *na, double *xi, double amp)
{
  int i, d;

  i = (int) (rand01() * na->n);
  for ( d = 0; d < D; d++ )
    xi[d] = na->x[i][d] + (rand01() * 2 - 1) * amp;
  return i;
}



/* compute pair energy */
static int na_pair(double *xi, double *xj, double l, double invl,
    double rc2, double *u)
{
  double dx[D], dr2, invr2, invr6;

  dr2 = na_dist2(dx, xi, xj);
  if (dr2 < rc2) {
    invr2 = 1 / dr2;
    invr6 = invr2 * invr2 * invr2;
    *u  = 4.f * invr6 * (invr6 - 1);
    return 1;
  } else {
    *u = 0;
    return 0;
  }
}



/* return the energy change from displacing x[i] to xi */
__inline static double na_depot(na_t *na, int i, double *xi)
{
  int j, n = na->n;
  double l = na->l, invl = 1/l, rc2 = na->rc2, u, du;

  u = 0;
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) continue;
    if ( na_pair(na->x[i], na->x[j], rc2, &du) ) {
      u -= du;
    }
    if ( na_pair(xi, na->x[j], rc2, &du) ) {
      u += du;
    }
  }
  return u;
}



/* commit a particle displacement */
__inline static void na_commit(na_t *na, int i,
    const double *xi, double du)
{
  vcopy(na->x[i], xi);
  na->ep0 += du;
  na->epot += du;
}



/* Metropolis algorithm */
__inline static int na_metro(na_t *na, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0.;

  i = na_randmv(na, xi, amp);
  du = na_depot(na, i, xi);
  if ( du < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -bet * du ) );
  }
  if ( acc ) {
    na_commit(na, i, xi, du);
    return 1;
  }
  return 0;
}
#endif


#endif /* NACORE_H__ */

#include "util.h"
#include "vct.h"
#include "mat.h"

typedef struct {
  double dist_oh; // in angstroms
  double ang_hoh; // in degrees
  double mass_o;
  double mass_h;
} settle_water_param_t;

settle_water_param_t spc_water_param = {
  .dist_oh = 1.0,
  .ang_hoh = 109.28,
  .mass_o = 16.00000,
  .mass_h =  1.00800,
};

settle_water_param_t tip3p_water_param = {
  .dist_oh = 0.9572,
  .ang_hoh = 104.52,
  .mass_o = 16.00000,
  .mass_h =  1.00800,
};

typedef struct {
  double dist_oh, dist_hh, dist_on, dist_hn;
  double ang_hoh; // in radians
  double mass_o, mass_h;
  double mass_total;
  double frac_mass_h;
  double invmass[3];
  double vinvmat[3][3];
} settle_water_t;

#define SETTLE_NPAIR 3

int settle_pair_ids_[SETTLE_NPAIR][2] = {
  {1, 2},
  {2, 0},
  {0, 1},
};


void settle_print_vec_(double *v, const char *name)
{
  printf("%s: %+8.5f %+8.5f %8.5f\n", name, v[0], v[1], v[2]);
}

void settle_print_mat_(double (*m)[3], const char *name)
{
  int i;
  printf("%s:\n", name);
  for (i = 0; i < 3; i++) {
    printf("    %+8.5f %+8.5f %8.5f\n", m[i][0], m[i][1], m[i][2]);
  }
  printf("\n");
}



// compute invmass, invmat used for velocities
void settle_water_prepare_v(settle_water_t *sw)
{
  sw->invmass[0] = 1.0 / sw->mass_o;
  sw->invmass[1] = 1.0 / sw->mass_h;
  sw->invmass[2] = 1.0 / sw->mass_h;

  double x[3][3] = {{0}};
  x[1][0] = -sw->dist_on;
  x[1][1] = +sw->dist_hn;
  x[2][0] = -sw->dist_on;
  x[2][1] = -sw->dist_hn;

  int ipr, jpr;
  int i0, i1, j0;
  double dxpr[SETTLE_NPAIR][3];

  // compute the pairwise displacement vectors
  for (ipr = 0; ipr < SETTLE_NPAIR; ipr++) {
    i0 = settle_pair_ids_[ipr][0];
    i1 = settle_pair_ids_[ipr][1];
    vdiff(dxpr[ipr], x[i0], x[i1]);
  }

  double mat[SETTLE_NPAIR][SETTLE_NPAIR];
  for (ipr = 0; ipr < SETTLE_NPAIR; ipr++) {
    i0 = settle_pair_ids_[ipr][0];
    i1 = settle_pair_ids_[ipr][1];
    for (jpr = 0; jpr < SETTLE_NPAIR; jpr++) {
      double imass;
      if (ipr == jpr) {
        imass = -(sw->invmass[i0] + sw->invmass[i1]);
      } else {
        //int j1 = settle_pair_ids_[jpr][1];
        j0 = settle_pair_ids_[jpr][0];
        if (i1 == j0) {
          imass = sw->invmass[i1];
        } else { // i0 == j1
          imass = sw->invmass[i0];
        }
      }
      //printf("ipr %d = (%d, %d), jpr %d = (%d, %d), imass %g\n", ipr, i0, i1, jpr, j0, j1, imass);
      mat[ipr][jpr] = imass * vdot(dxpr[ipr], dxpr[jpr]);
    }
  }
  //settle_print_mat_(mat, "Velocity matrix");

  // invert the matrix
  minv(sw->vinvmat, mat);
  //settle_print_mat_(sw->vinvmat, "Inverse velocity matrix");
}


settle_water_t *settle_water_open(settle_water_param_t *p)
{
  settle_water_t *sw;

  xnew(sw, 1);

  sw->dist_oh = p->dist_oh;
  sw->ang_hoh = p->ang_hoh * M_PI/180;
  sw->dist_on = sw->dist_oh * cos(sw->ang_hoh/2);
  sw->dist_hn = sw->dist_oh * sin(sw->ang_hoh/2);
  sw->dist_hh = sw->dist_hn * 2;
  sw->mass_o = p->mass_o;
  sw->mass_h = p->mass_h;
  sw->mass_total = sw->mass_o + sw->mass_h * 2;
  sw->frac_mass_h = sw->mass_h / sw->mass_total;
  settle_water_prepare_v(sw);
  return sw;
}

void settle_water_close(settle_water_t *sw)
{
  free(sw);
}


void settle_water_apply_coordinates(settle_water_t *sw,
    double (*x0)[3], double (*x1)[3], double (*xout)[3])
{
  double xab0[3], xac0[3];
  double vz0[3], z0[3], dz0;
  double xab1[3], xac1[3], zab, zac;
  double vaxis[3], daxis2;
  double frac_mass_h = sw->frac_mass_h;

  // 1. apply a rotation around an axis in the x-y plane
  //    to make the z-coordinates right

  // compute the coordinates
  vdiff(xab0, x0[0], x0[1]);
  vdiff(xac0, x0[0], x0[2]);
  vcross(vz0, xab0, xac0);
  dz0 = vnorm(vz0);
  vsmul2(z0, vz0, 1.0/dz0);
  printf("vz0 %g %g %g\n", vz0[0], vz0[1], vz0[2]);

  vdiff(xab1, x1[0], x1[1]);
  zab = vdot(xab1, z0);
  vdiff(xac1, x1[0], x1[2]);
  zac = vdot(xac1, z0);

  vlincomb(vaxis, xab0, xac0, zac, -zab);
  printf("vaxis %g %g %g; zab %g, zac %g\n", vaxis[0], vaxis[1], vaxis[2], zab, zac);
  daxis2 = vsqr(vaxis);
  double xab4[3], xac4[3];
  double sinth, sinth2, one_minus_costh;
  if (daxis2 < 1e-12) { // z1 == z2 == 0
    // no first-rotation is needed
    vcopy(xab4, xab0);
    vcopy(xac4, xac0);
  } else {
    double daxis = sqrt(daxis2);
    //vsmul(vaxis, 1.0/daxis);
    sinth = daxis/dz0;
    sinth2 = sinth*sinth;
    one_minus_costh = sinth2/(1+sqrt(1-sinth2));
    printf("daxis %g, sinth %.8f, costh %.8f, dz0 %g\n", daxis, sinth, 1-one_minus_costh, dz0);
    // compute the rotated xab0
    double xab0para[3], xab0perp[3], dxab0[3];
    vsmul2(xab0para, vaxis, vdot(vaxis, xab0)/daxis2);
    vdiff(xab0perp, xab0, xab0para);
    printf("xab %g %g %g\n", xab0[0], xab0[1], xab0[2]);
    printf("xabpara %g %g %g\n", xab0para[0], xab0para[1], xab0para[2]);
    printf("xabperp %g %g %g\n", xab0perp[0], xab0perp[1], xab0perp[2]);

    // vz0*zab
    vlincomb(dxab0, xab0perp, z0, -one_minus_costh, zab);
    vadd(xab4, xab0, dxab0);
    printf("xab0 %g, %g, %g, dist %g\n", xab0[0], xab0[1], xab0[2], vnorm(xab0));
    printf("xab4 %g, %g, %g, dist %g\n", xab4[0], xab4[1], xab4[2], vnorm(xab4));

    // compute the rotated xac0
    double xac0para[3], xac0perp[3], dxac0[3];
    vsmul2(xac0para, vaxis, vdot(vaxis, xac0)/daxis2);
    vdiff(xac0perp, xac0, xac0para);
    vlincomb(dxac0, xac0perp, z0, -one_minus_costh, zac);
    vadd(xac4, xac0, dxac0);
    printf("xac0 %g, %g, %g, dist %g\n", xac0[0], xac0[1], xac0[2], vnorm(xac0));
    printf("xac4 %g, %g, %g, dist %g\n", xac4[0], xac4[1], xac4[2], vnorm(xac4));
  }

  // 2. apply the x, y rotation
  //
  // Try to find a rotation matrix such that
  //
  //   mB vB0 x (R vBA4 - vBA1)
  // + mC vC0 x (R vCA4 - vCA1) = 0
  //
  // vB0 x (R - 1) vBA4 + vC0 x (R - 1) vCA4 = vB0 x (vBA1 - vBA4) + vC0 x (vCA1 - vCA4)

  // compute the center of mass
  // xa = xcom + frac_mass_h * [(xa - xb) + (xa - xc)]
  double xbo0[3], xco0[3], xbcsum0[3], xcom0[3];
  vadd(xbcsum0, xab0, xac0);
  vsadd(xcom0, x0[0], xbcsum0, -frac_mass_h);
  vdiff(xbo0, x0[1], xcom0);
  vdiff(xco0, x0[2], xcom0);

  double dxab[3], dxac[3];
  double vx1[3], vx2[3], rhs;
  // NOTE: perhaps casting to the x-y plane would accelerate the computation
  vdiff(dxab, xab1, xab4);
  vcross(vx1, xbo0, dxab);
  vdiff(dxac, xac1, xac4);
  vcross(vx2, xco0, dxac);
  //printf("%g %g %g\n", vx1[0], vx1[1], vx1[2]);
  //printf("%g %g %g\n", vx2[0], vx2[1], vx2[2]);
  vinc(vx1, vx2);
  rhs = vdot(vx1, z0);

  double va1[3], va2[3], a, b;
  a = vdot(vinc(vcross(va1, xbo0, xab4),
                vcross(va2, xco0, xac4)), z0);
  b = vdot(xbo0, xab4) + vdot(xco0, xac4);
  // solving the equation of theta
  // a * [cos(theta) - 1] + b * sin(theta) = rhs
  //printf("a %g, b %g, rhs %g\n", a, b, rhs);
  // solve the quadratic equation for sin(theta)
  double gam = rhs + a, a2b2 = a*a+b*b;
  sinth = (b*gam - a*sqrt(a2b2-gam*gam))/a2b2;
  sinth2 = sinth*sinth;
  one_minus_costh = sinth2/(1+sqrt(1-sinth2));

  // apply the rotation to xab and xac
  double xab4para[3], xab4perp[3], dxab4[3], uab4[3], xab5[3];
  vsmul2(xab4para, z0, vdot(z0, xab4));
  vdiff(xab4perp, xab4, xab4para);
  vcross(uab4, z0, xab4perp);
  vlincomb(dxab4, xab4perp, uab4, -one_minus_costh, sinth);
  vadd(xab5, xab4, dxab4);
  //printf("xab5 %+8.5f %+8.5f %+8.5f; %g\n", xab5[0], xab5[1], xab5[2], vnorm(xab5));

  double xac4para[3], xac4perp[3], dxac4[3], uac4[3], xac5[3];
  vsmul2(xac4para, z0, vdot(z0, xac4));
  vdiff(xac4perp, xac4, xac4para);
  vcross(uac4, z0, xac4perp);
  vlincomb(dxac4, xac4perp, uac4, -one_minus_costh, sinth);
  vadd(xac5, xac4, dxac4);
  //printf("xac5 %+8.5f %+8.5f %+8.5f; %g\n", xac5[0], xac5[1], xac5[2], vnorm(xac5));

  // coordinates of the center of mass satisfies
  //  xa = xcom + frac_mass_h * (xa - xb) + frac_mass_h * (xa - xc)
  //
  // since the center of mass is unchanged from x1 to x5
  //
  // xa^5 = xa^1 + frac_mass_h*((xab^5 + xac^5)
  //                           -(xab^1 + xac^1))
  double x1sum[3], x5sum[3], dx0[3];
  vadd(x1sum, xab1, xac1);
  vadd(x5sum, xab5, xac5);
  vdiff(dx0, x5sum, x1sum);
  vsadd(xout[0], x1[0], dx0, frac_mass_h);
  vdiff(xout[1], xout[0], xab5);
  vdiff(xout[2], xout[0], xac5);
}


void settle_water_apply_velocities(settle_water_t *sw,
    double (*x)[3], double (*v0)[3], double (*v)[3])
{
  int i, j, ipr;
  double dv[3][3] = {{0}};
  double dxpr[SETTLE_NPAIR][3], dvpr[SETTLE_NPAIR][3], xvpr[SETTLE_NPAIR];

  for (ipr = 0; ipr < SETTLE_NPAIR; ipr++) {
    i = settle_pair_ids_[ipr][0];
    j = settle_pair_ids_[ipr][1];
    vdiff(dxpr[ipr], x[i], x[j]);
    vdiff(dvpr[ipr], v0[i], v0[j]);
    xvpr[ipr] = vdot(dxpr[ipr], dvpr[ipr]);
  }

  for (ipr = 0; ipr < SETTLE_NPAIR; ipr++) {
    double lambda = vdot(sw->vinvmat[ipr], xvpr);
    //printf("xvpr %g %g %g; lambda %g\n", xvpr[0], xvpr[1], xvpr[2], lambda);
    i = settle_pair_ids_[ipr][0];
    j = settle_pair_ids_[ipr][1];
    vsinc(dv[i], dxpr[ipr], +lambda);
    vsinc(dv[j], dxpr[ipr], -lambda);
    //printf("%g %g %g; lambda %g\n", dv[i][0], dv[i][1], dv[i][2], lambda);
  }

  for (i = 0; i < 3; i++) {
    vsadd(v[i], v0[i], dv[i], sw->invmass[i]);
  }
  //exit(1);
}



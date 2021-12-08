#include "settle.h"
#include "mtrand.h"


void settle_water_get_com(double (*x)[3], double xcom[3], double frac_mass_h)
{
  double dx[2][3];

  vcopy(xcom, x[0]);
  vdiff(dx[0], x[1], x[0]);
  vsinc(xcom, dx[0], frac_mass_h);
  vdiff(dx[1], x[2], x[0]);
  vsinc(xcom, dx[1], frac_mass_h);
}


void settle_water_com_shift(double (*x)[3], double xcom[3], double (*xo)[3], double frac_mass_h)
{
  settle_water_get_com(x, xcom, frac_mass_h);
  vdiff(xo[0], x[0], xcom); // TODO: probably unnecessary
  vdiff(xo[1], x[1], xcom);
  vdiff(xo[2], x[2], xcom);
}



void test_single(settle_water_t *sw, double x1[3][3], double v0[3][3], double v[3][3])
{
  double x0[3][3] = {
    {0, 0, 0},
    {-sw->dist_on,  sw->dist_hh/2, 0},
    {-sw->dist_on, -sw->dist_hh/2, 0},
  };
  double xout[3][3];

  printf("BEFORE: %g %g %g\n", vdist(x1[0], x1[1]), vdist(x1[0], x1[2]), vdist(x1[1], x1[2])); 
  settle_print_vec_(x1[0], "x0");
  settle_print_vec_(x1[1], "x1");
  settle_print_vec_(x1[2], "x2");

  settle_water_apply_coordinates(sw, x0, x1, xout);

  printf("AFTER: %g %g %g\n", vdist(xout[0], xout[1]), vdist(xout[0], xout[2]), vdist(xout[1], xout[2]));
  settle_print_vec_(xout[0], "x0");
  settle_print_vec_(xout[1], "x1");
  settle_print_vec_(xout[2], "x2");

  settle_water_apply_velocities(sw, xout, v0, v);

  double dxpr[3][3], dvpr[3][3];
  printf("VELOCITY: %g %g %g\n",
      vdot(vdiff(dxpr[2], xout[0], xout[1]), vdiff(dvpr[2], v[0], v[1])),
      vdot(vdiff(dxpr[0], xout[1], xout[2]), vdiff(dvpr[0], v[1], v[2])),
      vdot(vdiff(dxpr[1], xout[2], xout[0]), vdiff(dvpr[1], v[2], v[0])));
  settle_print_vec_(v[0], "v0");
  settle_print_vec_(v[1], "v1");
  settle_print_vec_(v[2], "v2");
}



void test_singles()
{
  settle_water_t *sw = settle_water_open(&spc_water_param);
  double x1[][3][3] = {
    // symmetric model, the second rotation is not necessary 
    {{0.02, -0.0, 0.},
     {-sw->dist_on+0.02,  sw->dist_hh/2+0.05, 0.03},
     {-sw->dist_on+0.02, -sw->dist_hh/2-0.05, 0.03}},

    // a 90 deg rotation in the x-y plane
    {{0.0, 0.0, 0.},
     { sw->dist_hh/2,  sw->dist_on, 0.0},
     {-sw->dist_hh/2,  sw->dist_on, 0.0}},

    // generic example
    {{0.0, 0.0, 0.},
    {-sw->dist_on+0.01,  sw->dist_hh/2+0.02, 0.03},
    {-sw->dist_on+0.07, -sw->dist_hh/2+0.05, 0.06}},
  };
  double v0[3][3] = {
    { 1.01,  1.02, -1.03},
    {-1.04, -1.05,  1.06},
    {-1.07,  1.08, -1.09}}, v[3][3];
  int i;

  for ( i = 0; i < 3; i++ ) {
    printf("Model %d:\n", i + 1);
    test_single(sw, x1[i], v0, v);
    printf("\n\n");
  }

  settle_water_close(sw);
}


void apply_rotation(double *x, double *dx, double *u, double theta)
{
  double dot, x1[3], xpara[3], xperp[3], xcross[3], y[3], z[3];

  vdiff(x1, x, dx);
  //vcopy(x1, x);
  dot = vdot(x1, u);
  vsmul2(xpara, u, dot);
  vlincomb(xperp, x1, u, 1.0, -dot);
  vcross(xcross, u, x1);
  vlincomb(y, xperp, xcross, cos(theta), sin(theta));
  vadd(x, xpara, y);
}


void test_randoms()
{
  settle_water_t *sw = settle_water_open(&spc_water_param);
  double x0[3][3] = {
    {0, 0, 0},
    {-sw->dist_on,  sw->dist_hh/2, 0},
    {-sw->dist_on, -sw->dist_hh/2, 0},
  };
  double mass[3] = {sw->mass_o, sw->mass_h, sw->mass_h};
  double x0com[3], xo0[3][3];
  settle_water_get_com(x0, x0com, sw->frac_mass_h);
  vdiff(xo0[0], x0[0], x0com);
  vdiff(xo0[1], x0[1], x0com);
  vdiff(xo0[2], x0[2], x0com);

  double x1[3][3] = {{0}},
         xout[3][3] = {{0}};
  int trial;

  for (trial = 0; trial < 1; trial++) {
    printf("Random model %d:\n", trial);
    memcpy(x1, x0, sizeof(double)*9);

    // apply constraint forces
    int i, j, k, ids[][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (k = 0; k < 3; k++) {
      double vij[3];
      i = ids[k][0];
      j = ids[k][1];
      vdiff(vij, x0[i], x0[j]);
      vsmul(vij, 0.012+k*0.015);
      vsinc(x1[i], vij,  1.0/mass[i]);
      vsinc(x1[j], vij, -1.0/mass[j]);
    }

    for (i = 0; i < 3; i++) {
      printf("  %d: %+.6f %+.6f %+.6f\n", i, x1[i][0], x1[i][1], x1[i][2]);
    }
    printf("  AB: %6.4f vs %6.4f\n", vdist(x1[0], x1[1]), sw->dist_oh);
    printf("  AC: %6.4f vs %6.4f\n", vdist(x1[0], x1[2]), sw->dist_oh);
    printf("  BC: %6.4f vs %6.4f\n", vdist(x1[1], x1[2]), sw->dist_hh);
    double x0cross[3][3], cross_total = 0;
    for (i = 0; i < 3; i++) {
      vcross(x0cross[i], xo0[i], x1[i]);
      cross_total += x0cross[i][2] * mass[i];
    }
    printf("x-y plane cross product total %g\n\n", cross_total);

    // apply a random rotation
    double u[3] = {0}, dx[3] = {0}, theta = 0.05*(rand01() - 0.5);
    randdir(u);
    randdir(dx);
    for (i = 0; i < 3; i++ ) {
      apply_rotation(x1[i], dx, u, theta);
    }

    settle_water_apply_coordinates(sw, x0, x1, xout);

    for (i = 0; i < 3; i++) {
      printf("  %d: %+.6f %+.6f %+.6f\n", i, xout[i][0], xout[i][1], xout[i][2]);
    }
    printf("  AB: %6.4f vs %6.4f\n", vdist(xout[0], xout[1]), sw->dist_oh);
    printf("  AC: %6.4f vs %6.4f\n", vdist(xout[0], xout[2]), sw->dist_oh);
    printf("  BC: %6.4f vs %6.4f\n", vdist(xout[1], xout[2]), sw->dist_hh);
    double xout_cross[3][3];
    for (cross_total = 0, i = 0; i < 3; i++) {
      vcross(xout_cross[i], xo0[i], xout[i]);
      cross_total += xout_cross[i][2] * mass[i];
    }
    printf("x-y plane cross product total %g\n", cross_total);

    printf("\n\n");
  }

  settle_water_close(sw);
}


int main(void)
{
  test_singles();
  //test_randoms();
  return 0;
}




#include "mtrand.h"
#include "pca.h"
#include "argopt.h"

int chainlen = D + 1;
//int chainlen = 3;
//int chainlen = 4;
//int chainlen = 5;
double tp = 300;
long nsteps = 1000000;
double mddt = 0.002; /* in ps */
enum { NONE, VRESCALE, LANGEVIN, NHCHAIN, ANDERSEN, THERMOSTAT_COUNT };
const char *thermostats[] = {"none", "v-rescale", "Langevin", "NHChain", "Andersen", "count"};
int thstat = LANGEVIN;
//int thstat = ANDERSEN;
//int thstat = VRESCALE;
double thdt = 0.1; /* time step for velocity rescaling or NH-chain */
double langdt = 0.1;
int nstlog = 10;
const char *fnlog = "polymer.log";
const char *fnpos = "polymer.xyz";
const char *fnhis = "rad.his";

//int transform = TRANSFORM_NONE;
//int transform = TRANSFORM_COM;
//int transform = TRANSFORM_HEAD;
int transform = TRANSFORM_RMSD;

#define NNHC 5
double nhc_zeta[NNHC];
double nhc_zmass[NNHC] = {1, 1, 1, 1, 1};

double mass_u = 12;
double mass; /* to be computed from mass_u */
double bond0 = 3.79;
double kbond = 222.5;
double theta0_deg = 113;
//double theta0_deg = 90;
//double theta0_deg = 180;
double theta0; /* to be computed from theta0_deg */
double ktheta = 58.35;
//double ktheta = 1000;
double phi0_deg = 180;
double phi0; /* to be computed from phi0_deg */
double kdih1 = 1000;
double kdih3 = 0.0;



typedef struct {
  int n;
  int dof;
  double *m; /* mass */
  double (*x)[D];
  double (*v)[D];
  double (*f)[D];
  double epot, ekin;
  double (*xt)[D]; /* transformed coordinates */
  double (*xref)[D]; /* reference coordinates */
} polymer_t;



static void polymer_close(polymer_t *p)
{
  free(p->m);
  free(p->x);
  free(p->v);
  free(p->f);
  free(p->xt);
  free(p->xref);
  free(p);
}


/* print out the angular momentum */
__inline static void polymer_printang(polymer_t *p, const char *flag)
{
  int i, n = p->n;
#if D == 2
  double ang;
  for ( i = 0; i < n; i++ ) {
    ang += p->m[i] * vcross(p->xref[i], p->xt[i]);
  }
  printf("%s: %g\n", flag, ang);
#elif D == 3
  double x[D], ang[D];
  vzero(ang);
  for ( i = 0; i < n; i++ ) {
    vcross(x, p->xref[i], p->xt[i]);
    vsinc(ang, x, p->m[i]);
  }
  printf("%s: %g %g %g | norm %g\n", flag, ang[0], ang[1], ang[2], vnorm(ang));
#endif
}

/* to use the output:
 * splot "polymer.xyz" w lp pt 7 ps 5 */
static int polymer_write(polymer_t *p, double (*x)[D], const char *fn)
{
  FILE *fp;
  int i;

  if ((fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  if ( x == NULL ) {
    /* by default, write the transformed coordinates */
    x_transform(p->x, p->xt, p->xref, p->m, p->n, transform);
    x = p->xt;
  }

  fprintf(fp, "# %d %d\n", p->n, D);
  for ( i = 0; i < p->n; i++ ) {
#if D == 2
    fprintf(fp, "%.18f %.18f\n", x[i][0], x[i][1]);
#elif D == 3
    fprintf(fp, "%.18f %.18f %.18f\n", x[i][0], x[i][1], x[i][2]);
#endif
    //if ( i > 0 && i < p->n - 1 ) printf("i %d: %g\n", i, vang(x[i-1], x[i], x[i+1],NULL,NULL,NULL)*180/M_PI);
  }
  fclose(fp);
  return 0;
}



static void polymer_logpos(polymer_t *p, FILE *fp,
    double (*x)[D], long t, int header)
{
  int i, j, n = p->n;

  x_transform(x, p->xt, p->xref, p->m, p->n, transform);

  if ( fp == NULL ) return;

  if ( header ) {
    /* report the masses */
    fprintf(fp, "# %d %d", n, n * D);
    for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < D; j++ ) {
        fprintf(fp, " %g", p->m[i]);
      }
    }
    fprintf(fp, "\n");
  }

  /* report the transformed positions */
  fprintf(fp, "%ld", t);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < D; j++ )
      fprintf(fp, " %.6f", p->xt[i][j]);
  }
  fprintf(fp, "\n");
}

/* compute the force and energy */
double polymer_force(polymer_t *p)
{
  int i, n = p->n;

  p->epot = 0;
  for ( i = 0; i < n; i++ ) {
    vzero(p->f[i]);
  }

  /* bonds */
  for ( i = 0; i < n - 1; i++ ) {
    p->epot += md_potbond(p->x[i], p->x[i + 1], bond0,
            kbond, p->f[i], p->f[i + 1]);
  }

  /* angles */
  for ( i = 0; i < n - 2; i++ ) {
    p->epot += md_potang(p->x[i], p->x[i + 1], p->x[i + 2], theta0,
          ktheta, p->f[i], p->f[i + 1], p->f[i + 2]);
  }

#if D == 3
  /* dihedrals */
  for ( i = 0; i < n - 3; i++ ) {
    p->epot += md_potdih13(p->x[i], p->x[i + 1], p->x[i + 2], p->x[i + 3], phi0,
      kdih1, kdih3, p->f[i], p->f[i + 1], p->f[i + 2], p->f[i + 3]);
  }
#endif

  return p->epot;
}

/* remove center of mass motion, linear and angular */
static void polymer_rmcom(polymer_t *p, double (*x)[D], double (*v)[D])
{
  rmcom(x, p->m, p->n);
  rmcom(v, p->m, p->n);
  shiftang(x, v, p->m, p->n);
}

polymer_t *polymer_open(int n, double tp0)
{
  polymer_t *p;
  double dx[2][D] = {{1}, {1}}, s, ang;
  int i, j;

  ang = (PI - theta0) / 2;
  dx[0][0] = cos(ang);
  dx[0][1] = sin(ang);
  dx[1][0] = dx[0][0];
  dx[1][1] = -dx[0][1];
  xnew(p, 1);
  p->n = n;
  xnew(p->m, n);
  xnew(p->x, n);
  xnew(p->v, n);
  xnew(p->f, n);
  xnew(p->xt, n);
  xnew(p->xref, n);

  for ( i = 0; i < n; i++ ) {
    p->m[i] = mass;
  }
  // non-uniform mass
  //p->m[0] /= 10;
  //p->m[n-1] /= 10;

  /* generate the initial configuration as a zigzag */
  vzero(p->xref[0]);
  for ( i = 1; i < n; i++ ) {
    vsadd(p->xref[i], p->xref[i-1], dx[i % 2], bond0);
  }
  rmcom(p->xref, p->m, p->n);

  for ( i = 0; i < n; i++ ) {
    vcopy(p->x[i], p->xref[i]);
    /* add some random noise to the positions */
    s = 0.1 * sqrt( KB * tp0 / kbond );
    for ( j = 0; j < D; j++ ) {
      p->x[i][j] += s * randgaus();
    }
  }
  rmcom(p->x, p->m, p->n);

  /* initialize velocities */
  for (i = 0; i < n; i++) {
    s = sqrt( KB * tp0 / p->m[i] );
    for ( j = 0; j < D; j++ ) {
      p->v[i][j] = s * randgaus();
    }
  }
  if ( thstat != LANGEVIN && thstat != ANDERSEN ) {
    polymer_rmcom(p, p->x, p->v); /* remove center of mass motion */
    p->dof = n * D - D * (D + 1) / 2;
  } else {
    p->dof = n * D;
  }
  polymer_force(p);
  p->ekin = md_ekin(p->v, p->m, n);
  return p;
}

/* velocity Verlet */
int polymer_vv(polymer_t *p, double dt)
{
  int i, n = p->n;
  double dth = 0.5 * dt;

  for ( i = 0; i < n; i++ ) { /* VV part 1 */
    vsinc(p->v[i], p->f[i], dth / p->m[i]);
    vsinc(p->x[i], p->v[i], dt);
  }

  p->epot = polymer_force(p); /* calculate force */

  for ( i = 0; i < n; i++ ) { /* VV part 2 */
    vsinc(p->v[i], p->f[i], dth / p->m[i]);
  }

  return 0;
}

void polymer_thermostat(polymer_t *p)
{
  static int t;

  ++t;
  if ( thstat == VRESCALE ) {
    p->ekin = md_vrescale(p->v, p->m, p->n, p->dof,
        KB * tp, 0.5 * thdt);
  } else if ( thstat == LANGEVIN ) {
    p->ekin = md_langevin(p->v, p->m, p->n,
        KB * tp, 0.5 * langdt);
  } else if ( thstat == NHCHAIN ) {
    p->ekin = md_nhchain(p->v, p->m, p->n, p->dof,
        KB * tp, 0.5 * thdt, NNHC, nhc_zeta, nhc_zmass);
  } else if ( thstat == ANDERSEN ) {
    if ( t % 10 == 0 ) {
      int i, j;
      double s;
      i = (int) (rand01() * p->n);
      s = sqrt( KB * tp / p->m[i] );
      for ( j = 0; j < D; j++ ) {
        p->v[i][j] = s * randgaus();
      }
    }
    p->ekin = md_ekin(p->v, p->m, p->n);
  } else {
    p->ekin = md_ekin(p->v, p->m, p->n);
  }
}

static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &chainlen, "number of particles");
  argopt_addx(ao, "--trans", "%list", &transform, "coordinate transformation for PCA", transforms, TRANSFORM_COUNT);
  argopt_add(ao, "-t", "%ld", &nsteps, "number of MD steps");
  argopt_add(ao, "-T", "%lf", &tp, "temperature in K");
  argopt_add(ao, "-m", "%lf", &mass_u, "mass in atomic mass unit (u)");
  argopt_add(ao, "--dt", "%lf", &mddt, "MD time step");
  argopt_addx(ao, "--th", "%list", &thstat, "thermostat type", thermostats, THERMOSTAT_COUNT);
  argopt_add(ao, "--thdt", "%lf", &thdt, "thermostat time step");
  argopt_add(ao, "--langdt", "%lf", &langdt, "Langevin damping rate");
  argopt_add(ao, "-g", "%d", &nstlog, "logging period");
  argopt_add(ao, "-b", "%lf", &bond0, "equilibrium bond length");
  argopt_add(ao, "--Kb", "%lf", &kbond, "spring constant of bonds");
  argopt_add(ao, "-a", "%lf", &theta0_deg, "equilibrium angle in degrees");
  argopt_add(ao, "--Ka", "%lf", &ktheta, "spring constant of angles");
  argopt_add(ao, "-d", "%lf", &phi0_deg, "equilibrium dihedral in degrees");
  argopt_add(ao, "--Kd", "%lf", &kdih1, "spring constant of angles for the cos(phi-phi0) term");
  argopt_add(ao, "--Kt", "%lf", &kdih3, "spring constant of angles for the cos(3*(phi-phi0)) term");
  argopt_add(ao, "--log", NULL, &fnlog, "log file name");
  argopt_add(ao, "--pos", NULL, &fnpos, "position file name");
  argopt_add(ao, "--his", NULL, &fnhis, "histogram file");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

static double polymer_getrad(polymer_t *p, double (*x)[D])
{
  double miner[D][D], detmi, smass = 0;
  int i;

  detmi = mominer(miner, x, p->m, p->n, 0);
  for ( i = 0; i < p->n; i++ )
    smass += p->m[i];
#if D == 2
  return sqrt(detmi / smass);
#elif D == 3
  return sqrt(1.5 * pow(detmi, 1./3) / smass);
#endif
}

#define HIST_N  10000
#define HIST_DX 0.002
#define HIST_XMAX (HIST_N * HIST_DX)
double hist[HIST_N+1];

static int hist_save(double *h, int n, double dx, const char *fn)
{
  int i;
  double s = 0, norm;
  FILE *fp;

  for ( i = 0; i < n; i++ ) s += h[i];
  norm = 1/(s*dx);
  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ ) {
    if ( h[i] > 0 ) {
      fprintf(fp, "%g %g %g\n", (i+.5)*dx, h[i]*norm, h[i]);
    }
  }
  fclose(fp);
  return 0;
}

int main(int argc, char **argv)
{
  polymer_t *p;
  pca_t *pca;
  FILE *fplog = NULL;
  long t;
  double sumep = 0, sumek = 0;

  doargs(argc, argv);
  mass = mass_u / 418.4;
  theta0 = theta0_deg * PI / 180;
  phi0 = phi0_deg * PI / 180;

  p = polymer_open(chainlen, tp);
  if ( strcmpfuzzy(fnlog, "null") != 0 ) {
    if ( (fplog = fopen(fnlog, "w")) == NULL ) {
      fprintf(stderr, "cannot open log file\n");
      return -1;
    }
  }
  polymer_logpos(p, fplog, p->xref, 0, 1);
  pca = pca_open(p->n * D);
  pca_setmass(pca, p->m);
  pca_setxref(pca, (double *) p->xref);
  for ( t = 1; t <= nsteps; t++ ) {
    polymer_thermostat(p);
    polymer_vv(p, mddt);
    if ( t % 10 == 0 ) {
      if ( thstat != LANGEVIN && thstat != ANDERSEN ) {
        /* need to regularly remove COM motion */
        polymer_rmcom(p, p->x, p->v);
      } else {
        rmcom(p->x, p->m, p->n);
      }
    }
    polymer_thermostat(p);

    sumep += p->epot;
    sumek += p->ekin;
    if ( t % nstlog == 0 ) {
      double rad;
      polymer_logpos(p, fplog, p->x, t, 0);
      pca_add(pca, (double *) p->xt, TRANSFORM_NONE);
      rad = polymer_getrad(p, p->xt);
      if ( rad < HIST_XMAX ) {
        hist[(int)(rad/HIST_DX)] += 1;
      }
    }

    /* print out the progress */
    if ( t % 100000 == 0 ) {
      printf("t %ld, %g%%      \r", t, 100.*t/nsteps);
      fflush(stdout);
    }
    //printf("%ld %g %g %g\n", t, p->epot, p->ekin, p->epot + p->ekin);
  }
  pca_analyze(pca, KB * tp, transform);
  polymer_write(p, NULL, fnpos);
  hist_save(hist, HIST_N, HIST_DX, fnhis);
  printf("T %g, kT %g, epot %g, ekin %g, etot %g\n", tp, KB * tp, sumep/nsteps, sumek/nsteps, (sumep+sumek)/nsteps);

  if ( fplog != NULL ) fclose(fplog);
  polymer_close(p);
  pca_close(pca);
  return 0;
}

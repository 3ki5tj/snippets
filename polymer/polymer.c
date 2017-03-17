#define D 3
#include "vct.h"
#include "mat.h"
#include "mdutil.h"
#include "mtrand.h"

#define xnew(x, n) if ( (x = calloc((n), sizeof(*(x)))) == NULL ) exit(1);

const double KB = 0.0019872041; /* kcal/mol/K */

int chainlen = D + 1;
//int chainlen = 4;
double tp = 300;
long nsteps = 1000000;
double mddt = 0.002; /* in ps */
enum { NONE, VRESCALE, LANGEVIN, NHCHAIN, ANDERSEN};
int thstat = ANDERSEN;
//int thstat = VRESCALE;
double thdt = 0.1; /* time step for velocity rescaling or NH-chain */
double langdt = 0.1;
int nstlog = 20;

enum { TRANSFORM_NONE, TRANSFORM_COM, TRANSFORM_HEAD, TRANSFORM_RMSD };
//int transform = TRANSFORM_COM;
//int transform = TRANSFORM_HEAD;
int transform = TRANSFORM_RMSD;

#define NNHC 5
double nhc_zeta[NNHC];
double nhc_zmass[NNHC] = {1, 1, 1, 1, 1};

double mass = 12 / 418.4;
double bond0 = 3.79;
double kbond = 222.5;
double theta0 = 113 * PI / 180;
double ktheta = 58.35;
double phi0 = PI;
double kdih1 = 100.;
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
static void polymer_printang(polymer_t *p, const char *flag)
{
  int i, n = p->n;
#if D == 2
  double ang;
  for ( i = 0; i < n; i++ ) {
    ang += p->m[i] * vcross(p->xref[i], p->xt[i]);
  }
  printf("%s: %g\n", flag, ang);
#else
  double x[D], ang[D];
  vzero(ang);
  for ( i = 0; i < n; i++ ) {
    vcross(x, p->xref[i], p->xt[i]);
    vsinc(ang, x, p->m[i]);
  }
  printf("%s: %g %g %g\n", flag, ang[0], ang[1], ang[2]);
#endif
}

/* remove the center of mass motion as well as overall rotations */
static void polymer_transform(polymer_t *p, int type)
{
  int i, n = p->n;
  double rot[D][D], trans[D], v[D], newx[D];

  if ( type == TRANSFORM_NONE ) {
    for ( i = 0; i < n; i++ )
      vcopy(p->xt[i], p->x[i]);
  } else if ( type == TRANSFORM_COM ) {
    for ( i = 0; i < n; i++ ) {
      vcopy(p->xt[i], p->x[i]);
    }
    /* remove the center of mass motion */
    rmcom(p->xt, p->m, n);
    /* remove the overall rotation
     * by treating xt as velocity */
    //polymer_printang(p, "before");
    shiftangv(p->xref, p->xt, p->m, n);
    //polymer_printang(p, "after"); getchar();
  } else if ( type == TRANSFORM_HEAD ) {
    /* 1. displace x[0] at the origin */
    vcopy(trans, p->x[0]);
    for ( i = 0; i < n; i++ ) {
      vdiff(p->xt[i], p->x[i], trans);
    }
#if D == 2
    vdiff(v, p->xt[1], p->xt[0]);
    vnormalize(v);
    for ( i = 0; i < n; i++ ) {
      newx[0] =  v[0] * p->xt[i][0] + v[1] * p->xt[i][1];
      newx[1] = -v[1] * p->xt[i][0] + v[0] * p->xt[i][1];
      vcopy(p->xt[i], newx);
    }
#else
    {
      double x[D] = {1, 0, 0}, y[D] = {0, 1, 0}, rot[D][D];
      /* 2. rotate p->xt[1] to the x axis */
      mrotvv(rot, p->xt[1], x);
      for ( i = 1; i < n; i++ ) {
        mmxv(newx, rot, p->xt[i]);
        vcopy(p->xt[i], newx);
      }
      /* 3. rotate p->xt[2] on to the x-y plane */
      vdiff(v, p->xt[2], p->xt[1]);
      v[0] = 0; /* remove the x component */
      mrotvv(rot, v, y);
      for ( i = 2; i < n; i++ ) {
        mmxv(newx, rot, p->xt[i]);
        vcopy(p->xt[i], newx);
      }
    }
#endif
  } else if ( type == TRANSFORM_RMSD ) {
    vrmsd(p->x, p->xt, p->xref, p->m, n, 0, rot, trans);
  }
}

/* to use the output:
 * splot "polymer.xyz" w lp pt 7 ps 5 */
static int polymer_write(polymer_t *p, const char *fn)
{
  FILE *fp;
  int i;

  if ((fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  polymer_transform(p, transform);

  fprintf(fp, "# %d %d\n", p->n, D);
  for ( i = 0; i < p->n; i++ ) {
#if D == 2
    fprintf(fp, "%g %g\n", p->xt[i][0], p->xt[i][1]);
#else
    fprintf(fp, "%g %g %g\n", p->xt[i][0], p->xt[i][1], p->xt[i][2]);
#endif
  }
  fclose(fp);
  return 0;
}



void polymer_logpos(polymer_t *p, FILE *fp, long t, int header)
{
  int i, j, n = p->n;

  if ( header ) {
    /* report the masses */
    fprintf(fp, "# %d %d", n, n * D);
    for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ) {
        fprintf(fp, " %g", p->m[i]);
      }
    }
    fprintf(fp, "\n");
  }

  polymer_transform(p, transform);

  /* report the transformed positions */
  fprintf(fp, "%ld", t);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < D; j++ )
      fprintf(fp, " %.6f", p->xt[i][j]);
  }
  fprintf(fp, " %g %g\n", p->epot, p->ekin);
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
  rmcom(v, p->m, p->n);
  shiftangr(x, v, p->m, p->n);
}

polymer_t *polymer_open(int n, double tp0)
{
  polymer_t *p;
  double dx[2][D] = {{1}, {1}}, s, ang;
  int i, j;

  ang = (PI - theta0)/2;
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

int main(void)
{
  polymer_t *p;
  FILE *fplog;
  long t;

  p = polymer_open(chainlen, tp);
  if ( (fplog = fopen("polymer.log", "w")) == NULL ) {
    fprintf(stderr, "cannot open log file\n");
    return -1;
  }
  polymer_logpos(p, fplog, 0, 1);
  for ( t = 1; t <= nsteps; t++ ) {
    polymer_thermostat(p);
    polymer_vv(p, mddt);
    if ( thstat != LANGEVIN && thstat != ANDERSEN
      && t % 10 == 0 ) {
      /* need to regularly remove COM motion */
      polymer_rmcom(p, p->x, p->v);
    }
    polymer_thermostat(p);
#if 0
    if ( t % 1000 == 0 ) {
      double sx = p->x[0][0]+p->x[1][0]+p->x[2][0]+p->x[3][0];
      double sv = p->v[0][0]+p->v[1][0]+p->v[2][0]+p->v[3][0];
      double sf = p->f[0][0]+p->f[1][0]+p->f[2][0]+p->f[3][0];
      printf("%ld %g %g %g\n", t, sx, sv, sf);
      if (fabs(sx) > 1.1 || fabs(sv) > 0.1 || fabs(sf) > 0.1 ) {
        int i;
        polymer_write(p, "fail.xyz");
        for ( i = 0; i < p->n; i++ ) {
          printf("%d: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i,
              p->x[i][0], p->x[i][1], p->x[i][2],
              p->v[i][0], p->v[i][1], p->v[i][2],
              p->f[i][0], p->f[i][1], p->f[i][2]);
        }
        exit(1);
      }
    }
#endif

    if ( t % nstlog == 0 ) {
      polymer_logpos(p, fplog, t, 0);
    }
    //printf("%ld %g %g %g\n", t, p->epot, p->ekin, p->epot + p->ekin);
  }
  polymer_write(p, "polymer.xyz");
  fclose(fplog);
  polymer_close(p);
  return 0;
}

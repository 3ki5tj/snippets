#ifndef PCA_H__
#define PCA_H__



#include "mdutil.h"
#include "eig.h"

#define xnew(x, n) if ( (x = calloc((n), sizeof(*(x)))) == NULL ) exit(1);

#define KB        0.0019872041 /* kcal/mol/K */
#define AVOGADRO  6.022140857e23
#define HBAR      (1.054571800e-34*AVOGADRO*1e12/4184) /* kcal/mol*ps */

enum {
  TRANSFORM_NONE,
  TRANSFORM_COM, /* remove the center of mass motion and the overall rotation */
  TRANSFORM_HEAD, /* place the first atom at the origin,
                     the second on the x-axis, the third on the x-y plane */
  TRANSFORM_RMSD, /* best align with a reference structure,
                     equivalent to TRANSFORM_COM */
  TRANSFORM_COUNT
};

const char *transforms[] = {"none", "COM", "head", "RMSD", "count"};

typedef struct {
  int n;
  int np;
  int nin;
  long cnt; /* number of frames */
  double *mass; /* mass, array of length n/D */
  double *m; /* mass, array of length n */
  double *sqrtm;
  double *x;
  double *xt;
  double *xref;
  double *sx;
  double *sxx;
  double *ave;
  double *cov;
  double *sig;
  double *eval;
  double *evec;
  double *xin;
  double *xrefin;
  double *sxin;
  double *sxxin;
  double *avein;
  double *covin;
  double *sigin;
  double *evalin;
  double *evecin;
} pca_t;

static pca_t *pca_open(int n)
{
  pca_t *pca;
  int i, j, nin;

  xnew(pca, 1);
  pca->n = n;
  pca->np = n / D;
  pca->nin = nin = n - D*(D+1)/2;
  xnew(pca->mass, n);
  xnew(pca->m, n);
  xnew(pca->sqrtm, n);
  xnew(pca->x, n);
  xnew(pca->xt, n);
  xnew(pca->xref, n);
  xnew(pca->sx, n);
  xnew(pca->sxx, n * n);
  xnew(pca->ave, n);
  xnew(pca->cov, n * n);
  xnew(pca->sig, n);
  xnew(pca->eval, n);
  xnew(pca->evec, n * n);
  for ( i = 0; i < n; i++ ) { /* default mass */
    pca->m[i] = 12 / 418.4;
    if ( i % D == 0 ) pca->mass[i/D] = pca->m[i];
    pca->sqrtm[i] = sqrt( pca->m[i] );
  }
  for ( i = 0; i < n; i++ ) {
    pca->sx[i] = 0;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->sxx[i*n+j] = 0;
    }
  }
  xnew(pca->xin, nin);
  xnew(pca->xrefin, nin);
  xnew(pca->sxin, nin);
  xnew(pca->sxxin, nin * nin);
  xnew(pca->avein, nin);
  xnew(pca->covin, nin * nin);
  xnew(pca->sigin, nin);
  xnew(pca->evalin, nin);
  xnew(pca->evecin, nin * nin);
  for ( i = 0; i < nin; i++ ) {
    pca->sxin[i] = 0;
  }
  for ( i = 0; i < nin; i++ ) {
    for ( j = 0; j < nin; j++ ) {
      pca->sxxin[i*nin+j] = 0;
    }
  }
  pca->cnt = 0;
  return pca;
}

static void pca_close(pca_t *pca)
{
  free(pca->mass);
  free(pca->m);
  free(pca->sqrtm);
  free(pca->x);
  free(pca->xt);
  free(pca->xref);
  free(pca->sx);
  free(pca->sxx);
  free(pca->ave);
  free(pca->cov);
  free(pca->sig);
  free(pca->eval);
  free(pca->evec);
  free(pca->xin);
  free(pca->xrefin);
  free(pca->sxin);
  free(pca->sxxin);
  free(pca->avein);
  free(pca->covin);
  free(pca->sigin);
  free(pca->evalin);
  free(pca->evecin);
  free(pca);
}

/* remove the center of mass motion as well as overall rotations */
static void x_transform(double (*x)[D], double (*xt)[D],
    double (*xref)[D], const double *m, int n, int type)
{
  int i;
  double rot[D][D], trans[D], v[D], newx[D];

  if ( type == TRANSFORM_NONE ) {
    for ( i = 0; i < n; i++ )
      vcopy(xt[i], x[i]);
  } else if ( type == TRANSFORM_COM ) {
    /* 1. copy the coordinates */
    for ( i = 0; i < n; i++ )
      vcopy(xt[i], x[i]);
    /* 2. remove the center of mass motion */
    rmcom(xt, m, n);
    /* 3. remove the overall rotation by treating xt as velocity */
    //polymer_printang(p, "before");
    shiftangv(xref, xt, m, n);
    //polymer_printang(p, "after"); getchar();
  } else if ( type == TRANSFORM_HEAD ) {
    /* 1. displace x[0] at the origin */
    vcopy(trans, x[0]);
    for ( i = 0; i < n; i++ ) {
      vdiff(xt[i], x[i], trans);
    }
#if D == 2
    vdiff(v, xt[1], xt[0]);
    vnormalize(v);
    for ( i = 0; i < n; i++ ) {
      newx[0] =  v[0] * xt[i][0] + v[1] * xt[i][1];
      newx[1] = -v[1] * xt[i][0] + v[0] * xt[i][1];
      vcopy(xt[i], newx);
    }
#else
    {
      double x[D] = {1, 0, 0}, y[D] = {0, 1, 0}, rot[D][D];
      /* 2. rotate xt[1] to the x axis */
      mrotvv(rot, xt[1], x);
      for ( i = 1; i < n; i++ ) {
        mmxv(newx, rot, xt[i]);
        vcopy(xt[i], newx);
      }
      /* 3. rotate xt[2] on to the x-y plane */
      vdiff(v, xt[2], xt[1]);
      v[0] = 0; /* remove the x component */
      mrotvv(rot, v, y);
      for ( i = 2; i < n; i++ ) {
        mmxv(newx, rot, xt[i]);
        vcopy(xt[i], newx);
      }
    }
#endif
  } else if ( type == TRANSFORM_RMSD ) {
    vrmsd(x, xt, xref, m, n, 0, rot, trans);
  }
}

/* set the mass */
__inline static void pca_setmass(pca_t *pca, const double *m)
{
  int i, j, id;

  for ( i = 0; i < pca->np; i++ ) {
    pca->mass[i] = m[i];
    for ( j = 0; j < D; j++ ) {
      id = i * D + j;
      pca->m[id] = m[i];
      pca->sqrtm[id] = sqrt( pca->m[id] );
    }
  }
}

/* add a multiple of 2*PI to `ang` such that it is closest to `ang0` */
static double wrapang(double ang, double ang0)
{
  ang -= ang0;
  while ( ang >  PI ) ang -= PI * 2;
  while ( ang < -PI ) ang += PI * 2;
  /* now the difference between ang and ang0 should lie in [-PI, PI] */
  return ang + ang0;
}

/* compute the internal coordinates */
static void pca_getxin(pca_t *pca, double (*x)[D],
    double *xin, double *xrefin)
{
  int i, id = 0, np = pca->np;
  double dx[D], ang, ang0;

  /* bonds */
  for ( i = 0; i < np - 1; i++ ) {
    xin[id++] = vnorm( vdiff(dx, x[i], x[i+1]) );
  }

  /* angles */
  for ( i = 0; i < np - 2; i++ ) {
    ang = vang(x[i], x[i+1], x[i+2], NULL, NULL, NULL);
    ang0 = (xrefin != NULL) ? xrefin[id] : 0;
    xin[id++] = wrapang(ang, ang0);
  }

#if D == 3
  /* dihedrals */
  for ( i = 0; i < np - 3; i++ ) {
    ang = vdih(x[i], x[i+1], x[i+2], x[i+3], NULL, NULL, NULL, NULL);
    ang0 = (xrefin != NULL) ? xrefin[id] : 0;
    xin[id++] = wrapang(ang, ang0);
  }
#endif
}


/* set the reference structure */
static void pca_setxref(pca_t *pca, const double *x)
{
  int i;

  for ( i = 0; i < pca->n; i++ )
    pca->xref[i] = x[i];
  pca_getxin(pca, (vct *) pca->xref, pca->xrefin, NULL);
}

/* add a frame */
static void pca_add(pca_t *pca, const double *x,
    int transform)
{
  int i, j, n = pca->n, nin = pca->nin;

  /* transform coordinates */
  x_transform((vct *) x, (vct *) pca->xt, (vct *) pca->xref,
      pca->mass, pca->np, transform);

  for ( i = 0; i < n; i++ ) {
    pca->sx[i] += pca->xt[i];
  }

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->sxx[i*n+j] += pca->xt[i] * pca->xt[j];
    }
  }

  /* compute the internal coordinates */
  pca_getxin(pca, (vct *) pca->xt, pca->xin, pca->xrefin);

  for ( i = 0; i < nin; i++ ) {
    pca->sxin[i] += pca->xin[i];
  }

  for ( i = 0; i < nin; i++ ) {
    for ( j = 0; j < nin; j++ ) {
      pca->sxxin[i*nin+j] += pca->xin[i] * pca->xin[j];
    }
  }

  pca->cnt++;
}

__inline static pca_t *pca_load(const char *fn, long skip, int transform)
{
  pca_t *pca;
  FILE *fp;
  char buf[1024], *p;
  int i, dim, np, n, next, lines;
  long time;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return NULL;
  }
  /* determine the dimension */
  fgets(buf, sizeof buf, fp);
  sscanf(buf, "#%d%d%n", &np, &n, &next);
  dim = n / np;
  if ( dim != D ) {
    fprintf(stderr, "Dimension mismatch %d (file) vs %d (program)\n", dim, D);
    exit(1);
  }
  pca = pca_open(n);
  /* determine the masses */
  for ( p = buf + next, i = 0; i < n; i++ ) {
    sscanf(p, "%lf%n", &pca->m[i], &next);
    p += next;
  }
  for ( i = 0; i < n; i++ ) {
    pca->sqrtm[i] = sqrt( pca->m[i] );
  }

  lines = 0;
  while ( fgets(buf, sizeof buf, fp) ) {
    /* read in the coordinates */
    sscanf(buf, "%ld%n", &time, &next);
    //next = 0;
    for ( p = buf + next, i = 0; i < n; i++ ) {
      sscanf(p, "%lf%n", &pca->x[i], &next);
      p += next;
    }

    /* use the first frame as the reference */
    if ( lines == 0 ) pca_setxref(pca, pca->x);

    if ( ++lines <= skip ) continue;

    pca_add(pca, pca->x, transform);
  }
  fprintf(stderr, "Loaded %ld frames, skipped %ld\n",
      pca->cnt, skip);

  fclose(fp);
  return pca;
}

/* compute the mass-weighted covariance matrix */
static void pca_getcov(pca_t *pca)
{
  int i, j, n = pca->n;

  for ( i = 0; i < n; i++ ) {
    pca->ave[i] = pca->sx[i] / pca->cnt;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] = pca->sxx[i*n+j] / pca->cnt
        - pca->ave[i] * pca->ave[j];
    }
  }

  for ( i = 0; i < n; i++ ) {
    pca->ave[i] *= pca->sqrtm[i];
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] *= pca->sqrtm[i] * pca->sqrtm[j];
    }
  }

  for ( i = 0; i < n; i++ ) {
    pca->sig[i] = pca->cov[i*n+i] > 0 ? sqrt( pca->cov[i*n+i] ) : 1;
  }
}

/* compute the eigenvalues from the covariance matrix */
static double pca_coveig(pca_t *pca, double kT, int verbose)
{
  double del;
  int i, j, cnt, n = pca->n;

  if ( verbose ) {
    printf("Covariance matrix:\n");
    for ( i = 0; i < n; i++ ) {
      double ss = pca->cov[i*n+i] > 0 ? pca->sig[i] : 0;
      printf("%4d %11.7f: ", i, ss);
      for ( j = 0; j < n; j++ ) {
        printf(" %8.3f", pca->cov[i*n + j]/pca->sig[i]/pca->sig[j]);
      }
      printf("\n");
    }
  }

  /* determining the shift to improve stability */
  for ( del = 0, cnt = 0, i = 0; i < n; i++ ) {
    if ( pca->cov[i*n+i] > 0 ) {
      del += pca->cov[i*n+i];
      cnt += 1;
    }
  }
  del /= cnt;

  /* shift the diagonal */
  for ( i = 0; i < n; i++ ) {
    pca->cov[i*n+i] += del;
  }

  /* compute the eigenvalues */
  eigsym(pca->cov, pca->eval, pca->evec, n);

  /* unshift the eigenvalues */
  for ( i = 0; i < n; i++ ) {
    pca->eval[i] -= del;
    if ( pca->eval[i] < 0 ) {
      pca->eval[i] = 0;
    }
  }

  printf("1/omega in fs\n");

  for ( i = 0; i < n; i++ ) {
    double x = sqrt(pca->eval[i]/kT);
    x *= 1000;
    printf("%4d: %10.4f: ", i + 1, x);
    for ( j = 0; j < n; j++ ) {
      printf(" %7.3f", pca->evec[j*n + i]);
    }
    printf("\n");
  }

  return del;
}

/* compute the moment of inertia */
static double pca_miner(pca_t *pca, double miner[D][D], double (*x)[D])
{
  int i, j, k;
  double sxy[D][D], tr;

  /* compute the correction factor for the overall translation and rotation */
  /* compute the moment of inertia */
  mzero(sxy);
  for ( i = 0; i < pca->np; i++ ) {
    for ( j = 0; j < D; j++ ) {
      for ( k = 0; k < D; k++ ) {
        sxy[j][k] += pca->mass[i] * x[i][j] * x[i][k];
      }
    }
  }
  for ( tr = 0, j = 0; j < D; j++ ) {
    tr += sxy[j][j];
  }
  for ( j = 0; j < D; j++ ) {
    for ( k = 0; k < D; k++ ) {
      miner[j][k] = -sxy[j][k];
    }
    miner[j][j] += tr;
  }

  printf("Moment of inertia:\n");
  for ( j = 0; j < D; j++ ) {
    for ( k = 0; k < D; k++ ) {
      printf(" %10.6f", miner[j][k]);
    }
    printf("\n");
  }
#if D == 2
  return miner[0][0] + miner[1][1];
#elif D == 3
  return mdet(miner);
#endif
}

/* `vol` in 2D means the area of the simulation box */
static int pca_entxyz(pca_t *pca, double kT, double vol, int transform)
{
  int i, n = pca->n, nmodes = 0;
  double alpha, entc = 0, entq = 0;
  double smass, miner[D][D], det, enttr = 0;

  /* get the mass weighted covariance matrix */
  pca_getcov(pca);

  /* compute the eigenvalues */
  pca_coveig(pca, kT, 1);

  /* we compute the correction factor from the reference structure
   * but ideally it should be computed from the average of each frame */
  for ( smass = 0, i = 0; i < pca->np; i++ ) {
    smass += pca->mass[i];
  }
  det = pca_miner(pca, miner, (vct *) pca->xref);
#if D == 2
  enttr = log(2*PI*vol);
#elif D == 3
  enttr = log(8*PI*PI*vol);
#endif
  if ( transform == TRANSFORM_HEAD ) {
    double (*x)[D] = (vct *) pca->xref, b;
    enttr += 0.5*D*log(pca->mass[0]);
#if D >= 2
    b = vdist(x[0], x[1]);
    enttr += 0.5*(D-1)*log(pca->mass[1]*b*b);
#endif
#if D == 3
    {
      double ang;
      b = vdist(x[1], x[2]);
      ang = vang(x[0], x[1], x[2], NULL, NULL, NULL);
      enttr += 0.5*log(pca->mass[2]) + log(b*sin(ang));
    }
#endif
  } else {
    enttr += 0.5*(D*log(smass) + log(det));
  }
  /* The first 0.5 term in the following formula is from the kinetic
   * energetic of the overall linear and angular momenta */
  enttr += D*(D+1)/2*(0.5 + 0.5*log(kT/(2*PI)) - log(HBAR));

  /* compute the entropy */
  nmodes = 0;
  printf("hbar*omega/kT: ");
  for ( i = 0; i < n - D*(D+1)/2; i++ ) {
    alpha = pca->eval[i]/kT;
    //if ( alpha <= 0 ) continue;
    alpha = HBAR/sqrt(alpha)/kT;
    entc += 1 - log(alpha);
    entq += alpha/(exp(alpha)-1) - log(1-exp(-alpha));
    nmodes += 1;
    printf(" %g", alpha);
  }
  printf("\nEntropy from Cartesian coordinates PCA (%d modes, volume %g, detI %g):\n", nmodes, vol, det);
  printf("Classical: %12.7f kcal/mol/K = %10.6f kB | corrected %12.7f kcal/mol/K = %10.6f kB\n",
      entc * KB, entc, (entc + enttr) * KB, entc + enttr);
  printf("Quantum:   %12.7f kcal/mol/K = %10.6f kB | corrected %12.7f kcal/mol/K = %10.6f kB\n",
      entq * KB, entq, (entq + enttr) * KB, entq + enttr);

  return 0;
}

/* compute the covariance matrix */
static void pca_getcovin(pca_t *pca)
{
  int i, j, nin = pca->nin;
  double x;

  for ( i = 0; i < nin; i++ ) {
    pca->avein[i] = pca->sxin[i] / pca->cnt;
  }
  for ( i = 0; i < nin; i++ ) {
    for ( j = 0; j < nin; j++ ) {
      pca->covin[i*nin+j] = pca->sxxin[i*nin+j] / pca->cnt
        - pca->avein[i] * pca->avein[j];
    }
  }

  for ( i = 0; i < nin; i++ ) {
    x = pca->covin[i*nin+i];
    pca->sigin[i] = x > 0 ? sqrt(x) : 1;
  }
}

/* compute the eigenvalues from the covariance matrix */
static void pca_covineig(pca_t *pca, int verbose)
{
  int i, j, nin = pca->nin;

  /* normalize the covariance matrix by the standard deviations */
  for ( i = 0; i < nin; i++ ) {
    for ( j = 0; j < nin; j++ ) {
      pca->covin[i*nin + j] /= pca->sigin[i] * pca->sigin[j];
    }
  }

  if ( verbose ) {
    printf("Covariance matrix:\n");
    for ( i = 0; i < nin; i++ ) {
      printf("%4d %11.7f: ", i, pca->sigin[i]);
      for ( j = 0; j < nin; j++ ) {
        printf(" %8.3f", pca->covin[i*nin+j]);
      }
      printf("\n");
    }
  }

  /* compute the eigenvalues */
  eigsym(pca->covin, pca->evalin, pca->evecin, nin);

  for ( i = 0; i < nin; i++ ) {
    if ( pca->evalin[i] < 0 ) {
      pca->evalin[i] = 0;
    }
  }

  if ( verbose ) {
    printf("Eigenvalues x 1000: eigenvectors\n");
    for ( i = 0; i < nin; i++ ) {
      printf("%4d %10g: ", i + 1, sqrt(pca->evalin[i]) * 1000);
      for ( j = 0; j < nin; j++ ) {
        printf(" %7.3f", pca->evecin[j*nin + i]);
      }
      printf("\n");
    }
  }
}


/* entropy from internal coordinates */
static int pca_entint(pca_t *pca, double kT, double vol)
{
  int i, nin = pca->nin, nmodes = 0;
  double detcov, entc = 0;
  double mprod, enttr = 0;

  /* get the mass weighted covariance matrix */
  pca_getcovin(pca);

  /* compute the eigenvalues */
  pca_covineig(pca, 0);

  /* we compute the correction factor from the reference structure
   * but ideally it should be computed from the average of each frame */
  for ( mprod = 1, i = 0; i < pca->np; i++ ) {
    mprod *= pca->mass[i];
  }
  /* note: pca->n == pca->np * D */
#if D == 2
  enttr = log(2*PI*vol);
#else
  enttr = log(8*PI*PI*vol);
#endif
  enttr += 0.5*D*log(mprod);
  /* The first 0.5 term in the following formula is from the kinetic energetic
   * contribution from all momenta */
  enttr += pca->n * (0.5 + 0.5 * log(kT/(2*PI)) - log(HBAR));

  {
    double (*x)[D] = (vct *) pca->xref;
#if D == 2
    for ( i = 0; i < pca->np - 1; i++ ) {
      enttr += 0.5 * log( vdist2(x[i], x[i+1]) );
    }
#else
    for ( i = 0; i < pca->np - 1; i++ ) {
      enttr += log( vdist2(x[i], x[i+1]) );
    }
    for ( i = 0; i < pca->np - 2; i++ ) {
      double ang = vang(x[i], x[i+1], x[i+2], NULL, NULL, NULL);
      enttr += log( sin(ang) );
    }
#endif
  }

  /* compute the entropy */
  entc = 0;
  nmodes = 0;
  detcov = 1;
  for ( i = 0; i < nin; i++ ) {
    double ev = pca->evalin[i];
    detcov *= pca->sigin[i] * sqrt(ev);
    if ( ev <= 0 ) continue;
    entc += 0.5 * (1 + log(2*PI));
    printf(" %g/%g", pca->sigin[i], sqrt(ev));
    nmodes += 1;
  }
  entc += log(detcov);
  printf("\nEntropy from internal coordinates PCA (%d modes, volume %g, det(oov) = %g):\n", nmodes, vol, detcov);
  printf("Classical: %12.7f kcal/mol/K = %10.6f kB | corrected %12.7f kcal/mol/K = %10.6f kB\n",
      entc * KB, entc, (entc + enttr) * KB, entc + enttr);

  return 0;
}

static int pca_analyze(pca_t *pca, double kT, int transform)
{
  pca_entxyz(pca, kT, 1.0, transform);
  pca_entint(pca, kT, 1.0);
  return 0;
}


#endif /* PCA_H__ */

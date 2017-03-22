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

typedef struct {
  int n;
  int dim;
  long cnt; /* number of frames */
  double *m;
  double *sqrtm;
  double *x;
  double *xt;
  double *xref;
  double *ave;
  double *cov;
  double *eval;
  double *evec;
} pca_t;

static pca_t *pca_open(int n)
{
  pca_t *pca;
  int i, j;

  xnew(pca, 1);
  pca->n = n;
  xnew(pca->m, n);
  xnew(pca->sqrtm, n);
  xnew(pca->x, n);
  xnew(pca->xt, n);
  xnew(pca->xref, n);
  xnew(pca->ave, n);
  xnew(pca->cov, n * n);
  xnew(pca->eval, n);
  xnew(pca->evec, n * n);
  for ( i = 0; i < n; i++ ) {
    pca->m[i] = 12 / 418.4;
    pca->sqrtm[i] = sqrt( pca->m[i] );
  }
  for ( i = 0; i < n; i++ ) {
    pca->ave[i] = 0;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] = 0;
    }
  }
  pca->cnt = 0;
  return pca;
}

static void pca_close(pca_t *pca)
{
  free(pca->m);
  free(pca->sqrtm);
  free(pca->x);
  free(pca->xt);
  free(pca->xref);
  free(pca->ave);
  free(pca->cov);
  free(pca->eval);
  free(pca->evec);
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

/* set the reference structure */
static void pca_setxref(pca_t *pca, const double *x)
{
  int i;

  for ( i = 0; i < pca->n; i++ )
    pca->xref[i] = x[i];
}

/* add a frame */
static void pca_add(pca_t *pca, const double *x, int transform)
{
  int i, j, n = pca->n;

  /* transform coordinates */
  x_transform((vct *) x, (vct *) pca->xt, (vct *) pca->xref,
      pca->m, n / D, transform);

  for ( i = 0; i < n; i++ ) {
    pca->xt[i] *= pca->sqrtm[i];
    pca->ave[i] += pca->xt[i];
  }

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] += pca->xt[i] * pca->xt[j];
    }
  }

  pca->cnt++;
}

static void pca_getave(pca_t *pca)
{
  int i, j, n = pca->n;

  for ( i = 0; i < n; i++ ) {
    pca->ave[i] /= pca->cnt;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] = pca->cov[i*n+j]/pca->cnt
        - pca->ave[i] * pca->ave[j];
    }
  }
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


static int pca_analyze(pca_t *pca, double kT)
{
  int i, j, cnt = 0, n = pca->n, nmodes = 0;
  double del = 0, alpha, entc = 0, entq = 0;
  double entm = 0;
  double *sig;

  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ )
      pca->cov[i*n+j] /= kT;

  xnew(sig, n);
  for ( i = 0; i < n; i++ ) {
    sig[i] = pca->cov[i*n+i] > 0 ? sqrt( pca->cov[i*n+i] ) : 1;
  }

  for ( i = 0; i < n; i++ ) {
    double ss = pca->cov[i*n+i] > 0 ? 1000*sig[i] : 0;
    printf("%4d %8.3f: ", i, ss);
    for ( j = 0; j < n; j++ ) {
      printf(" %8.3f", pca->cov[i*n + j]/sig[i]/sig[j]);
    }
    printf("\n");
  }

  printf("n %d\n", n);
  /* determining the shift to improve stability*/
  for ( i = 0; i < n; i++ ) {
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

  eigsym(pca->cov, pca->eval, pca->evec, n);
#if 0
  /* recover */
  for ( i = 0; i < n; i++ ) {
    pca->cov[i*n+i] -= del;
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] *= kT;
    }
  }
#endif
  /* unshift the eigenvalues */
  for ( i = 0; i < n; i++ ) {
    pca->eval[i] -= del;
    if ( pca->eval[i] < 0 ) {
      pca->eval[i] = 0;
    }
  }

  /* compute the entropy */
  nmodes = 0;
  printf("hbar*omega/kT: ");
  for ( i = 0; i < n - D*(D+1)/2; i++ ) {
    alpha = pca->eval[i];
    //if ( alpha <= 0 ) continue;
    alpha = HBAR/sqrt(alpha)/kT;
    entc += 1 - log(alpha);
    entq += alpha/(exp(alpha)-1) - log(1-exp(-alpha));
    nmodes += 1;
    printf(" %g", alpha);
  }
  printf("\nEntropy: %g kcal/mol/K, %g kB (classical) %g kcal/mol/K, %g kB (quantum); %d none-zero modes\n",
      entc * KB, entc, entq * KB, entq, nmodes);

  printf("1/omega in fs\n");
  for ( i = 0; i < n; i++ ) {
    printf("ev %4d: %10.4f: ", i + 1, sqrt(pca->eval[i])*1000);
    for ( j = 0; j < n; j++ ) {
      printf(" %7.3f", pca->evec[j*n + i]);
    }
    printf("\n");
  }
  free(sig);

  return 0;
}


#endif /* PCA_H__ */

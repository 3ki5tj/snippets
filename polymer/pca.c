#include <stdio.h>
#include "eig.h"

#define xnew(x, n) if ( (x = calloc((n), sizeof(*(x)))) == NULL ) exit(1);

char *fnin = "polymer.log";

typedef struct {
  int n;
  int cnt;
  double *mass;
  double *sqrtm;
  double *x;
  double *ave;
  double *cov;
  double *eval;
  double *evec;
} pca_t;


pca_t *pca_load(const char *fn, int skip)
{
  pca_t *pca;
  FILE *fp;
  char buf[1024], *p;
  int i, j, np, n, next, lines;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return NULL;
  }
  /* determine the dimension */
  fgets(buf, sizeof buf, fp);
  sscanf(buf, "#%d%d%n", &np, &n, &next);
  xnew(pca, 1);
  pca->n = n;
  xnew(pca->mass, n);
  xnew(pca->sqrtm, n);
  xnew(pca->x, n);
  xnew(pca->ave, n);
  xnew(pca->cov, n * n);
  xnew(pca->eval, n);
  xnew(pca->evec, n);
  /* determine the masses */
  for ( p = buf + next, i = 0; i < np; i++ ) {
    sscanf(p, "%lf%n", &pca->mass[i*3], &next);
    pca->mass[i*3 + 1] = pca->mass[i*3];
    pca->mass[i*3 + 2] = pca->mass[i*3];
    p += next;
  }
  for ( i = 0; i < n; i++ ) {
    pca->sqrtm[i] = sqrt(pca->mass[i]);
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
  lines = 0;
  while ( fgets(buf, sizeof buf, fp) ) {
    if ( ++lines <= skip ) continue;
    /* read in the coordinates */
    for ( p = buf + next, i = 0; i < n; i++ ) {
      sscanf(p, "%lf%n", &pca->x[i], &next);
      p += next;
    }
    for ( i = 0; i < n; i++ ) {
      pca->x[i] *= pca->sqrtm[i];
      pca->ave[i] += pca->x[i];
    }
    for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ) {
        pca->cov[i*n+j] += pca->x[i] * pca->x[j];
      }
    }
    pca->cnt += 1;
  }
  fprintf(stderr, "done loading %d frames, skipped %d\n",
      pca->cnt, skip);

  for ( i = 0; i < n; i++ ) {
    pca->ave[i] /= pca->cnt;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] = pca->cov[i*n + j] / pca->cnt
        - pca->ave[i] * pca->ave[j];
    }
  }

  fclose(fp);
  return pca;
}


int pca_analyze(pca_t *pca, double kT)
{
  int i, j, n = pca->n;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] /= kT;
      printf(" %8.3f", pca->cov[i*n + j]);
    }
    printf("\n");
  }

  printf("n %d\n", n);

  eigsym(pca->cov, pca->eval, pca->evec, n);
  
  for ( i = 0; i < n; i++ ) {
    printf("ev %4d: %8.4f: ", i + 1, pca->eval[i]);
    for ( j = 0; j < n; j++ ) {
      printf(" %7.3f", pca->evec[j*n + i]);
    }
    printf("\n");
  }
  return 0;
}

int main(void)
{
  pca_t *pca;
  pca = pca_load(fnin, 100000);
  pca_analyze(pca, 0.6);
  return 0;
}

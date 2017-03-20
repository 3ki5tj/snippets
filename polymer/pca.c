#include <stdio.h>
#include "eig.h"

#define xnew(x, n) if ( (x = calloc((n), sizeof(*(x)))) == NULL ) exit(1);

const double KB = 0.0019872041; /* kcal/mol/K */
#define AVOGADRO 6.022140857e23
const double hbar = 1.054571800e-34*AVOGADRO*1e12/4184; /* kcal/mol*ps */

char *fnin = "polymer.log";
double tp = 300;
long nstskip = 10;

typedef struct {
  int n;
  int dim;
  int cnt;
  double *mass;
  double *sqrtm;
  double *x;
  double *ave;
  double *cov;
  double *eval;
  double *evec;
} pca_t;

void pca_close(pca_t *pca)
{
  free(pca->mass);
  free(pca->sqrtm);
  free(pca->x);
  free(pca->ave);
  free(pca->cov);
  free(pca->eval);
  free(pca->evec);
  free(pca);
}

pca_t *pca_load(const char *fn, long skip)
{
  pca_t *pca;
  FILE *fp;
  char buf[1024], *p;
  int i, j, np, n, next, lines;
  long time;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return NULL;
  }
  /* determine the dimension */
  fgets(buf, sizeof buf, fp);
  sscanf(buf, "#%d%d%n", &np, &n, &next);
  xnew(pca, 1);
  pca->n = n;
  pca->dim = n / np;
  xnew(pca->mass, n);
  xnew(pca->sqrtm, n);
  xnew(pca->x, n);
  xnew(pca->ave, n);
  xnew(pca->cov, n * n);
  xnew(pca->eval, n);
  xnew(pca->evec, n * n);
  /* determine the masses */
  for ( p = buf + next, i = 0; i < n; i++ ) {
    sscanf(p, "%lf%n", &pca->mass[i], &next);
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
    sscanf(buf, "%ld%n", &time, &next);
    //next = 0;
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
  fprintf(stderr, "Loaded %d frames, skipped %d\n",
      pca->cnt, skip);

  for ( i = 0; i < n; i++ ) {
    pca->ave[i] /= pca->cnt;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] = pca->cov[i*n + j] / pca->cnt
        - pca->ave[i] * pca->ave[j];
      //printf("%d %d %g %g %g\n", i, j, pca->cov[i*n+j], pca->ave[i], pca->ave[j]);
    }
  }

  fclose(fp);
  return pca;
}


int pca_analyze(pca_t *pca, double kT)
{
  int i, j, cnt = 0, n = pca->n, nmodes = 0;
  double del = 0, alpha, entc = 0, entq = 0;
  double *sig;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      pca->cov[i*n+j] /= kT;
    }
  }

  xnew(sig, n);
  for ( i = 0; i < n; i++ )
    sig[i] = pca->cov[i*n+i] > 0 ? sqrt( pca->cov[i*n+i] ) : 1;

  for ( i = 0; i < n; i++ ) {
    printf("%4d %8.3f: ", i, sig[i]*1000);
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
  for ( i = 0; i < n - pca->dim*(pca->dim+1)/2; i++ ) {
    alpha = pca->eval[i];
    //if ( alpha <= 0 ) continue;
    alpha = hbar/sqrt(alpha)/kT;
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

int main(int argc, char **argv)
{
  pca_t *pca;
  if ( argc > 1 ) fnin = argv[1];
  pca = pca_load(fnin, nstskip);
  pca_analyze(pca, KB * tp);
  pca_close(pca);
  return 0;
}

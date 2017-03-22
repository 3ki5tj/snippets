#include <stdio.h>
#include "pca.h"

char *fnin = "polymer.log";
double tp = 300;
long nstskip = 10;
int transform = TRANSFORM_RMSD;

int main(int argc, char **argv)
{
  pca_t *pca;
  if ( argc > 1 ) fnin = argv[1];
  pca = pca_load(fnin, nstskip, transform);
  pca_analyze(pca, KB * tp);
  pca_close(pca);
  return 0;
}

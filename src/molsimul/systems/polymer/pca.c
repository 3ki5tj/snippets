#include <stdio.h>
#include "pca.h"
#include "argopt.h"

char *fnin = "polymer.log";
double tp = 300;
long nstskip = 10;
int transform = TRANSFORM_RMSD;

static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  ao = argopt_open(0);
  argopt_add(ao, NULL, NULL, &fnin, "input file");
  argopt_addx(ao, "--trans", "%list", &transform, "coordinate transformation for PCA", transforms, TRANSFORM_COUNT);
  argopt_add(ao, "--skip", "%d", &nstskip, "number of steps to skip");
  argopt_add(ao, "-T", "%lf", &tp, "temperature in K");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

int main(int argc, char **argv)
{
  pca_t *pca;

  doargs(argc, argv);
  pca = pca_load(fnin, nstskip, transform);
  if ( pca == NULL ) return -1;
  pca_analyze(pca, KB * tp, transform);
  pca_close(pca);
  return 0;
}

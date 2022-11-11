#include "mat.h"

#define N 4

static int load(double (*x)[D], const char *fn)
{
  FILE *fp;
  char s[1024];
  int i, n = 0;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    return 0;
  }

  if ( fgets(s, sizeof s, fp) == NULL ) {
    return 0;
  }

  sscanf(s, "# %d", &n);

  for ( i = 0; i < n; i++ ) {
    if (fgets(s, sizeof s, fp) == NULL) {
      break;
    }
    sscanf(s, "%lf%lf%lf", &x[i][0], &x[i][1], &x[i][2]);
  }

  fclose(fp);

  return n;
}


int main(void)
{
  double xt[N][D], xref[N][D], x0[N][D], m[N], rmsd;
  int i, n;

  n = load(xref, "xref.xyz");
  load(x0, "x1.xyz");
  for ( i = 0; i < n; i++ ) {
    m[i] = 12/418.4;
  }

  rmsd = vrmsd(x0, xt, xref, m, n, 0, NULL, NULL);

  fprintf(stderr, "RMSD after alignment is %lf\n", rmsd);

  return 0;
}

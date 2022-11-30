/* molecular dynamics of hard-sphere fluid */

#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "hsmd.h"

int n = 108;
int d = 3;
real l = 6;
real ntime = 100;


static void domd(void)
{
  hsmd_t *m;
  real *box;
  int k;

  xnew(box, d);
  for ( k = 0; k < d; k++ ) box[k] = l;
  m = hsmd_open(n, d, box);
  while ( m->time < ntime ) {
    hsmd_step(m);
  }
  hsmd_close(m);
  free(box);
}



int main(void)
{
  domd();
  return 0;
}


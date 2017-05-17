#include "is2exact.h"


int main(int argc, char **argv)
{
  int n = 4, m = 4;
  double tp = 2.3;
  double *p;

  if ( argc == 2 ) {
    n = m = atoi(argv[1]);
  } else if ( argc >= 3 ) {
    n = atoi(argv[1]);
    m = atoi(argv[2]);
    if ( argc >= 4 ) tp = atof(argv[3]);
  }
  {
    double lnz1, eav1, cv1, lnz2 = 0, eav2 = 0, cv2 = 0;
    lnz1 = is2exact(n, m, 1./tp, &eav1, &cv1);
    lnz2 = is2_exact(n, m, 1./tp, &eav2, &cv2);
    printf("lnz %g,%g, eav %g,%g, cv %g,%g\n", lnz1, lnz2, eav1, eav2, cv1, cv2);
  }
  p = is2dos(n, m);
  is2dos_save(p, n, m);
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include "cisi.h"

int main(int argc, char **argv)
{
  double x = 3, ci, si;

  if ( argc > 1 ) x = atof(argv[1]);
  cisi(x, &ci, &si);
  printf("x %g, ci %.14f, si %.14f\n", x, ci, si);
  return 0;
}

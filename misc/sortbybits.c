/* sort 1 .. 2^n by the number of nonempty bits */
#include <stdio.h>
#include <stdlib.h>



#define NMAX 32



/* enumerate numbers from 1 to 2^n with m bits */
static void enumerate(int n, int m, int *st)
{
  int top, v, k;
  unsigned x = 0; /* the number represented by the stack */

  st[top = 0] = 0;
  while (top >= 0) {
    if ( st[top] < n ) { /* push */
      v = st[top++];
      x ^= 1u << v;
      /* push if we still have numbers to fill */
      if (top < m) {
        /* the number on the next level must be greater than
         * the number on this level */
        st[top] = v + 1;
        continue;
      } else { /* we are at the top, print and fall through to pop */
        printf("%#8x == %-8d : ", x, x);
        for (k = 0; k < m; k++) printf("%d ", st[k]);
        printf("\n");
      }
    }
    /* pop */
    v = st[--top];
    st[top] = v + 1; /* v is no longer available on top */
    x ^= 1u << v;
  }
}



int main(int argc, char **argv)
{
  int st[NMAX + 1];
  int n = 4, m;

  if (argc >= 2) n = atoi(argv[1]);
  for (m = 1; m <= n; m++)
    enumerate(n, m, st);
  return 0;
}


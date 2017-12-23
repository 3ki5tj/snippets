#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *mygetline(char **s, size_t *n, FILE *fp)
{
  const int bufsz = 10;
  char *p;
  size_t cnt, sz;

  if ( *s == NULL && *n == 0 ) {
    *n = bufsz;
    if ( (*s = calloc(*n, sizeof(char))) == NULL ) exit(-1);
  }
  p = *s;
  sz = *n;
  while ( 1 ) {
    if ( fgets(p, sz, fp) == NULL ) return NULL;
    cnt = strlen(*s);
    if ( (*s)[cnt-1] == '\n' ) {
      break;
    } else { /* line too long, expand the buffer */
      *n += bufsz;
      if ( (*s = realloc(*s, (*n)*sizeof(char))) == NULL ) exit(-1);
      p = *s + cnt;
      sz = bufsz;
    }
  }
  return *s;
}

int main(void)
{
  char *s = NULL;
  size_t n = 0;
  FILE *fp;

  fp = fopen("a.dat", "w");
  fprintf(fp, "01234567890123456789 a\n01234567890123456789 b\n");
  fclose(fp);
  fp = fopen("a.dat", "r");
  while (mygetline(&s, &n, fp)) {
    puts(s);
  }
  fclose(fp);
  return 0;
}

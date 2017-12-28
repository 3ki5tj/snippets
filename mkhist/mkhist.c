#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, DATA, VARIABLE };

/* supported operators:
 * + and - (binary or unary), *, /, % (remainder), ^ or ** (power),
 * <, <=, >, >=, ==, !=, && (logical AND), || (logical OR) and parentheses () */
const char *ops = ",+-*/%^_<>!=&|()";

#define VARNAME_MAX 128

typedef struct {
  int type;
  char s[VARNAME_MAX]; /* operator/variable/function name */
  double val;
  int col; /* column of the data */
} token_t;


/* get a token from the string and return a pointer
 * after the token */
static char *gettoken(token_t *t, char *s)
{
  char *p;
  int i;

  /* skip leading spaces */
  while ( *s && isspace(*s) ) s++;
  if ( *s == '\0' ) return NULL;

  t->s[0] = '\0';
  t->val = 0;
  t->col = 0;

  if ( strchr(ops, *s) != NULL ) {
    t->type = OPERATOR;
    t->s[0] = *s;
    t->s[1] = '\0';
    p = s + 1;
    if ( *s == '*' && s[1] == '*' ) {
      t->s[0] = '^'; /* convert "**" to "^" */
      p++;
    } else if ( (strchr("<>!=", *s) != NULL && s[1] == '=')
       || ( *s == '<' && s[1] == '>' )
       || ( *s == '&' && s[1] == '&' )
       || ( *s == '|' && s[1] == '|' ) ) {
      t->s[1] = s[1];
      t->s[2] = '\0';
      p++;
    }
  } else if ( isdigit(*s) || *s == '.' ) { /* .5 is understood as a number */
    t->type = NUMBER;
    t->val = strtod(s, &p);
  } else if ( isalpha(*s) ) {
    for ( p = s, i = 0; isalnum(*p) && i < VARNAME_MAX - 1; p++, i++ ) {
      t->s[i] = *p;
    }
    t->s[i] = '\0';
    while ( *p && isspace(*p) ) p++;
    if ( *p == '(' ) {
      t->type = FUNCTION;
      p++;
    } else {
      t->type = VARIABLE;
    }
  } else if ( *s == '$' ) {
    t->type = DATA;
    t->col = (int) strtol(s + 1, &p, 10);
  } else {
    fprintf(stderr, "Error: unknown token [%s]\n", s);
    exit(1);
  }
  return p;
}

/* string representation of token */
__inline static char *token2str(char *s, const token_t *tok)
{
  if ( tok->type == NUMBER ) {
    sprintf(s, "%g", tok->val);
  } else if ( tok->type == OPERATOR || tok->type == FUNCTION || tok->type == VARIABLE ) {
    sprintf(s, "%s", tok->s);
  } else if ( tok->type == DATA ) {
    sprintf(s, "$%d", tok->col);
  }
  return s;
}

static void copytoken(token_t *a, token_t *b)
{
  a->type = b->type;
  strncpy(a->s, b->s, VARNAME_MAX);
  a->val = b->val;
  a->col = b->col;
}

/* return the precedence of the operator */
static int preced(const char *s)
{
  static char prtable[256];
  static struct { const char *s; int p; } prlist[] = {
    { "<=", 9 }, { ">=", 9 },
    { "!=", 8 }, { "==", 8 }, { "<>", 8 },
    { "&&", 5 },
    { "||", 4 },
    { NULL, 0 } }, *prptr;

  if ( prtable['+'] == 0 ) { /* initialize the precedence table */
    prtable['('] = prtable[')'] = 16;
    prtable['_'] = prtable['!'] = 15;
    prtable['^'] = 13;
    prtable['*'] = prtable['/'] = prtable['%'] = 12;
    prtable['+'] = prtable['-'] = 11;
    prtable['<'] = prtable['>'] = 9;
    prtable[','] = 1;
  }

  if isalpha(*s) { /* functions take the highest precedence */
    return 1000;
  } else if ( s[1] == '\0' ) {
    return prtable[(int)(*s)];
  } else { /* search the list */
    for ( prptr = prlist; prptr->s != NULL; prptr++ ) {
      if ( strcmp(prptr->s, s) == 0 ) {
        return prptr->p;
      }
    }
  }
  fprintf(stderr, "unknown operator [%s]\n", s);
  return 0;
}

static int isleftassoc(const char *s)
{
  if ( *s == '^' ) return 0;
  return 1;
}

/* convert an expression to a stack of postfix expression
 * The shunting-yard algorithm:
 * https://en.wikipedia.org/wiki/Shunting-yard_algorithm
 * */
static token_t *parse2postfix(char *s)
{
  token_t *que; /* output queue */
  token_t *ost; /* operator stack */
  token_t *pos, *top, tok[2];
  char *p;
  int n;

  n = strlen(s);

  /* initialize the output queue */
  if ( (que = calloc(n + 1, sizeof(*que))) == NULL ) exit(-1);
  pos = que;

  /* initialize the operator stack */
  if ( (ost = calloc(n + 1, sizeof(*ost))) == NULL ) exit(-1);
  top = ost;
  strcpy(ost[0].s, ""); /* special operator */

  tok[1].type = NULLTYPE; /* for the previous token */

  /* when there are tokens to be read, read a token */
  for ( p = s; (p = gettoken(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == DATA || tok->type == VARIABLE ) {
      /* if the token is a number, then push it to the output queue */
      copytoken(pos++, tok);
    } else if ( tok->s[0] == '(' || tok->type == FUNCTION ) {
      /* push "(" or a function onto the operator stack */
      copytoken(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a "("
       * or a function */
      while ( top->s[0] != '(' && top->type != FUNCTION ) {
        /* pop operators from the operator stack onto the output queue */
        copytoken(pos++, top--);
        if ( top <= ost ) break;
      }
      if ( top->s[0] == '(' ) {
        top--; /* pop the left bracket from the stack */
      } else {
        copytoken(pos++, top--);
      }
    } else if ( tok->type == OPERATOR ) {
      if ( (tok->s[0] == '+' || tok->s[0] == '-') &&
           ( tok[1].type == NULLTYPE || tok[1].type == FUNCTION
          || (tok[1].type == OPERATOR && tok[1].s[0] != ')') ) ) {
        /* handle a unary operator */
        if ( tok->s[0] == '-' ) {
          tok->s[0] = '_'; /* change it to negation */
          copytoken(++top, tok); /* top = token */
        } /* skip the unary '+' */
      } else {
        while ( ( (preced(top->s) > preced(tok->s))
              || (preced(top->s) == preced(tok->s) && isleftassoc(top->s)) )
                 && top->s[0] != '(' && top->type != FUNCTION ) {
          /* pop operators from the operator stack onto the output queue */
          copytoken(pos++, top--); /* pos = top */
          if ( top <= ost ) break;
        }
        /* push the read operator onto the operator stack */
        copytoken(++top, tok); /* top = token */
      }
    }
    copytoken(tok+1, tok); /* make a copy */
  }

  /* while there are still operator tokens on the stack,
   * pop the operator onto the output queue */
  while ( top > ost )
    copytoken(pos++, top--);

  pos->type = NULLTYPE;
  free(ost);
  return que;
}

static double max(double a, double b) { return (a > b) ? a : b; }
static double min(double a, double b) { return (a < b) ? a : b; }
static double iif(double c, double a, double b) { return ( c != 0 ) ? a : b; }

typedef struct {
  char s[VARNAME_MAX];
  double (*f)();
  int narg;
} funcmap_t;

typedef struct {
  char s[VARNAME_MAX];
  double val;
} varmap_t;

static funcmap_t funcmap[] = {
  {"fabs", fabs, 1},
  {"abs", fabs, 1},
  {"sqrt", sqrt, 1},
  {"exp", exp, 1},
  {"sinh", sinh, 1},
  {"cosh", cosh, 1},
  {"tanh", tanh, 1},
  {"log", log, 1},
  {"ln", log, 1},
  {"log10", log10, 1},
  {"pow", pow, 2},
  {"sin", sin, 1},
  {"cos", cos, 1},
  {"tan", tan, 1},
  {"asin", asin, 1},
  {"acos", acos, 1},
  {"atan", atan, 1},
  {"atan2", atan2, 2},
  {"ceil", ceil, 1},
  {"floor", floor, 1},
  {"fmod", fmod, 2},
  {"min", min, 2},
  {"max", max, 2},
  {"if", iif, 3},
  {"iif", iif, 3},
  {"", NULL, 0}
};

#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif

#ifndef M_E
#define M_E    2.71828182845904523536028747135266250
#endif

static varmap_t varmap[] = {
  {"pi", M_PI},
  {"Pi", M_PI},
  {"PI", M_PI},
  {"e",  M_E},
  {"E", M_E},
  {"", 0}
};

/* evaluate the postfix expression */
static double evalpostfix(const token_t *que, const double *arr, int narr)
{
  int i, n, top;
  double *st, ans;
  const token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ( (st = calloc(n, sizeof(*st))) == NULL ) exit(-1);
  top = 0;

  for ( pos = que; pos->type != NULLTYPE; pos++ ) {
    if ( pos->type == NUMBER ) {
      st[top++] = pos->val;
    } else if ( pos->type == DATA ) {
      st[top++] = ( pos->col <= narr ) ? arr[pos->col - 1] : 0.0; /* array index off-by-one */
    } else if ( pos->type == VARIABLE ) {
      for ( i = 0; varmap[i].s[0] != '\0'; i++ ) {
        if ( strcmp(varmap[i].s, pos->s) == 0 ) {
          st[top++] = varmap[i].val;
          break;
        }
      }
    } else if ( pos->type == OPERATOR ) {
      if ( pos->s[0] == '_' ) { /* unary operators */
        st[top-1] = -st[top-1];
      } else if ( pos->s[0] == '!' && pos->s[1] == '\0' ) {
        st[top-1] = !( st[top-1] != 0 );
      } else if ( pos->s[0] == ',' ) { /* do nothing for comma */
        ;
      } else { /* binary operators */
        --top;
        if ( pos->s[0] == '+' ) {
          st[top-1] += st[top];
        } else if ( pos->s[0] == '-' ) {
          st[top-1] -= st[top];
        } else if ( pos->s[0] == '*' ) {
          st[top-1] *= st[top];
        } else if ( pos->s[0] == '/' ) {
          st[top-1] /= st[top];
        } else if ( pos->s[0] == '%' ) {
          st[top-1] = fmod(st[top-1], st[top]);
        } else if ( pos->s[0] == '^' ) {
          st[top-1] = pow(st[top-1], st[top]);
        } else if ( strcmp(pos->s, "<") == 0 ) {
          st[top-1] = ( st[top-1] < st[top] );
        } else if ( strcmp(pos->s, ">") == 0 ) {
          st[top-1] = ( st[top-1] > st[top] );
        } else if ( strcmp(pos->s, "<=") == 0 ) {
          st[top-1] = ( st[top-1] <= st[top] );
        } else if ( strcmp(pos->s, ">=") == 0 ) {
          st[top-1] = ( st[top-1] >= st[top] );
        } else if ( strcmp(pos->s, "!=") == 0 || strcmp(pos->s, "<>") == 0 ) {
          st[top-1] = ( st[top-1] != st[top] );
        } else if ( strcmp(pos->s, "==") == 0 ) {
          st[top-1] = ( st[top-1] == st[top] );
        } else if ( strcmp(pos->s, "&&") == 0 ) {
          st[top-1] = ( (st[top-1] != 0) && (st[top] != 0) );
        } else if ( strcmp(pos->s, "||") == 0 ) {
          st[top-1] = ( (st[top-1] != 0) || (st[top] != 0) );
        } else {
          fprintf(stderr, "Error: uknown operator [%s]\n", pos->s);
          exit(1);
        }
      }
    } else if ( pos->type == FUNCTION ) {
      for ( i = 0; funcmap[i].f != NULL; i++ ) {
        if ( strcmp(funcmap[i].s, pos->s) == 0 ) {
          if ( funcmap[i].narg == 1 ) {
            st[top-1] = (*funcmap[i].f)(st[top-1]);
          } else if ( funcmap[i].narg == 2 ) {
            --top;
            st[top-1] = (*funcmap[i].f)(st[top-1], st[top]);
          } else if ( funcmap[i].narg == 3 ) {
            top -= 2;
            st[top-1] = (*funcmap[i].f)(st[top-1], st[top], st[top+1]);
          }
          break;
        }
      }
    }
  }
  ans = st[0];
  free(st);
  return ans;
}

/* read a long line from file */
static char *mygetline(char **s, size_t *n, FILE *fp)
{
  const int bufsz = 1024;
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

/* get the maximal number of required data columns */
static int getmaxcol(token_t *p)
{
  int maxcol = 0;

  for ( ; p->type != NULLTYPE; p++ )
    if ( p->type == DATA && p->col > maxcol )
      maxcol = p->col;
  return maxcol;
}

/* display a 1D histogram on the terminal */
static void showhist(const double *h, double xmin, double dx, int xn)
{
  int i0 = 0, i, j;
  double m = 0;
  i0 = 0;
  for ( ; i0 < xn; i0++ ) if ( h[i0] > 0 ) break;
  for ( ; xn > i0; xn-- ) if ( h[xn-1] > 0 ) break;
  for ( i = i0; i < xn; i++ ) if ( h[i] > m ) m = h[i];
  for ( i = i0; i < xn; i++ ) {
    fprintf(stderr, "%10g: ", xmin + (i+0.5)*dx);
    for ( j = 0; j < 60*h[i]/m; j++ ) fprintf(stderr, "*");
    fprintf(stderr, "\n");
  }
}

/* generate a histogram from the command and the input file
 * write the output to stdout */
static int mkhist(const char *command, const char *fnin)
{
  token_t *px, *py = NULL, *pz = NULL, *pw = NULL;
  char *sx, *sy = NULL, *sz = NULL, *sw = NULL;
  char *sxyz, *shist = NULL, *p;
  int dim = 1, i, n, col, maxcol = 0, ix, iy, iz, xn, yn, zn, hn;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  double z, zmin, zmax, dz;
  double w, wtot = 0, norm, *hist = NULL;
  FILE *fp;
  double *data = NULL;
  size_t bufsz = 0;
  char *buf = NULL;
  const char *delims = " \t\r\n,";

  /* make a copy of the command */
  n = strlen(command);
  if ( (sxyz = calloc(n + 1, 1)) == NULL ) exit(-1);
  strcpy(sxyz, command);

  /* get the histogram string */
  if ( (p = strstr(sxyz, ">>")) != NULL ) {
    *p = '\0';
    shist = p + 2;
  } else {
    fprintf(stderr, "no >> in [%s]\n", command);
    return -1;
  }

  /* get the weight expression */
  if ( (p = strstr(sxyz, "::")) != NULL ) {
    *p = '\0';
    sw = p + 2;
    pw = parse2postfix(sw);
    if ( (i = getmaxcol(pw)) > maxcol ) maxcol = i;
  } else {
    fprintf(stderr, "no :: for weight in [%s], assuming 1.0\n", sxyz);
  }

  /* get the x expression */
  sx = sxyz;
  /* get the y expression */
  if ( (p = strchr(sx, ':')) != NULL ) {
    dim++;
    *p = '\0';
    sy = p + 1;
    /* get the z expression */
    if ( (p = strchr(sy, ':')) != NULL ) {
      dim++;
      *p = '\0';
      sz = p + 1;
    }
  }

  /* get the postfix expressions for x, y, z */
  px = parse2postfix(sx);
  if ( (i = getmaxcol(px)) > maxcol ) maxcol = i;
  if ( dim >= 2 ) {
    py = parse2postfix(sy);
    if ( (i = getmaxcol(py)) > maxcol ) maxcol = i;
    if ( dim >= 3 ) {
      pz = parse2postfix(sz);
      if ( (i = getmaxcol(pz)) > maxcol ) maxcol = i;
    }
  }

  /* allocate the data array */
  if ( (data = calloc(maxcol + 1, sizeof(*data))) == NULL ) exit(-1);

  /* get histogram parameters */
  if (3 != sscanf(shist, " h ( %d, %lf, %lf%n", &xn, &xmin, &xmax, &i)) {
    fprintf(stderr, "invalid histogram x-parameters [%s]\n", p);
    goto END;
  }
  p = shist + i;
  dx = (xmax - xmin) / xn;
  /* set the default parameters for y and z the same as those for x */
  zn   = yn   = xn;
  zmin = ymin = xmin;
  zmax = ymax = xmax;
  dz   = dy   = dx;

  if ( dim >= 2 ) { /* scan histogram y parameters */
    while ( *p && isspace(*p) ) p++;
    if ( *p != ')' ) {
      p++; /* skip a , or ; */
      if (3 != sscanf(p, "%d, %lf, %lf%n", &yn, &ymin, &ymax, &i)) {
        fprintf(stderr, "invalid histogram y-parameters [%s]\n", p);
        goto END;
      }
      p += i;
      dy = (ymax - ymin) / yn;
      /* set the default parameters for z the same as those for y */
      zn   = yn;
      zmin = ymin;
      zmax = ymax;
      dz   = dy;
      if ( dim >= 3 ) { /* scan histogram z parameters */
        while ( *p && isspace(*p) ) p++;
        if ( *p != ')' ) {
          p++; /* skip a , or ; */
          if (3 != sscanf(p, "%d, %lf, %lf%n", &zn, &zmin, &zmax, &i)) {
            fprintf(stderr, "invalid histogram z-parameters [%s]\n", p);
            goto END;
          }
          p += i;
          dz = (zmax - zmin) / zn;
        }
      }
    }
  }

  if ( dim == 1 ) {
    fprintf(stderr, "1D: x %s, weight %s, x %d bins (%g, %g)\n",
        sx, (sw ? sw : "1"), xn, xmin, xmax);
  } else if ( dim == 2 ) {
    fprintf(stderr, "2D: x %s, y %s, weight %s, x %d bins (%g, %g); y %d bins (%g, %g)\n",
        sx, sy, (sw ? sw : "1"), xn, xmin, xmax, yn, ymin, ymax);
  } else if ( dim == 3 ) {
    fprintf(stderr, "3D: x %s, y %s, z %s, weight %s, x %d bins (%g, %g); y %d bins (%g, %g); z %d bins (%g, %g)\n",
        sx, sy, sz, (sw ? sw : "1"), xn, xmin, xmax, yn, ymin, ymax, zn, zmin, zmax);
  }

  hn = xn;
  if ( dim >= 2 ) hn *= yn;
  if ( dim >= 3 ) hn *= zn;
  if ( (hist = calloc(hn, sizeof(*hist))) == NULL) exit(-1);
  for ( i = 0; i < hn; i++ ) hist[i] = 0;

  if ((fp = fopen(fnin, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnin);
    goto END;
  }

  /* read file line by line */
  while ( mygetline(&buf, &bufsz, fp) ) {
    if (buf[0] == '#' || buf[0] == '!') continue;
    p = strtok(buf, delims);
    for ( col = 0; col < maxcol && p != NULL; col++ ) {
      data[col] = atof(p);
      p = strtok(NULL, delims);
      //fprintf(stderr, "%g ", data[col]);
    }

    /* compute the index */
    x = evalpostfix(px, data, col);
    ix = ( x >= xmin ) ? ( x - xmin ) / dx : -1;
    i = ix;
    if ( ix < xn ) {
      if ( dim >= 2 ) {
        i *= yn;
        y = evalpostfix(py, data, col);
        iy = ( y >= ymin ) ? ( y - ymin ) / dy : -1;
        if ( iy < yn ) {
          i += iy;
          if ( dim >= 3 ) {
            i *= zn;
            z = evalpostfix(pz, data, col);
            iz = ( z >= zmin ) ? ( z - zmin ) / dz : -1;
            i = ( iz < zn ) ? i + iz : -1;
          }
        } else { i = -1; }
      }
    } else { i = -1; }

    if ( i >= 0 ) {
      w = ( pw != NULL ) ? evalpostfix(pw, data, col) : 1.0;
      //fprintf(stderr, " |  %d %g\n", i, w);
      hist[i] += w;
      wtot += w;
    }
  }

  /* print the histogram */
  if ( dim == 1 ) { /* 1D histogram */
    norm = 1.0/(wtot*dx);
    /* truncate the edges */
    for ( i = 0; i < xn; i++ ) {
      printf("%g %g\n", xmin + (i + 0.5) * dx, hist[i]*norm);
    }
    showhist(hist, xmin, dx, xn);
  } else if ( dim == 2 ) { /* 2D histogram */
    norm = 1.0/(wtot*dx*dy);
    for ( i = 0, ix = 0; ix < xn; ix++ ) {
      for ( iy = 0; iy < yn; iy++, i++ ) {
        printf("%g %g %g\n", xmin + (ix + 0.5) * dx,
            ymin + (iy + 0.5) * dy, hist[i]*norm);
      }
      printf("\n");
    }
  } else if ( dim == 3 ) { /* 3D histogram */
    norm = 1.0/(wtot*dx*dy*dz);
    for ( i = 0, ix = 0; ix < xn; ix++ ) {
      for ( iy = 0; iy < yn; iy++, i++ ) {
        for ( iz = 0; iz < zn; iz++, i++ ) {
          printf("%g %g %g %g\n", xmin + (ix + 0.5) * dx,
              ymin + (iy + 0.5) * dy, zmin + (iz + 0.5) * dz,
              hist[i]*norm);
        }
      }
      printf("\n");
    }
  }

  fclose(fp);
END:
  free(hist);
  free(buf);
  free(data);
  free(sxyz);
  free(px); free(py); free(pz); free(pw);
  return 0;
}


/* generate a random test file, if not exists */
static void genrandfile(const char *fn) {
  FILE *fp;
  int i, j;
  double x, y;

  if ( (fp = fopen(fn, "r")) != NULL ) {
    fclose(fp);
    return;
  }
  srand(clock());
  if ((fp = fopen(fn, "w")) != NULL) {
    for ( i = 0; i < 10000; i++ ) {
      for ( x = 0, j = 0; j < 10; j++ )
        x += 1.0*rand()/RAND_MAX;
      x /= 10;
      for ( y = 0, j = 0; j < 10; j++ )
        y += 1.0*rand()/RAND_MAX;
      y /= 10;
      fprintf(fp, "%g %g %g\n", x*2.19, x*y*6.19, y*y*9.19);
    }
    fclose(fp);
  }
}

int main(int argc, char **argv)
{
  const char *command = "sin(-$1+.05e1*(if($2>$3,$2,$3)))::1.0>>h(100, -1.0, 1.0)";
  //const char *command = "sin(-$1+$2):cos($3)>>h(100, -1.0, 1.0)";
  const char *fnin = "data.dat";

  if ( argc >= 2 ) {
    command = argv[1];
  }

  if ( argc >= 3 ) {
    fnin = argv[2];
  } else { /* otherwise generate a test data file */
    genrandfile(fnin);
  }
  mkhist(command, fnin);
  return 0;
}

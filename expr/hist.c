#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, DATA, VARIABLE };

const char *ops = ",+-*/%^_()";

#define VARNAME_MAX 128

typedef struct {
  int type;
  char s[VARNAME_MAX]; /* operator/function name */
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
    if ( *s == '*' && s[1] == '*' )
      t->s[0] = '^'; /* convert "**" to "^" */
    p = s + 1;
  } else if ( isdigit(*s) ) {
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
__inline static char *token2str(char *s, token_t *tok)
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

  if ( prtable['+'] == 0 ) { /* initialize the precedence table */
    prtable['('] = prtable[')'] = 16;
    prtable['_'] = 15;
    prtable['^'] = 13;
    prtable['*'] = prtable['/'] = prtable['%'] = 12;
    prtable['+'] = prtable['-'] = 11;
    prtable[','] = 1;
  }

  /* functions take the highest precedence */
  return isalpha(*s) ? 1000 : prtable[(int)(*s)];
}

static int isleftassoc(const char *s)
{
  if ( *s == '^' ) return 0;
  return 1;
}

/* convert an expression to a stack of postfix expression
 * Shunting Yard algorithm:
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

  /* where there are tokens to be read, read a token */
  for ( p = s; (p = gettoken(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == DATA || tok->type == VARIABLE ) {
      copytoken(pos++, tok);
    } else if ( tok->s[0] == ',' ) {
      /* do nothing for the comma */
    } else if ( tok->s[0] == '(' || tok->type == FUNCTION ) {
      copytoken(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a left bracket
       * and is not a function */
      while ( top->s[0] != '(' && !isalpha(top->s[0]) ) {
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
                 && top->s[0] != '(' && !isalpha(top->s[0]) ) {
          /* pop operators from the operator stack onto the output queue */
          copytoken(pos++, top--); /* pos = top */
          if ( top <= ost ) break;
        }
        //printf("pushing tok %s\n", tok->s);
        copytoken(++top, tok); /* top = token */
      }
    }
    copytoken(tok+1, tok); /* make a copy */
  }

  /* where there are still operator tokens on the stack
   * pop the operator onto the output queue */
  while ( top > ost )
    copytoken(pos++, top--);

  pos->type = NULLTYPE;
  free(ost);
  return que;
}

static double max(double a, double b) { return (a > b) ? a : b; }
static double min(double a, double b) { return (a < b) ? a : b; }
static double iif(double c, double a, double b) { return (c == 0) ? a : b; }

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
static double evalpostfix(token_t *que, const double *arr)
{
  int i, n, top;
  double *st, ans;
  token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ((st = calloc(n, sizeof(*st))) == NULL) exit(-1);
  top = 0;

  for ( pos = que; pos->type != NULLTYPE; pos++ ) {
    if ( pos->type == NUMBER ) {
      st[top++] = pos->val;
    } else if ( pos->type == DATA ) {
      st[top++] = arr[pos->col - 1]; /* array index off-by-one */
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

static int mkhist(const char *command, const char *fnin)
{
  token_t *px, *pw;
  char *sx, *sw, *p, *q;
  int i, n, maxcol, xn;
  double x, w, xmin, xmax, dx, wtot = 0, norm, *hist;
  FILE *fp;
  double *data;
  size_t bufsz = 0;
  char *buf = NULL;
  const char *delims = " \t\r\n";

  /* get the x value */
  if ( (p = strstr(command, ":")) == NULL ) {
    fprintf(stderr, "no : in [%s]\n", command);
    return -1;
  }
  n = p - command;
  if ( (sx = calloc(n + 1, 1)) == NULL ) exit(-1);
  strncpy(sx, command, n + 1);
  sx[n] = '\0';
  px = parse2postfix(sx);
  maxcol = getmaxcol(px);

  /* get the weight */
  if ( (q = strstr(p, "::")) != NULL ) {
    p = q + 2;
    q = strstr(p, ">>");
    n = q - p;
    if ( (sw = calloc(n + 1, 1)) == NULL ) exit(-1);
    strncpy(sw, p, n + 1);
    sw[n] = '\0';
    pw = parse2postfix(sw);
    if ((i = getmaxcol(pw)) > maxcol) maxcol = i;
  } else {
    fprintf(stderr, "no :: for weight in [%s], assuming 1.0\n", command);
    sw = "1.0";
    pw = NULL;
  }

  /* allocate the data array */
  if ( (data = calloc(maxcol + 1, sizeof(*data))) == NULL ) exit(-1);

  /* get histogram parameters */
  if ( (p = strstr(p, "h(")) == NULL ) {
    fprintf(stderr, "no histogram function h(...) in [%s]\n", command);
    return -1;
  }
  if (3 != sscanf(p, " h ( %d, %lf, %lf )", &xn, &xmin, &xmax)) {
    fprintf(stderr, "invalid histogram string [%s]\n", p);
    goto END;
  }

  fprintf(stderr, "data %s, weight %s, # of bins %d, min %g, max %g\n",
      sx, sw, xn, xmin, xmax);
  if ( (hist = calloc(xn, sizeof(*hist))) == NULL) exit(-1);
  for ( i = 0; i < xn; i++ ) hist[i] = 0;
  dx = (xmax - xmin) / xn;

  if ((fp = fopen(fnin, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnin);
    goto END;
  }

  /* read file line by line */
  while ( mygetline(&buf, &bufsz, fp) ) {
    if (buf[0] == '#' || buf[0] == '!') continue;
    p = strtok(buf, delims);
    for ( i = 0; i < maxcol && p != NULL; i++ ) {
      data[i] = atof(p);
      p = strtok(NULL, delims);
      //fprintf(stderr, "%g ", data[i]);
    }
    x = evalpostfix(px, data);
    w = (pw != NULL) ? evalpostfix(pw, data) : 1.0;
    //fprintf(stderr, " |  %g %g\n", x, w);
    if ( x >= xmin ) {
      i = ( x - xmin ) / dx;
      if ( i < xn ) {
        hist[i] += w;
        wtot += w;
      }
    }
  }

  norm = 1.0/(wtot*dx);
  /* truncate the edges */
  i = 0;
  for ( ; i < xn; i++ ) if ( hist[i] > 0 ) break;
  for ( ; xn > i; xn-- ) if ( hist[xn-1] > 0 ) break;
  for ( ; i < xn; i++ ) {
    printf("%g %g\n", xmin + (i + 0.5) * dx, hist[i]*norm);
  }

  fclose(fp);
  free(hist);
  free(buf);
END:
  free(px);
  free(pw);
  return 0;
}


/* generate a random test file, if not exists */
static void genrandfile(const char *fn) {
  FILE *fp;
  int i, j;
  double x;

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
      fprintf(fp, "%g %g %g\n", x*2.19, x*x*6.19, x*x*x*9.19);
    }
    fclose(fp);
  }
}

int main(int argc, char **argv)
{
  const char *command = "sin(-$1+$2)::1.0>>h(100, -1.0, 1.0)";
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

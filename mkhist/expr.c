#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, DATA, VARIABLE };

const char *ops = ",+-*/%^_()";

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
    }
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

__inline static void test_gettoken(char *s)
{
  char *p, *q;
  token_t tok[1];
  int i;

  printf("Expression: %s\n", s);
  for ( p = s, i = 0; ; i++ ) {
    q = gettoken(tok, p);
    if ( q == NULL ) break;
    p = q;
    printf("%4d: %4s, %8g, %4d\n", i, tok->s, tok->val, tok->col);
  }
}


static void copytoken(token_t *a, token_t *b)
{
  a->type = b->type;
  strncpy(a->s, b->s, VARNAME_MAX);
  a->val = b->val;
  a->col = b->col;
}

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

#ifdef DEBUG
__inline static void sy_print(token_t *tok, token_t *que, token_t *pos, token_t *ost, token_t *top)
{
  token_t *t;
  char s[VARNAME_MAX];

  if ( tok != NULL ) {
    printf("Token: %s\n", token2str(s, tok));
  }
  printf("Queue(%ld): ", pos - que);
  for ( t = que; t < pos; t++ ) printf("%s ", token2str(s, t));
  printf("\nStack(%ld): ", top - ost);
  for ( t = top; t > ost; t-- ) printf("%s ", token2str(s, t));
  printf("\n\n");
  getchar();
}
#endif

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

  /* when there are tokens to be read, read a token */
  for ( p = s; (p = gettoken(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == DATA || tok->type == VARIABLE ) {
      copytoken(pos++, tok);
    } else if ( tok->s[0] == '(' || tok->type == FUNCTION ) {
      copytoken(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a "("
       * or function */
      while ( top->s[0] != '(' && tok->type != FUNCTION ) {
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
#ifdef DEBUG
    sy_print(tok, que, pos, ost, top);
#endif
  }

  /* where there are still operator tokens on the stack
   * pop the operator onto the output queue */
  while ( top > ost )
    copytoken(pos++, top--);
#ifdef DEBUG
  sy_print(NULL, que, pos, ost, top);
#endif

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
static double evalpostfix(token_t *que, const double *arr)
{
  int i, n, top;
  double *st, ans;
  token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ( (st = calloc(n, sizeof(*st))) == NULL ) exit(-1);
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
      } else if ( pos->s[0] == ',' ) {
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
#ifdef DEBUG
    {
      int j;
      char s[VARNAME_MAX];
      printf("%10s: ", token2str(s, pos));
      for ( j = 0; j < top; j++ ) printf("%g ", st[j]);
      printf("\n");
      //getchar();
    }
#endif
  }
  ans = st[0];
  free(st);
  return ans;
}

int main(void)
{
  char *expr1 = "pow(1.0+2.0, 3)"; // "sin(1.0e-3 * $1) -3 + 4 * -2 / ( -1 - 1 ) ^ -2 ^ +2";
  char *expr2 = "sin ( max ( 2, -$1 + 4 ) / 3 * 3.1415 )";
  token_t *pexpr;
  double data[] = {1.1, 2.2, 3.3}, y;

  //test_gettoken(expr2);
#if 1
  pexpr = parse2postfix(expr1);
  y = evalpostfix(pexpr, data);
  printf("%s = %g\n", expr1, y);
  free(pexpr);
#endif

#if 1
  /* construct a postfix expression from the input */
  pexpr = parse2postfix(expr2);
  /* evaluate the postfix expression */
  y = evalpostfix(pexpr, data);
  printf("%s = %g\n", expr2, y);
  free(pexpr);
#endif

  return 0;
}

#define DEBUG 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define VARNAME_MAX 128

static double dblneg(double a) { return -a; }
static double dbladd(double a, double b) { return a + b; }
static double dblsub(double a, double b) { return a - b; }
static double dblmul(double a, double b) { return a * b; }
static double dbldiv(double a, double b) { return a / b; }
static double dblnot(double a) { return a == 0; }
static double dbllt(double a, double b) { return a < b; }
static double dblle(double a, double b) { return a <= b; }
static double dblgt(double a, double b) { return a > b; }
static double dblge(double a, double b) { return a >= b; }
static double dbleq(double a, double b) { return a == b; }
static double dblneq(double a, double b) { return a != b; }
static double dbland(double a, double b) { return (a != 0) && (b != 0); }
static double dblor(double a, double b)  { return (a != 0) || (b != 0); }
static double dblmax(double a, double b) { return ( a > b ) ? a : b; }
static double dblmin(double a, double b) { return ( a < b ) ? a : b; }
static double iif(double c, double a, double b) { return ( c != 0 ) ? a : b; }

enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, COLUMN, VARIABLE };

typedef struct {
  const char *s;
  double (*f)(); /* corresponding function */
  int nops; /* number of operands */
  int preced; /* precedence */
  int rightassoc; /* right associative */
} opinfo_t;

/* supported operators: */
opinfo_t oplist[] = {
  { "",      NULL, 0,  0, 0 }, /* root of the operator stack */
  { "(",     NULL, 0, 16, 0 }, { ")",     NULL, 0, 16, 0 },
  { "_",   dblneg, 1, 15, 0 }, /* negative */
  { "!",   dblnot, 1, 15, 0 }, { "not", dblnot, 1, 15, 0 },
  { "^",      pow, 2, 13, 1 }, { "**",     pow, 2, 13, 1 }, /* power */
  { "*",   dblmul, 2, 12, 0 }, { "/",   dbldiv, 2, 12, 0 }, { "%",     fmod, 2, 12, 0 },
  { "+",   dbladd, 2, 11, 0 }, { "-",   dblsub, 2, 11, 0 },
  { "<",    dbllt, 2,  9, 0 }, { ">",    dblgt, 2,  9, 0 },
  { "lt",   dbllt, 2,  9, 0 }, { "gt",   dblgt, 2,  9, 0 },
  { "<=",   dblle, 2,  9, 0 }, { ">=",   dblge, 2,  9, 0 },
  { "le",   dblle, 2,  9, 0 }, { "ge",   dblge, 2,  9, 0 },
  { "==",   dbleq, 2,  8, 0 }, { "eq",   dbleq, 2,  8, 0 },
  { "!=",  dblneq, 2,  8, 0 }, { "<>",  dblneq, 2,  8, 0 }, { "neq", dblneq, 2,  8, 0 },
  { "&&",  dbland, 2,  5, 0 }, { "and", dbland, 2,  5, 0 },
  { "||",   dblor, 2,  4, 0 }, { "or",   dblor, 2,  4, 0 },
  { ",",     NULL, 0,  1, 0 },
  { NULL,    NULL, 0,  0, 0 } };

typedef struct {
  const char *s;
  double (*f)();
  int narg;
} funcmap_t;

funcmap_t funclist[] = {
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
  {"min", dblmin, 2},
  {"max", dblmax, 2},
  {"if", iif, 3},
  {"iif", iif, 3},
  {NULL, NULL, 0}
};

typedef struct token_s {
  int type;
  char s[VARNAME_MAX]; /**< operator/function name */
  double val;
  int col;             /**< column id in the data */
  opinfo_t *op;        /**< pointer to operator information */
  funcmap_t *func;     /**< pointer to function information */
} token_t;


/* get a token from the string and return a pointer
 * after the token */
static const char *gettoken(token_t *t, const char *s)
{
  char *p;
  int i;
  size_t oplen;
  opinfo_t *op;

  /* skip leading spaces */
  while ( *s && isspace(*s) ) s++;
  if ( *s == '\0' ) return NULL;

  t->s[0] = '\0';
  t->val = 0;
  t->col = 0;

  if ( isdigit(*s) || (*s == '.' && isdigit(s[1])) ) { /* .5 is understood as a number */
    t->type = NUMBER;
    t->val = strtod(s, &p);
  } else if ( *s == '$' ) {
    t->type = COLUMN;
    t->col = (int)strtol(s + 1, &p, 10);
  } else {
    /* scan for operators */
    for ( oplen = 3; oplen > 0; oplen-- ) {
      for ( op = oplist + 1; op->s != NULL; op++ ) {
        if ( strlen(op->s) == oplen && strncmp(op->s, s, oplen) == 0 ) {
          t->type = OPERATOR;
          t->op = op;
          strcpy(t->s, op->s);
          return s + oplen;
        }
      }
    }
    if ( isalpha(*s) ) {
      for ( i = 0; isalnum(s[i]) && i < VARNAME_MAX - 1; i++ )
        t->s[i] = s[i];
      t->s[i] = '\0';
      t->type = VARIABLE;
      return s + i;
    } else {
      fprintf(stderr, "Error: unknown token [%s]\n", s);
      exit(1);
    }
  }
  return (const char *) p;
}

/* string representation of token */
__inline static char *token2str(char *s, const token_t *tok)
{
  if ( tok->type == NUMBER ) {
    sprintf(s, "%g", tok->val);
  } else if ( tok->type == OPERATOR || tok->type == FUNCTION || tok->type == VARIABLE ) {
    sprintf(s, "%s", tok->s);
  } else if ( tok->type == COLUMN ) {
    sprintf(s, "$%d", tok->col);
  }
  return s;
}

__inline static void test_gettoken(const char *s)
{
  const char *p, *q;
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


/** a = b */
static void token_copy(token_t *a, const token_t *b)
{
  a->type = b->type;
  strncpy(a->s, b->s, VARNAME_MAX);
  a->val = b->val;
  a->col = b->col;
  a->op = b->op;
  a->func = b->func;
}

/** set token `a' as operator `s' */
static void set_op(token_t *a, const char *s)
{
  opinfo_t *op;
  a->type = OPERATOR;
  for ( op = oplist; op->s != NULL; op++ ) {
    if ( strcmp(op->s, s) == 0 ) {
      strcpy(a->s, s);
      a->op = op;
      break;
    }
  }
}

/** convert a variable to a function */
static void make_func(token_t *t)
{
  funcmap_t *func;

  t->type = FUNCTION;
  for ( func = funclist; func->s != NULL; func++ ) {
    if ( strcmp(func->s, t->s) == 0 ) {
      t->func = func;
      return;
    }
  }
  fprintf(stderr, "Error: unknown function [%s]\n", t->s);
  exit(1);
}

#ifdef DEBUG
__inline static void sy_print(const token_t *tok, token_t *que, token_t *pos, token_t *ost, token_t *top)
{
  const token_t *t;
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


/** Convert an expression to a stack of postfix expression.
 * The shunting-yard algorithm:
 * https://en.wikipedia.org/wiki/Shunting-yard_algorithm
 *
 * @return a malloc-ed queue stack of tokens.
 */
static token_t *parse2postfix(const char *s)
{
  token_t *que; /* output queue */
  token_t *ost; /* operator stack */
  token_t *pos, *top, tok[2];
  const char *p;
  int n;

  n = strlen(s);

  /* initialize the output queue */
  if ( (que = calloc(n + 1, sizeof(*que))) == NULL ) exit(-1);
  pos = que;

  /* initialize the operator stack */
  if ( (ost = calloc(n + 1, sizeof(*ost))) == NULL ) exit(-1);
  top = ost;
  /* initialize the operator stack with a null operator */
  set_op(top, "");
  token_copy(&tok[1], top); /* tok[1] is the previous token */

  /* when there are tokens to be read, read a token */
  for ( p = s; (p = gettoken(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == COLUMN || tok->type == VARIABLE ) {
      /* if the token is a number, then push it to the output queue */
      token_copy(pos++, tok);
    } else if ( tok->s[0] == '(' ) {
      if ( tok[1].type == VARIABLE ) {
        /* convert the preceding variable to a function */
        make_func(--pos);
        /* move the function onto the operator stack */
        token_copy(++top, pos);
      }
      /* push the "(" onto the operator stack */
      token_copy(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a "(" */
      while ( top->s[0] != '(' ) {
        /* pop operators from the operator stack onto the output queue */
        token_copy(pos++, top--);
        if ( top <= ost ) break;
      }
      if ( top->s[0] == '(' ) {
        top--; /* pop the left bracket from the stack */
        if ( top->type == FUNCTION ) {
          /* pop the function from the operator stack onto the output queue */
          token_copy(pos++, top--);
        }
      } else {
        token_copy(pos++, top--);
      }
    } else if ( tok->type == OPERATOR ) {
      if ( (tok->s[0] == '+' || tok->s[0] == '-') && tok->s[1] == '\0'
           && ( tok[1].type == OPERATOR && tok[1].s[0] != ')' ) ) {
        /* handle a unary operator */
        if ( tok->s[0] == '-' ) {
          set_op(tok, "_"); /* change it to negation */
          token_copy(++top, tok); /* top = token */
        } /* skip the unary '+' */
      } else {
        while ( ((top->op->preced > tok->op->preced)
                 || (top->op->preced == tok->op->preced && !top->op->rightassoc))
                && top->s[0] != '(' ) {
          /* pop operators from the operator stack onto the output queue */
          token_copy(pos++, top--); /* pos = top */
          if ( top <= ost ) break;
        }
        /* push the read operator onto the operator stack */
        token_copy(++top, tok); /* top = token */
      }
    }
#ifdef DEBUG
    sy_print(tok, que, pos, ost, top);
#endif
    token_copy(tok+1, tok); /* make a copy */
  }

  /* while there are still operator tokens on the stack,
   * pop the operator onto the output queue */
  while ( top > ost )
    token_copy(pos++, top--);
#ifdef DEBUG
    sy_print(tok, que, pos, ost, top);
#endif

  pos->type = NULLTYPE;
  free(ost);
  return que;
}

typedef struct {
  char s[VARNAME_MAX];
  double val;
} varmap_t;

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
  int i, n, top, narg;
  double *st, ans, (*f)();
  const token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ( (st = calloc(n, sizeof(*st))) == NULL ) exit(-1);
  top = 0;

  for ( pos = que; pos->type != NULLTYPE; pos++ ) {
    if ( pos->type == NUMBER ) {
      st[top++] = pos->val;
    } else if ( pos->type == COLUMN ) {
      st[top++] = ( pos->col <= narr ) ? arr[pos->col] : 0.0;
    } else if ( pos->type == VARIABLE ) {
      for ( i = 0; varmap[i].s[0] != '\0'; i++ ) {
        if ( strcmp(varmap[i].s, pos->s) == 0 ) {
          st[top++] = varmap[i].val;
          break;
        }
      }
    } else if ( pos->type == OPERATOR || pos->type == FUNCTION ) {
      if ( pos->type == OPERATOR ) {
        narg = pos->op->nops;
        f = pos->op->f;
      } else {
        narg = pos->func->narg;
        f = pos->func->f;
      }

      if ( f != NULL ) { /* NOTE: operators like "," have no function */
        if ( top < narg ) {
          fprintf(stderr, "Error: insufficient arguments for [%s], has %d, require %d\n",
                  pos->s, top, narg);
          exit(1);
        }
        top -= narg - 1;
        if ( narg == 0 ) {
          st[top-1] = (*f)();
        } else if ( narg == 1 ) {
          st[top-1] = (*f)(st[top-1]);
        } else if ( narg == 2 ) {
          st[top-1] = (*f)(st[top-1], st[top]);
        } else if ( narg == 3 ) {
          st[top-1] = (*f)(st[top-1], st[top], st[top+1]);
        }
      }
    }
#ifdef DEBUG
    { /* print out the evaluation stack */
      int j;
      char s[VARNAME_MAX];
      fprintf(stderr, "%10s: ", token2str(s, pos));
      for ( j = 0; j < top; j++ ) fprintf(stderr, "%g ", st[j]);
      fprintf(stderr, "\n");
      //getchar();
    }
#endif
  }
  ans = ( top > 0 ) ? st[top-1] : 0.0;
  free(st);
  return ans;
}

int main(void)
{
  char *expr1 = "max (sin($1<$2&&$2<$3),0.5)";
  //char *expr1 = "(1+3 != 2) && !(2 == 3)"; // "pow(1.0+2.0, 3)"; // "sin(1.0e-3 * $1) -3 + 4 * -2 / ( -1 - 1 ) ^ -2 ^ +2";
  //char *expr2 = "sin ( max ( 2, -$1 + 4 ) / 3 * 3.1415 )";
  char *expr2 = "sin ( if ( 2 > -$1 +4, 2, -$1 ** -1 ) / 3 * 3.1415 )";
  token_t *pexpr;
  double data[] = {1, 1.1, 2.2, 3.3}, y;

  //test_gettoken(expr2);
#if 0
  pexpr = parse2postfix(expr1);
  y = evalpostfix(pexpr, data, 3);
  printf("%s = %g\n", expr1, y);
  free(pexpr);
#endif

#if 1
  /* construct a postfix expression from the input */
  pexpr = parse2postfix(expr2);
  /* evaluate the postfix expression */
  y = evalpostfix(pexpr, data, 3);
  printf("%s = %g\n", expr2, y);
  free(pexpr);
#endif

  return 0;
}

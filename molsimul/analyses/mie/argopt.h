#ifndef ARGOPT_H__
#define ARGOPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "util.h"


/* compare two strings ignoring cases */
__inline static int strcmpnc(const char *s, const char *t)
{
  int i, cs, ct;

  if (s == NULL || t == NULL) return 0;
  for (i = 0; ; i++) {
    cs = s[i];
    ct = t[i];
    if (cs == 0 || ct == 0) break;
    cs = toupper( (unsigned char) cs );
    ct = toupper( (unsigned char) ct );
    if (cs != ct) break;
  }
  return cs - ct;
}



enum { OPT_ARG, OPT_OPT, OPT_CFG, OPT_COUNT };



/* option either from arguments or configuration */
typedef struct {
  int isopt; /* one of OPT_xxx values */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *key; /* key, for cfg files as in `key = val' */
  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  int initval; /* initial value, for a switch option */
  const char **sarr; /* array of string values, for a list option */
  int scnt; /* length of the array, for a list option */
  unsigned flags;
} opt_t;


/* return the index of string from a predefined array
 * using fuzzy string comparison */
__inline static int opt_select(const char *s, const char **sarr, int n)
{
  int i;

  for ( i = 0; i < n; i++ )
    if ( strcmpfuzzy(sarr[i], s) == 0 )
      return i;
  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i >= 0 && i < n ) return i;
  }
  fprintf(stderr, "Error: cannot select %s from the array of %d items:",
      s, n);
  for ( i = 0; i < n; i++ )
    fprintf(stderr, "%s, ", sarr[i]);
  fprintf(stderr, "\n");
  exit(1);
  return 0;
}



#define OPT_MUST     0x0001  /* a mandatory argument or option */
#define OPT_SWITCH   0x0002  /* an option is a switch */
#define OPT_SET      0x0004  /* an argument/option is set */

/* translate string value in `o->val' into
 * actual ones (o->ptr) through sscanf(), etc */
__inline static int opt_getval(opt_t *o)
{
  const char *fmt = o->fmt;

  /* for a string argument, it can be obtained from NULL or "%s"
   * NULL: string memory from command-line */
  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
    *((const char **) o->ptr) = o->val;
  } else if (strcmpnc(fmt, "%list") == 0
          || strcmpnc(fmt, "%enum") == 0) {
    *((int *) o->ptr) = opt_select(o->val, o->sarr, o->scnt);
  } else if (strcmp(fmt, "%b") == 0) { /* switch */
    /* switch the default value */
    if (o->flags & OPT_SET) return !o->initval;
    else return o->initval;
  } else if (strcmp(fmt, "%+") == 0) { /* incremental, like -vv */
    *((int *) o->ptr) += 1;
  } else if (strcmp(fmt, "%-") == 0) { /* decremental, like -vv */
    *((int *) o->ptr) -= 1;
  } else { /* call sscanf */
    if (1 != sscanf(o->val, fmt, o->ptr)) {
      fprintf(stderr, "Error: unable to convert a value for [%s] as fmt [%s], raw string: [%s]\n",
          o->desc, fmt, o->val);
      return 1;
    }
  }
  return 0;
}



/* register an option
 *
 * for a configure entry, set `key' and leave `sflag' = NULL
 * for a command-line option, set `sflag' and leave `key' = NULL
 * `fmt' is the sscanf() format string
 * `*ptr' is the target variable
 * `fmt' can "%b" for a switch (like an command-line option "-v")
 * `fmt' can have a prefix `!' to mean a mandatory option
 * both NULL and "%s" of `fmt' mean string values, the type of
 *  `ptr' should be `char **', the difference is that `*ptr'
 *  is directly assigned to `o->val' during opt_getval() in the
 *  former case, but extra memory is allocated to copy `o->val'
 *  in the latter case */
__inline static void opt_set(opt_t *o, const char *sflag, const char *key,
    const char *fmt, void *ptr, const char *desc,
    const char **sarr, int scnt)
{
  o->ch = '\0';
  if ( key != NULL ) { /* cfg file `key = val', not a command-line argument */
    o->isopt = OPT_CFG;
  } else if ( sflag != NULL ) { /* option */
    o->isopt = OPT_OPT;
    o->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    o->isopt = OPT_ARG;
  }
  o->sflag = sflag;
  o->key = key;
  o->flags = 0;
  if (ptr == NULL) {
    fprintf(stderr, "null pass to opt with %s: %s\n", sflag, desc);
    exit(1);
  }
  o->ptr = ptr;
  if (fmt == NULL) fmt = "";
  if (fmt[0] == '!') {
    fmt++;
    o->flags |= OPT_MUST;
  }
  if (fmt[0] != '\0' && fmt[0] != '%') {
    fprintf(stderr, "unknown format (missing `%%') flag `%s\', fmt `%s', description: %s\n",
        sflag, fmt, desc);
    exit(1);
  }
  if (strcmp(fmt, "%b") == 0) {
    o->flags |= OPT_SWITCH;
    o->initval = *((int *) ptr); /* save the initial value */
  }
  o->fmt = fmt;
  o->pfmt = NULL;
  o->desc = desc;
  o->sarr = sarr;
  o->scnt = scnt;
}



/* print the value of o->ptr */
#define opt_printptr(o) opt_fprintptr(stderr, o)
__inline static void opt_fprintptr(FILE *fp, opt_t *o)
{
  const char *fmt;

  for (fmt = o->fmt; *fmt && *fmt != '%'; fmt++) ;

#define ELIF_PF_(fm, fmp, type) \
  else if (strcmp(fmt, fm) == 0) \
  { fprintf(fp, (o->pfmt ? o->pfmt : fmp), *(type *)o->ptr); }

  if (fmt == NULL || *fmt == '\0' || strcmp(fmt, "%s") == 0) {
    fprintf(fp, "%s", (*(char **) o->ptr) ? (*(char **) o->ptr) : "NULL");
  } else if (strcmp(fmt, "%b") == 0) {
    fprintf(fp, "%s", (*(int *)o->ptr) ? "true" : "false");
  } else if (strcmpfuzzy(fmt, "%list") == 0
          || strcmpfuzzy(fmt, "%enum") == 0) {
    int ival = *((int *) o->ptr);
    fprintf(fp, "%d (%s)", ival, o->sarr[ival]);
  }
  ELIF_PF_("%+", "%d", int)
  ELIF_PF_("%-", "%d", int)
  ELIF_PF_("%d", "%d", int)
  ELIF_PF_("%u", "%u", unsigned)
  ELIF_PF_("%x", "0x%x", unsigned)
  ELIF_PF_("%ld", "%ld", long)
  ELIF_PF_("%lo", "%lo", long)
  ELIF_PF_("%lu", "%lu", unsigned long)
  ELIF_PF_("%lx", "0x%lx", unsigned long)
#if defined(HAVELONGLONG) /* C99 or GCC extension */
  ELIF_PF_("%lld", "%lld", long long)
  ELIF_PF_("%llo", "%llo", long long)
  ELIF_PF_("%llu", "%llu", unsigned long long)
  ELIF_PF_("%llx", "0x%llx", unsigned long long)
#endif
  ELIF_PF_("%f", "%g", float)
  ELIF_PF_("%lf", "%g", double)
  ELIF_PF_("%Lf", "%Lg", long double)
  else {
    fprintf(fp, "unknown %s-->%%d: %d", fmt, *(int *) o->ptr);
  }
}



/* search an option list, return an option whose variable address is p */
__inline static opt_t *opt_find(opt_t *ls, int n, const void *p)
{
  int i;

  for (i = 0; i < n; i++) if (ls[i].ptr == p) return ls + i;
  return NULL;
}



/* search an option list to see if an option is explicitly set */
__inline static int opt_isset(opt_t *ls, int n, const void *p, const char *var)
{
  opt_t *o = opt_find(ls, n, p);
  if ( o == NULL ) {
    fprintf(stderr, "cannot find var %s, ptr %p\n", var, p);
    exit(1);
  }
  return o->flags & OPT_SET ? 1 : 0;
}


typedef struct {
  int nopt;
  opt_t *opts;
  const char *prog;
  const char *desc;
  int version;
  unsigned flags;
  int dum_[4]; /* space holder */
} argopt_t;

#define ARGOPT_MUST     OPT_MUST    /* mandatory argument or option, format starts with ! */
#define ARGOPT_SWITCH   OPT_SWITCH  /* format "%b" */
#define ARGOPT_SET      OPT_SET



#define argopt_getopt(ao, p) opt_find(ao->opts, ao->nopt, p)
#define argopt_getarg argopt_getopt

/* test if argument/option is explicitly set */
#define argopt_isset(ao, var) opt_isset(ao->opts, ao->nopt, &var, #var)



/* initialize the argument structure */
__inline static argopt_t *argopt_open(unsigned flags)
{
  argopt_t *ao;

  if ( (ao = calloc(1, sizeof(*ao))) == NULL ) return NULL;
  ao->flags = flags;
  ao->nopt = 0;
  ao->opts = NULL;
  memset(ao->dum_, '\0', sizeof(ao->dum_));
  return ao;
}



__inline static void argopt_close(argopt_t *ao)
{
  if (ao->opts) free(ao->opts);
  free(ao);
}



/* print version and die */
__inline static void argopt_version(argopt_t *ao)
{
  fprintf(stderr, "%s: %s, version %d\n",
      ao->prog, ao->desc ? ao->desc : "", ao->version);
  argopt_close(ao);
  exit(1);
}



/* print help message and die */
__inline static void argopt_help(argopt_t *ao)
{
  int i, j, len, maxlen;
  opt_t *o;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  fprintf(stderr, "%s, version %d\n",
      ao->desc ? ao->desc : ao->prog, ao->version);
  fprintf(stderr, "USAGE\n  %s {OPTIONS}", ao->prog);
  for (i = 0; i < ao->nopt; i++) {
    const char *bra = "", *ket = "";
    o = ao->opts + i;
    if (o->isopt != OPT_ARG) continue;
    if (o->flags & OPT_MUST) {
      if (strchr(o->desc, ' '))
        bra = "[", ket = "]";
    } else
      bra = "{", ket = "}";
    fprintf(stderr, " %s%s%s", bra, o->desc, ket);
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "OPTIONS:\n") ;
  /* compute the width of the longest option */
  for (maxlen = 0, i = 0; i < ao->nopt; i++) {
    if (ao->opts[i].isopt == OPT_ARG) continue;
    len = strlen(ao->opts[i].sflag);
    if (len > maxlen) maxlen = len;
  }
  for (i = 0; i < ao->nopt; i++) {
    o = ao->opts + i;
    if (o->isopt == OPT_ARG) continue;
    desc = o->desc;
    if (strcmp(desc, "$HELP") == 0)
      desc = sysopt[0];
    else if (strcmp(desc, "$VERSION") == 0)
      desc = sysopt[1];
    fprintf(stderr, "  %-*s : %s%s%s", maxlen, o->sflag,
        ((o->flags & OPT_MUST) ? "[MUST] " : ""),
        (!(o->flags & OPT_SWITCH) ? "followed by " : ""), desc);
    if ( o->sarr != NULL ) {
      fprintf(stderr, ", options: ");
      for ( j = 0; j < o->scnt; j++ ) {
        fprintf(stderr, "%s", o->sarr[j]);
        if ( j < o->scnt - 1 ) fprintf(stderr, ", ");
      }
    }
    if (o->ptr && o->ptr != ao->dum_) { /* print default values */
      fprintf(stderr, ", default: ");
      opt_fprintptr(stderr, o);
    }
    fprintf(stderr, "\n");
  }
  argopt_close(ao);
  exit(1);
}



#define argopt_addhelp(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$HELP")
#define argopt_addversion(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$VERSION")
#define argopt_add(ao, sflag, fmt, ptr, desc) \
  argopt_addx(ao, sflag, fmt, ptr, desc, NULL, 0)

/* register an argument or option
 * sflag: string flag, or NULL for an argument
 * fmt: sscanf() format string, "%b" for a switch
 * return the index */
__inline static int argopt_addx(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc,
    const char **sarr, int scnt)
{
  opt_t *o;
  int n;

  n = ao->nopt++;
  if ( (ao->opts = realloc(ao->opts, sizeof(*ao->opts) * ao->nopt)) == NULL ) {
    exit(1);
  }
  o = ao->opts + n;
  opt_set(o, sflag, NULL, fmt, ptr, desc, sarr, scnt);
  return n;
}



/* main parser of arguments */
__inline static void argopt_parse(argopt_t *ao, int argc, char **argv)
{
  int i, j, k, ch, acnt = 0;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      /* find the `acnt`th argument from the `ol` list */
      for ( k = 0, j = 0; j < ao->nopt; j++ ) {
        if ( ol[j].isopt ) continue;
        if ( ++k > acnt ) break;
      }
      if (j >= ao->nopt) argopt_help(ao);
      ol[j].val = argv[i];
      ol[j].flags |= OPT_SET;
      if (0 != opt_getval(ol + j))
        argopt_help(ao);
      acnt++;
      continue;
    }

    /* it's an option, loop for abbreviated form "-abc" == "-a -b -c" */
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      int islong = (j == 1 && argv[i][1] == '-');

      if (islong) { /* match against long options */
        for (k = 0; k < ao->nopt; k++) {
          int lenf;
          if (ol[k].sflag == NULL) continue;
          lenf = strlen(ol[k].sflag);
          if (ol[k].isopt &&
              strncmp(argv[i], ol[k].sflag, lenf) == 0 &&
              ( (ol[k].flags & OPT_SWITCH)
              || strchr("= ", argv[i][lenf]) ) ) /* followed by a space or "=" */
            break;
        }
      } else { /* match against short options */
        for (k = 0; k < ao->nopt; k++)
          if (ol[k].isopt && ch == ol[k].ch)
            break;
      }
      if (k >= ao->nopt) {
        fprintf(stderr, "cannot handle option [%s]\n", argv[i]);
        argopt_help(ao);
      }

      if (ol[k].desc[0] == '$') { /* system commands */
        if (strcmp(ol[k].desc, "$HELP") == 0)
          argopt_help(ao);
        else if (strcmp(ol[k].desc, "$VERSION") == 0)
          argopt_version(ao);
      }

      if (ol[k].flags & OPT_SWITCH) { /* handle switches "%b" */
        ol[k].flags |= OPT_SET;
        /* switch the default value, note this flag may be passed
         * several times, so we don't want to flip around */
        *((int *) ol[k].ptr) = !ol[k].initval;
        if (islong) break; /* go to the next argument argv[i+1] */
      } else if (strcmp(ol[k].fmt, "%+") == 0
              || strcmp(ol[k].fmt, "%-") == 0) {
        if (0 != opt_getval(ol + k)) argopt_help(ao);
      } else { /* look for the argument for this option */
        int hasv = 0;
        if (islong) { /* e.g., --version=11 */
          j = strlen(ol[k].sflag);
          if (argv[i][j] == '=') {
            ol[k].val = argv[i] + j + 1;
            hasv = 1;
          }
        } else { /* e.g., -n8 */
          if (argv[i][++j]) {
            ol[k].val = argv[i] + j;
            hasv = 1;
          }
        }

        if (!hasv) { /* --version 11 or -n 8 */
          if (++i >= argc) {
            fprintf(stderr, "%s(%s) requires an argument!\n", ol[k].sflag, argv[i - 1]);
            argopt_help(ao);
          }
          ol[k].val = argv[i];
        }
        ol[k].flags |= OPT_SET;
        if (0 != opt_getval(ol + k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & OPT_MUST) && !(ol[i].flags & OPT_SET)) {
      fprintf(stderr, "Error: missing %s %s: %s\n\n",
          ol[i].isopt ? "option" : "argument", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}



/* dump the current values */
__inline static void argopt_dump(const argopt_t *ao)
{
  int i, j, len = 2;
  opt_t *ol = ao->opts;

  /* get the width of the widest option */
  for (i = 0; i < ao->nopt; i++)
    if ( ol[i].sflag && (j = strlen(ol[i].sflag)) > len )
      len = j;

  /* print values of all options */
  for (i = 0; i < ao->nopt; i++) {
    const char *sflag = ol[i].sflag;

    if (sflag == NULL) sflag = "arg";
    fprintf(stderr, "%*s: ", len + 1, sflag);
    opt_fprintptr(stderr, ol + i);
    fprintf(stderr, ",  %s", ol[i].desc);
    if (ol[i].sarr != NULL) {
      fprintf(stderr, ", options: ");
      for ( j = 0; j < ol[i].scnt; j++ ) {
        fprintf(stderr, "%s", ol[i].sarr[j]);
        if ( j < ol[i].scnt - 1 ) fprintf(stderr, ", ");
      }
    }
    fprintf(stderr, "\n");
  }
}


#endif

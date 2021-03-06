#ifndef OPT_H__
#define OPT_H__



/* option either from arguments or configuration */
typedef struct {
  int isopt; /* is option (1) or argument (0) or cfg file entry (2) */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *key; /* key, for cfg files as in `key = val' */
  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  int ival; /* initial value, for switches */
  unsigned flags;
} opt_t;



#define OPT_MUST     0x0001  /* a mandatory argument or option */
#define OPT_SWITCH   0x0002  /* an option is a switch */
#define OPT_SET      0x0004  /* an argument/option is set */

/* translate string value in `o->val' into
 * actual ones through sscanf(), etc */
static int opt_getval(opt_t *o)
{
  const char *fmt = o->fmt;

  /* for a string argument, it can be obtained from NULL or "%s"
   * NULL: string memory from command-line
   * "%s":  string memory from sscpy()
   * if the string is read-only, like a output file name,
   * but the string is to be changed, use "%s" */
  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
    *((const char **) o->ptr) = o->val;
  } else if (strcmp(fmt, "%s") == 0) { /* copy the string */
    sscpy( *((char **) o->ptr), o->val);
  } else if (strcmp(fmt, "%b") == 0) { /* switch */
    /* switch the default value */
    if (o->flags & OPT_SET) return !o->ival;
    else return o->ival;
  } else { /* call sscanf */
    if (strcmp(fmt, "%r") == 0) /* real */
      fmt = (sizeof(real) == sizeof(float)) ? "%f" : "%lf";
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
static void opt_set(opt_t *o, const char *sflag, const char *key,
    const char *fmt, void *ptr, const char *desc)
{
  o->ch = '\0';
  if (key) { /* cfg file `key = val', not a command-line argument */
    o->isopt = 2;
  } else if (sflag) { /* option */
    o->isopt = 1;
    o->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    o->isopt = 0;
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
    fprintf(stderr, "unknown format (missing `%%') "
        "flag `%s\', fmt `%s', description: %s\n",
        sflag, fmt, desc);
  }
  if (strcmp(fmt, "%b") == 0) {
    fmt = "%d";
    o->flags |= OPT_SWITCH;
    o->ival = *((int *) ptr); /* save the initial value */
  }
  o->fmt = fmt;
  o->pfmt = NULL;
  o->desc = desc;
}



/* print the value of o->ptr */
#define opt_printptr(o) opt_fprintptr(stderr, o)
static void opt_fprintptr(FILE *fp, opt_t *o)
{
  const char *fmt;

  for (fmt = o->fmt; *fmt && *fmt != '%'; fmt++) ;
#define ELIF_PF_(fm, fmp, type) \
  else if (strcmp(fmt, fm) == 0) \
    fprintf(fp, (o->pfmt ? o->pfmt : fmp), *(type *)o->ptr)

  if (fmt == NULL || *fmt == '\0' || strcmp(fmt, "%s") == 0)
    fprintf(fp, "%s", (*(char **) o->ptr) ? (*(char **) o->ptr) : "NULL");
  ELIF_PF_("%b", "%d", int); /* switch */
  ELIF_PF_("%d", "%d", int);
  ELIF_PF_("%u", "%u", unsigned);
  ELIF_PF_("%x", "0x%x", unsigned);
  ELIF_PF_("%ld", "%ld", long);
  ELIF_PF_("%lo", "%lo", long);
  ELIF_PF_("%lu", "%lu", unsigned long);
  ELIF_PF_("%lx", "0x%lx", unsigned long);
#if 1  /* C99 or GCC extension */
  ELIF_PF_("%lld", "%lld", long long);
  ELIF_PF_("%llo", "%llo", long long);
  ELIF_PF_("%llu", "%llu", unsigned long long);
  ELIF_PF_("%llx", "0x%llx", unsigned long long);
#endif
  ELIF_PF_("%f", "%g", float);
  ELIF_PF_("%lf", "%g", double);
  ELIF_PF_("%Lf", "%Lg", long double);
  ELIF_PF_("%r", "%g", real);
  else fprintf(fp, "unknown %s-->%%d: %d", fmt, *(int *) o->ptr);
#undef ELIF_PF_
}



/* search an option list, return an option whose variable address is p */
static opt_t *opt_find(opt_t *ls, int n, const void *p)
{
  int i;

  for (i = 0; i < n; i++) if (ls[i].ptr == p) return ls + i;
  return NULL;
}



/* search an option list to see if an option is explicitly set */
static int opt_isset(opt_t *ls, int n, const void *p, const char *var)
{
  opt_t *o = opt_find(ls, n, p);
  if (o == NULL) {
    fprintf(stderr, "opt error: cannot find var %s, ptr %p\n", var, p);
    return 0;
  }
  return o->flags & OPT_SET ? 1 : 0;
}



#endif /* OPT_H__ */


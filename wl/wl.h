#ifndef WL_H__
#define WL_H__



#include <stdio.h>
#include <stdlib.h>
#include <float.h>





#define WL_FLATNESSABS 0x0010 /* use (hmax-hmin)/(hmax+hmin) for flatness */

typedef struct {
  int n0, n1, n;
  int isinvt; /* has entered the 1/t stage */
  unsigned flags;
  double *h; /* histogram */
  double *v; /* bias potential */
  double tot;
  double lnf; /* current updating factor lnf */
  double flatness; /* Wang-Landau threshold for the histogram flatness */
  double frac; /* Wang-Landau reduction factor for lnf */
  double c; /* c/t */
} wl_t;



__inline static wl_t *wl_open(int n0, int n1,
    double lnf0, double flatness, double frac,
    double c, unsigned flags)
{
  wl_t *wl;
  int i;

  if ( (wl = calloc(1, sizeof(*wl))) == NULL ) {
    fprintf(stderr, "no memory for WL\n");
    return NULL;
  }
  wl->n0 = n0;
  wl->n1 = n1;
  wl->n = n1 - n0;
  if ( (wl->h = calloc(1, sizeof(double) * wl->n)) == NULL ) {
    fprintf(stderr, "no memory for the WL histogram\n");
    free(wl);
    return NULL;
  }
  if ( (wl->v = calloc(1, sizeof(double) * wl->n)) == NULL ) {
    fprintf(stderr, "no memory for the WL potential\n");
    free(wl);
    return NULL;
  }
  for ( i = 0; i < wl->n; i++ ) {
    wl->h[i] = 0.0;
    wl->v[i] = 0.0;
  }
  wl->lnf = lnf0;
  wl->tot = 0;
  wl->isinvt = 0;
  wl->flatness = flatness;
  wl->frac = frac;
  wl->c = c;
  wl->flags = flags;
  return wl;
}



__inline static void wl_close(wl_t *wl)
{
  free(wl->h);
  free(wl->v);
  free(wl);
}



/* compute the total of the histogram */
__inline static double wl_gethtot(const double *h, int n)
{
  double htot = 0;
  int i;

  for ( i = 0; i < n; i++ ) {
    htot += h[i];
  }
  return htot;
}



/* clear the histogram */
__inline static void wl_clearh(double *h, int n)
{
  int i;

  for ( i = 0; i < n; i++ ) {
    h[i] = 0.0;
  }
}



/* trim the bottom of the potential */
__inline static void wl_trimv(double *v, int n)
{
  double vmin = v[0];
  int i;

  for ( i = 1; i < n; i++ ) {
    if ( v[i] < vmin ) {
      vmin = v[i];
    }
  }
  for ( i = 0; i < n; i++ ) {
    v[i] -= vmin;
  }
}



/* add an entry, update the histogram and potential */
__inline static int wl_add(wl_t *wl, int i)
{
  i -= wl->n0;
  if ( i < 0 || i >= wl->n ) {
    fprintf(stderr, "wl: out of range i %d, n %d\n", i, wl->n);
    return -1;
  }
  wl->h[i] += 1.0;
  wl->v[i] += wl->lnf;
  wl->tot += 1.0;
  return 0;
}



/* compute the histogram flatness from the absolute value */
__inline static double wl_getflatnessabs(const double *h, int n)
{
  double hmin, hmax;
  int i;

  hmin = hmax = h[0];
  for ( i = 1; i < n; i++ ) {
    if ( h[i] > hmax ) {
      hmax = h[i];
    } else if ( h[i] < hmin ) {
      hmin = h[i];
    }
  }
  return hmax > hmin ? (hmax - hmin) / (hmax + hmin) : 1.0;
}



/* compute the histogram flatness from the standard deviation */
__inline static double wl_getflatnessstd(const double *h, int n)
{
  double y, sh = 0, shh = 0;
  int i;

  for ( i = 0; i < n; i++ ) {
    y = h[i];
    sh += y;
    shh += y * y;
  }
  if ( sh <= 0 ) return 1.0;
  sh /= n;
  shh = shh / n - sh * sh;
  if ( shh < 0 ) shh = 0;
  return sqrt(shh) / sh;
}



__inline static double wl_getflatness(const wl_t *wl)
{
  if ( wl->flags & WL_FLATNESSABS ) {
    return wl_getflatnessabs(wl->h, wl->n);
  } else {
    return wl_getflatnessstd(wl->h, wl->n);
  }
}



/* lnf = 1/t */
__inline static double wl_lnfinvt(const wl_t *wl)
{
  return wl->c * wl->n / wl->tot;
}



/* update lnf, return 1 if the Wang-Landau stage is switched */
__inline static int wl_updatelnf(wl_t *wl)
{
  double flatness, nlnf, lnfinvt;

  if ( wl->isinvt ) {
    wl->lnf = wl_lnfinvt(wl);
    return 0;
  }

  flatness = wl_getflatness(wl);
  if ( flatness < wl->flatness ) {
    nlnf = wl->lnf * wl->frac;
    lnfinvt = wl_lnfinvt(wl);
    if ( nlnf < lnfinvt ) {
      fprintf(stderr, "changing lnf from %g to %g(1/t), flatness %g%%\n",
        wl->lnf, lnfinvt, flatness*100);
      wl->isinvt = 1;
      wl->lnf = lnfinvt;
    } else {
      fprintf(stderr, "changing lnf from %g to %g (1/t %g), flatness %g%%\n",
        wl->lnf, nlnf, lnfinvt, flatness*100);
      wl->lnf = nlnf;
      wl_clearh(wl->h, wl->n);
    }
    return 1;
  }
  return 0;
}



/* save data to file */
__inline static int wl_save(wl_t *wl, const char *fn)
{
  FILE *fp;
  int i;
  double htot;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  htot = wl_gethtot(wl->h, wl->n);
  wl_trimv(wl->v, wl->n);
  fprintf(fp, "# %d %d %g\n", wl->n0, wl->n1, wl->tot);
  for ( i = 0; i < wl->n; i++ ) {
    fprintf(fp, "%d %g %g %g\n",
        i + wl->n0, wl->v[i], wl->h[i]/htot, wl->h[i]);
  }
  fclose(fp);
  return 0;
}



/* load data from file */
__inline static int wl_load(wl_t *wl, const char *fn)
{
  FILE *fp;
  int i, i1, n, n0, n1;
  double x, y, v, tot;
  char ln[64000] = "";

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  if ( fgets(ln, sizeof ln, fp) == NULL
    || ln[0] != '#'
    || sscanf(ln + 1, "%d%d%lf", &n0, &n1, &tot) != 3 ) {
    fprintf(fp, "%s: bad information line!\n%s", fn, ln);
    fclose(fp);
    return -1;
  }
  if ( n0 != wl->n0 || n1 != wl->n1 ) {
    fprintf(stderr, "%s: dimensions mismatch [%d, %d) vs [%d, %d)\n",
        fn, n0, wl->n0, n1, wl->n1);
    fclose(fp);
    return -1;
  }
  wl->tot = tot;
  n = n1 - n0;

  for ( i = 0; i < n; i++ ) {
    if ( fgets(ln, sizeof ln, fp) == NULL ) {
      fprintf(stderr, "%s: cannot read for n %d\n", fn, i);
      fclose(fp);
      return -1;
    }
    if ( 4 != sscanf(ln, "%d %lf %lf %lf", &i1, &v, &x, &y)
      || i + n0 != i1 ) {
      fprintf(stderr, "%s: bad line %d\n%s", fn, i, ln);
      fclose(fp);
      return -1;
    }
    wl->v[i] = v;
    wl->h[i] = y;
  }
  fclose(fp);
  return 0;
}






#endif /* WL_H__ */


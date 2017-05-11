#ifndef POTTS2_H__
#define POTTS2_H__



/* two-dimensional Potts model */



#include "mtrand.h"



typedef struct {
  int l, n;
  int q; /* number of states */
  int E;
  int *s; /* 0 ... q - 1 */
  unsigned uproba[5]; /* temporary probability for MC transitions */
  /* for Wolff's algorithm */
  int *queue;
  int *used;
} potts2_t;



/* initialize an lxl Potts model */
__inline static potts2_t *potts2_open(int l, int q)
{
  int i, n;
  potts2_t *pt;

  if ( (pt = calloc(1, sizeof(*pt))) == NULL ) {
    fprintf(stderr, "no memory for is2\n");
    return NULL;
  }
  pt->l = l;
  pt->n = n = l * l;
  pt->q = q;
  if ( (pt->s = calloc(n, sizeof(*pt->s))) == NULL ) {
    fprintf(stderr, "no memory for pt->n\n");
    return NULL;
  }
  for (i = 0; i < n; i++) pt->s[i] = 0;
  pt->E = -2*n;
  pt->uproba[0] = 0xffffffff;
  if ( (pt->queue = calloc(n, sizeof(*pt->queue))) == NULL ) {
    fprintf(stderr, "no memory for pt->queue\n");
    return NULL;
  }
  if ( (pt->used = calloc(n, sizeof(*pt->used))) == NULL ) {
    fprintf(stderr, "no memory for pt->used\n");
    return NULL;
  }
  return pt;
}



__inline static void potts2_close(potts2_t *pt)
{
  free(pt->s);
  free(pt->queue);
  free(pt->used);
  free(pt);
}



/* set transition probability */
__inline static void potts2_setuproba(double bet, unsigned *p)
{
  double x = exp(-bet), y = x;
  p[1] = (unsigned) ((double)(0xffffffff) * y); y *= x;
  p[2] = (unsigned) ((double)(0xffffffff) * y); y *= x;
  p[3] = (unsigned) ((double)(0xffffffff) * y); y *= x;
  p[4] = (unsigned) ((double)(0xffffffff) * y);
}



/* pick a random site, count neighbors with different spins */
__inline static int potts2_pick(const potts2_t *pt, int *h, int *sn)
{
  int ix, ixp, ixm, iy, iyp, iym, id, s, so;
  int l = pt->l, n = pt->n;

  id = (int) ( rand01() * n );
  ix = id % l;
  iy = id - ix;
  ixp = ( ix + 1 ) % l;
  ixm = ( ix + l - 1 ) % l;
  iyp = ( iy + l ) % n;
  iym = ( iy + n - l ) % n;
  so = pt->s[id];
  *sn = (so + 1 + (int)((pt->q - 1) * rand01())) % pt->q;
  *h = 0;
  *h += ((s = pt->s[iy  + ixp]) == so); *h -= (s == *sn);
  *h += ((s = pt->s[iy  + ixm]) == so); *h -= (s == *sn);
  *h += ((s = pt->s[iyp +  ix]) == so); *h -= (s == *sn);
  *h += ((s = pt->s[iym +  ix]) == so); *h -= (s == *sn);
  return id;
}



/* flip site id, with (-h) pt the energy before the flip */
__inline static int potts2_flip(potts2_t *pt, int id, int h, int sn)
{
  pt->s[id] = sn;
  return pt->E += h;
}


/* faster macros for systems with fixed (upon compiling) size
 * to use them one must define POTTS2_LB before including
 * POTTS2_PICK()/POTTS2_PSEQ() and POTTS2_FLIP() */
#ifdef  POTTS2_LB  /* L = 2^LB, N = L*L */
#define POTTS2_L   (1 << POTTS2_LB)
#define POTTS2_N   (POTTS2_L * POTTS2_L)

#define POTTS2_GETH(pt, id, h, sn) { \
  unsigned ix, ixp, ixm, iy, iyp, iym; \
  int s, so; \
  ix = id % POTTS2_L; \
  iy = id - ix; \
  ixp = (ix + 1) % POTTS2_L; \
  ixm = (ix + (POTTS2_L - 1)) % POTTS2_L; \
  iyp = (iy + POTTS2_L) % POTTS2_N; \
  iym = (iy + (POTTS2_N - POTTS2_L)) % POTTS2_N; \
  so = pt->s[id]; \
  sn = (so + 1 + (int) ((pt->q - 1) * rand01())) % pt->q; \
  h  = ((s = pt->s[iy  + ixp]) == so); h -= (s == sn); \
  h += ((s = pt->s[iy  + ixm]) == so); h -= (s == sn); \
  h += ((s = pt->s[iyp +  ix]) == so); h -= (s == sn); \
  h += ((s = pt->s[iym +  ix]) == so); h -= (s == sn); }
#define POTTS2_IRND(pt, id)  id = mtrand() >> (32 - 2*POTTS2_LB);
/* random picking */
#define POTTS2_PICK(pt, id, h, sn) { POTTS2_IRND(pt, id); POTTS2_GETH(pt, id, h, sn); }
#define POTTS2_ISEQ(pt, id)  id = (id + 1) % POTTS2_N;
/* sequential picking */
#define POTTS2_PSEQ(pt, id, h, sn) { POTTS2_ISEQ(pt, id); POTTS2_GETH(pt, id, h, sn); }

#define POTTS2_FLIP(pt, id, h, sn) { \
  pt->E += h; pt->s[id] = sn; }

#else

#define POTTS2_PICK(pt, id, h, sn)  id = potts2_pick(pt, &h, &sn)
#define POTTS2_FLIP(pt, id, h, sn)  potts2_flip(pt, id, h, sn)

#endif



/* compute total energy */
__inline static int potts2_energy(potts2_t *pt)
{
  int l, n, i, j, e;

  e = 0;
  l = pt->l;
  n = l * l;
  for ( i = 0; i < n; i += l ) {
    for ( j = 0; j < l; j++ ) {
      int id = i + j;
      int idr = i + (j + 1) % l;
      int idu = (i + l) % n + j;
      int s = pt->s[id];
      int su = pt->s[idu];
      int sr = pt->s[idr];
      e += (s == su) + (s == sr);
    }
  }
  return pt->E = -e;
}



/* add spin j to the queue if s[j] pt different from s
 * return the spin */
__inline static int potts2_addtoqueue(potts2_t *pt, int j,
    int so, int sn, double r, int *cnt)
{
  int sj = pt->s[j];

  if ( sj == so && !pt->used[j] && rand01() < r ) {
    pt->queue[ (*cnt)++ ] = j;
    pt->used[j] = (char) 1;
  }
  return (sj == so) - (sj == sn);
}



/* Wolff algorithm */
__inline static int potts2_wolff(potts2_t *pt, double padd)
{
  int l = pt->l, n = pt->n, i, ix, iy, id, so, sn, cnt = 0, h = 0;

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  so = pt->s[id];
  sn = (so + 1 + (int) (rand01() * pt->q)) % pt->q;
  pt->queue[ cnt++ ] = id;
  for ( i = 0; i < n; i++ )
    pt->used[i] = 0;
  pt->used[id] = (char) 1;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = pt->queue[i];
    pt->s[id] = sn;
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    h += potts2_addtoqueue(pt, iy + (ix + 1) % l,     so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, iy + (ix + l - 1) % l, so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, (iy + l) % n + ix,     so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, (iy + n - l) % n + ix, so, sn, padd, &cnt);
  }

  pt->E += h;
  return 0;
}



__inline static int potts2_save(const potts2_t *pt, const char *fname)
{
  FILE *fp;
  int i, j, l, *p;

  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fname);
    return -1;
  }
  l = pt->l;
  fprintf(fp, "2 %d %d %d\n", l, l, pt->n);
  for (p = pt->s, i = 0; i < l; i++) {
    for (j = 0; j < l; j++, p++)
      fprintf(fp, "%c", (*p + '0'));
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



__inline static int potts2_load(potts2_t *pt, const char *fname)
{
  FILE *fp;
  int i, lx, ly, n, c;
  char s[80];

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fname);
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "missing first line %s\n", fname);
    fclose(fp);
    return -1;
  }
  if (4 != sscanf(s, "%d%d%d%d", &i, &lx, &ly, &n)
      || i != 2 || lx != ly || lx != pt->l || n != pt->n) {
    fprintf(stderr, "bad setting: %dD, %dx%d = %d\n", i, lx, ly, n);
    return -1;
  }
  for (i = 0; i < n; i++) {
    while ((c = fgetc(fp)) != EOF && c == '\n') ;
    if (c == EOF) break;
    pt->s[i] = (int) (c - '0');
  }
  if (i < n)
    fprintf(stderr, "%s: data stopped at i = %d\n", fname, i);
  fclose(fp);
  potts2_energy(pt);
  return 0;
}



#endif /* POTTS2_H__ */


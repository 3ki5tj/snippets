#include "mtrand.h"
#include "ave.h"
#include "argopt.h"
#include <time.h>

int q = 10000;
int nsys = 1000;
double mean = 1.0; /* mean of the relative probability */
double width = 0.0; /* width of the relative probability */
long nsteps = 10000;
int npart = 2;
int npart2 = 4;
double gam = 1.0;
long nstrep = 100;
char *fnlog = "walk.log";
int randomize = 0; // random seed from the current time


static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "Random walk";
  argopt_add(ao, "-q", "%d",  &q,         "number of states");
  argopt_add(ao, "-M", "%d",  &nsys,      "number of replica systems (walkers)");
  argopt_add(ao, "-m", "%lf", &mean,      "mean of the relative probability");
  argopt_add(ao, "-w", "%lf", &width,     "width of the relative probability");
  argopt_add(ao, "-t", "%ld", &nsteps,    "number of steps");
  argopt_add(ao, "-r", "%ld", &nstrep,    "number of steps to report");
  argopt_add(ao, "-P", "%d",  &npart,     "number of partitions for the block method");
  argopt_add(ao, "-Q", "%d",  &npart2,    "number of partitions for the block method (2)");
  argopt_add(ao, "-g", "%lf", &gam,       "gamma factor for the exponential extrapolation");
  argopt_add(ao, "-R", "%b",  &randomize, "using the current time as the seed of the random number generator");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  if ( npart2 == npart ) npart2 = npart * 2;
  argopt_dump(ao);
  argopt_close(ao);
}

typedef struct {
  int q;
  double mean;
  double width;
  double entref;
  double *proba, *Proba;
} distr_t;

typedef struct {
  int s;
  int q;
  int trjn, trji, *trj;
  long tot, *cnt; /* histogram data */
  double ent, entc, entcb, ents, entsb;
  double mest;
} walker_t;

static distr_t *distr_open(int q, double mean, double width)
{
  distr_t *d;
  int s;
  double x;

  xnew(d, 1);
  d->q = q;
  d->mean = mean;
  d->width = width;

  /* initialize the probabilities of the states */
  xnew(d->proba, q);
  for ( s = 0; s < q; s++ ) {
    d->proba[s] = mean - width + width * 2 * rand01();
  }
  /* compute the cumulative probabilities */
  xnew(d->Proba, q + 1);
  d->Proba[0] = 0;
  for ( s = 0; s < q; s++ ) {
    d->Proba[s+1] = d->Proba[s] + d->proba[s];
    printf("%4d: %8.3f %8.3f\n", s, d->proba[s], d->Proba[s+1]);
  }
  x = d->Proba[q];
  for ( s = 0; s < q; s++ ) {
    d->proba[s] /= x;
    d->Proba[s+1] /= x;
  }
  /* compute the entropy */
  d->entref = 0;
  for ( s = 0; s < q; s++ ) {
    x = d->proba[s];
    d->entref += -log(x) * x;
  }
  printf("%d states, entropy %g\n", q, d->entref);
  //getchar();
  return d;
}

static int distr_sample(distr_t *d)
{
  int q = d->q;
  double x;

  if ( d->width <= 0 ) {
    return (int) (rand01() * q);
  } else {
    x = rand01();
    /* binary search to find the brackets that contain x */
    int l = 0, r = q, m;
    while ( l + 1 < r ) {
      m = (l + r)/2;
      //printf("l %d, r %d, m %d, x %g\n", l, r, m, x);
      if ( d->Proba[m] > x ) {
        r = m;
      } else {
        l = m;
      }
    }
    //getchar();
    return l;
  }
}


static void distr_close(distr_t *d)
{
  free(d->proba);
  free(d->Proba);
  free(d);
}


static walker_t *walker_open(int q, int trjlen)
{
  walker_t *w;
  int s;

  xnew(w, 1);
  w->q = q;
  w->s = 0;

  w->tot = 0;
  xnew(w->cnt, q);
  for ( s = 0; s < q; s++ )
    w->cnt[s] = 0;

  /* initialize the trajectory data */
  w->trjn = trjlen;
  w->trji = 0;
  xnew(w->trj, w->trjn);

  return w;
}

/* count occurrences from frame start to frame end */
static void walker_count(walker_t *w, int start, int end)
{
  int i, s;

  for ( s = 0; s < w->q; s++ )
    w->cnt[s] = 0;
  for ( i = start; i < end; i++ ) {
    w->cnt[ w->trj[i] ] += 1;
  }
}

/* compute the entropy from a segment [start, end) */
static double walker_ent_seg(walker_t *w, int start, int end)
{
  double ent = 0, tot = (double) (end - start), pr; // ent1 = log(w->q);
  int s;

  walker_count(w, start, end);
  for ( s = 0; s < w->q; s++ ) {
    long c = w->cnt[s];
    if ( c <= 0 ) continue;
    pr = 1.0*c/tot;
    ent += -pr*log(pr);
  }
  return ent;
}

static double walker_entropy(walker_t *w, int npart, int npart2, int verbose)
{
  int trjn = w->trji, ip, blksz, blksz2;
  double entp = 0, entp2 = 0, entt, entc2;

  /* entropy estimate from the entire trajectory */
  entt = walker_ent_seg(w, 0, trjn);

  /* entropy estimated from trajectory block(s) */
  blksz = trjn / npart;
  for ( ip = 0; ip < npart; ip++ ) {
    entp += walker_ent_seg(w, ip * blksz, (ip + 1) * blksz);
  }
  entp /= npart;

  blksz2 = trjn / npart2;
  for ( ip = 0; ip < npart2; ip++ ) {
    entp2 += walker_ent_seg(w, ip * blksz2, (ip + 1) * blksz2);
  }
  entp2 /= npart2;

  /* since we have St = S - a/t, Sp = S - a/blksz
   * S = (St*t - Sp*blksz)/(t-blksz); */
  w->ent = entt;
  w->entc = (entt*trjn - entp*blksz)/(trjn - blksz);
  {
    double x = (entt - entp)/(trjn - blksz);
    w->entcb = w->entc - trjn*blksz/3*x*x;
  }

  /* fitting to log(1+gam*c/t)/gam */
  {
    double ds = entt - entp, xp = exp(gam*ds);
    double ds2 = entt - entp2, xp2 = exp(gam*ds2);
    double xpc, xpc2, xpcb;

    {
      if ( xp >= .99*trjn/blksz ) xp = .99*trjn/blksz;
      xpc = (trjn - xp*blksz)/(trjn - blksz);
      w->ents = entt - log(xpc)/gam;
      w->mest = 2*trjn*(xp - 1)/(1.*trjn/blksz - xp);

      if ( xp2 >= .99*trjn/blksz2 ) xp2 = .99*trjn/blksz2;
      xpc2 = (trjn - xp2*blksz2)/(trjn - blksz2);
      //ents2 = entt - log(xpc2);

      xpcb = (xpc*blksz - xpc2*blksz2)/(blksz - blksz2);
      if ( xpcb < 0 ) xpcb = xpc;
      w->entsb = entt - log(xpcb)/gam;

      w->entsb = 2*w->entc - w->ents;
    }
  }

  return w->ent;
}

static void walker_close(walker_t *w)
{
  free(w->cnt);
  free(w->trj);
  free(w);
}

static void work(void)
{
  walker_t **w;
  av_t avent[1], aventc[1], aventcb[1], avents[1], aventsb[1], avmest[1];
  double ent, entc, entcb, ents, entsb, mest;
  double var, varc, varcb, vars, varsb, mvar;
  distr_t *distr = NULL;
  int i;
  long t;
  FILE *fplog;

  if ( (fplog = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open log file %s\n", fnlog);
    fplog = stderr;
  }

  distr = distr_open(q, mean, width);

  xnew(w, nsys);
  for ( i = 0; i < nsys; i++ ) {
    w[i] = walker_open(q, nsteps);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      w[i]->s = distr_sample(distr);
      w[i]->trj[ w[i]->trji++ ] = w[i]->s;
    }

    if ( t % nstrep == 0 ) {
      av_clear(avent);
      av_clear(aventc);
      av_clear(aventcb);
      av_clear(avents);
      av_clear(aventsb);
      av_clear(avmest);
      for ( i = 0; i < nsys; i++ ) {
        walker_entropy(w[i], npart, npart2, i <= 0);
        av_add(avent,  w[i]->ent);
        av_add(aventc, w[i]->entc);
        av_add(aventcb, w[i]->entcb);
        av_add(avents, w[i]->ents);
        av_add(aventsb, w[i]->entsb);
        av_add(avmest, w[i]->mest);
      }
      ent  = av_getave(avent,  &var);
      entc = av_getave(aventc, &varc);
      entcb = av_getave(aventcb, &varcb);
      ents = av_getave(avents, &vars);
      entsb = av_getave(aventsb, &varsb);
      mest = av_getave(avmest, &mvar);
      printf("%9ld: entropy %8.4f(%.4f), %8.4f, %8.4f; %8.4f, %8.4f(%8.4f), m %6.0f(%5.0f)\n",
          t, ent, sqrt(var), entc, entcb, ents, entsb, distr->entref, mest, sqrt(mvar));
      fprintf(fplog, "%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t,
          ent, sqrt(var), distr->entref,
          entc, sqrt(varc), entcb, sqrt(varcb),
          ents, sqrt(vars), entsb, sqrt(varsb),
          mest, sqrt(mest));
      fflush(fplog);
    }
  }

  distr_close(distr);
  for ( i = 0; i < nsys; i++ )
    walker_close(w[i]);
  free(w);
  fclose(fplog);
}

int main(int argc, char **argv)
{
  doargs(argc, argv);
  if ( randomize ) {
    mtscramble(time(NULL) + clock());
  }
  work();
  return 0;
}

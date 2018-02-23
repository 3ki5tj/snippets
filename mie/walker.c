#include "mtrand.h"
#include "ave.h"
#include "argopt.h"
#include <time.h>

int q = 10000;
int nsys = 1000;
long nsteps = 10000;
int npart = 2;
long nstrep = 100;
char *fnlog = "walk.log";

static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "Random walk";
  argopt_add(ao, "-q", "%d", &q, "number of states");
  argopt_add(ao, "-M", "%d", &nsys, "number of replica systems (walkers)");
  argopt_add(ao, "-t", "%ld", &nsteps, "number of steps");
  argopt_add(ao, "-r", "%ld", &nstrep, "number of steps to report");
  argopt_add(ao, "-P", "%d", &npart, "number of partitions for the block method");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

typedef struct {
  int s;
  int q;
  int trjn, trji, *trj;
  long tot, *cnt; /* histogram data */
  double ent, entc, ents;
  double mest;
} walker_t;

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

static double walker_ent_seg(walker_t *w, int start, int end)
{
  double ent = 0, tot = (double) (end - start), pr;
  int s;

  walker_count(w, start, end);
  for ( s = 0; s < w->q; s++ ) {
    long c = w->cnt[s];
    if ( c <= 0 ) continue;
    pr = 1.0*c/tot;
    ent += -pr*log(pr);
    //printf("state %4d: %ld, tot %g, ent %g\n", s, c, tot, ent);
  }
  //printf("ent %g\n", ent); getchar();
  return ent;
}

static double fc(double m, double t)
{
  return (t/m + 1) * log(1 + m/t) - 1 - 0.5/t;
}

static double walker_entropy(walker_t *w, int npart)
{
  int trjn = w->trji, ip, blksz;
  double entp = 0, entt;

  /* entropy estimate from the entire trajectory */
  entt = walker_ent_seg(w, 0, trjn);

  /* entropy estimated from trajectory block(s) */
  blksz = trjn / npart;
  for ( ip = 0; ip < npart; ip++ ) {
    entp += walker_ent_seg(w, ip * blksz, (ip + 1) * blksz);
  }
  entp /= npart;

  /* since we have St = S - a/t, Sp = S - a*npart/t
   * a/t = (St - Sp)/(npart-1)
   * S = (npart*St - Sp)/(npart-1); */
  w->ent = entt;
  w->entc = (entt*npart - entp)/(npart - 1);

  int tr;
  double del, del0 = entt - entp;
  w->mest = (entt - entp)*blksz*npart*2/(npart - 1) + 1;

  for ( tr = 1; tr <= 100; tr++ ) {
    del = fc(w->mest, blksz) - fc(w->mest, blksz*npart);
    printf("%5d: m %g, del = %g, %g, corr %g\n", tr, w->mest, del0, del, fc(w->mest, trjn));
    w->mest *= del0/del;
    if ( fabs(del0 - del) < 1e-5 ) break;
  }
  getchar();

  w->ents = entt + fc(w->mest, trjn);
  //w->ents = entt + 0.5*(w->mest - 1)/trjn;
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
  av_t avent[1], aventc[1], avents[1];
  double ent, entc, ents;
  double var, varc, vars;
  int i;
  long t;
  FILE *fplog;

  if ( (fplog = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open log file %s\n", fnlog);
    fplog = stderr;
  }

  xnew(w, nsys);
  for ( i = 0; i < nsys; i++ ) {
    w[i] = walker_open(q, nsteps);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      w[i]->s = (int) (rand01() * q);
      w[i]->trj[ w[i]->trji++ ] = w[i]->s;
    }

    if ( t % nstrep == 0 ) {
      av_clear(avent);
      av_clear(aventc);
      av_clear(avents);
      for ( i = 0; i < nsys; i++ ) {
        walker_entropy(w[i], npart);
        av_add(avent,  w[i]->ent);
        av_add(aventc, w[i]->entc);
        av_add(avents, w[i]->ents);
      }
      ent  = av_getave(avent,  &var);
      entc = av_getave(aventc, &varc);
      ents = av_getave(avents, &vars);
      printf("%9ld: entropy %8.4f, %8.4f, %8.4f(%8.4f), %g\n",
          t, ent, entc, ents, log(q), w[0]->mest);
      fprintf(fplog, "%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t,
          ent, sqrt(var), entc, sqrt(varc), log(q), ents, sqrt(vars));
    }
  }

  for ( i = 0; i < nsys; i++ )
    walker_close(w[i]);
  free(w);
  fclose(fplog);
}

int main(int argc, char **argv)
{
  doargs(argc, argv);
  mtscramble(time(NULL) + clock());
  work();
  return 0;
}

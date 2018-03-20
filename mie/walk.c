#include "mtrand.h"
#include "ave.h"
#include "argopt.h"
#include <time.h>

int q = 10000;
int nsys = 1000;
long nsteps = 10000;
int npart = 2;
int npart2 = 4;
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
  argopt_add(ao, "-Q", "%d", &npart2, "number of partitions for the block method (2)");
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
  double ent, entc, entcb, ents, entsb;
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

/*
// for brute-force fitting, deprecated
static double ftail(double t, double m)
{
  // goodness
  return log(1+0.5*m/t);  // goodness 20
  //return (1+t/m)*log(1+m/t) - 1; // func1, goodness 10
  //return 2*log(1+0.25*m/t);  // goodness 5
  //return 0.5*m/t; // goodness 0;
}
*/

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
  entc2 = (entt*trjn - entp2*blksz2)/(trjn - blksz2);
  w->entcb = (w->entc*blksz - entc2*blksz2)/(blksz - blksz2);

  /* fitting to log(1+m/2/t) */
  {
    double ds = entt - entp, xp = exp(ds);
    double ds2 = entt - entp2, xp2 = exp(ds2);
    double xpc, xpc2, xpcb;

    //double alpha = 0.5, threshold = log((1 + alpha*trjn/blksz)/(1 + alpha));
    //printf("ds %g vs %g\n", ds, threshold); getchar();
    //if ( ds < threshold ) {
    //  w->ents = w->entc;
    //  w->entsb = w->entcb;
    //} else
    {
      if ( xp >= .99*trjn/blksz ) xp = .99*trjn/blksz;
      xpc = (trjn - xp*blksz)/(trjn - blksz);
      w->ents = entt - log(xpc);
      w->mest = 2*trjn*(xp - 1)/(1.*trjn/blksz - xp);

      if ( xp2 >= .99*trjn/blksz2 ) xp2 = .99*trjn/blksz2;
      xpc2 = (trjn - xp2*blksz2)/(trjn - blksz2);
      //ents2 = entt - log(xpc2);

      xpcb = (xpc*blksz - xpc2*blksz2)/(blksz - blksz2);
      if ( xpcb < 0 ) xpcb = xpc;
      w->entsb = entt - log(xpcb);
      //printf("m %g, ds %g,%g xp %g,%g xpc %g,%g S %g -> %g, %g\n", w->mest, ds, ds2, xp, xp2, xpc, xpc2, entt, w->ents, w->entsb);
      //getchar();
    }
  }

#if 0
  { // working out the polynomial interpolation
    double a2 = (entt - entp)/(1.0/blksz - 1.0/trjn);
    double S2 = entt + a2/trjn;
    double a3 = (entt - entp2)/(1.0/blksz2 - 1.0/trjn);
    double S3 = entt + a3/trjn;
    double S = (S2*blksz - S3*blksz2)/(blksz - blksz2);
    double b = (S2 - S3)*trjn*blksz*blksz2/(blksz - blksz2);
    double a = (S - entt)*trjn + b/trjn;
    printf("2: a %g, S %g,%g, entt %g, %g; entp %g, %g\n", a2, S2, w->entc, entt, S2-a2/trjn, entp, S2-a2/blksz);
    printf("3: a %g, S %g,%g, entt %g, %g; entp %g, %g\n", a3, S3, w->entcb, entt, S3-a3/trjn, entp2, S3-a3/blksz2);
    printf("S: %g/%g, a %g, b %g, entt %g, %g\n", S, w->ents, a, b, entt, S-a/trjn+b/trjn/trjn);
  }
#endif

#if 0 /* fitting to S = S0*(1 + a/t)/(1 + c/t) */
  {
    double k2 = (entt - entp)/(trjn - blksz);
    double k3 = (entt - entp2)/(trjn - blksz2);
    double c = (w->entc - w->entcb)/(k3 - k2);
    double S = w->entc + c*k2;
    double a = (S - entt)*(trjn + c);
    w->entsb = S;
    w->mest = 2 * (a - c);
    //printf("S %g, a %g, c %g, entt %g, %g; entp %g, %g\n", S, a, c, entt, S-a/trjn/(1+c/trjn), entp, S-a/blksz/(1+c/blksz));
  }
#endif

#if 0 /* brute-force extrapolation */
  {
    double del, delr = entt - entp, del2, del2r = entt - entp2;
    double tailt = 0, tailp, tailp2;
    int i;

    w->mest = (entt - entp)/(1.0/blksz - 1.0/trjn);
    for ( i = 0; i < 1000; i++ ) {
      tailt = ftail(trjn, w->mest);
      tailp = ftail(blksz, w->mest);
      tailp2 = ftail(blksz2, w->mest);
      del = tailp - tailt;
      del2 = tailp2 - tailt;
      if (verbose ) {
        printf("%d: t %d, S %g, %g/%g, S %g(%g) mest %g, del %g, %g\n",
          i, trjn, entt, entp, entp2, entt + tailt, log(q), w->mest, del/delr, del2/del2r);
      }
      w->mest *= sqrt(delr/del);
      if ( fabs(del/delr - 1) < 0.001 ) break;
      //w->mest *= delr/del*del2r/del2;
      //if ( fabs(del/delr*del2/del2r - 1) < 0.01 ) break;
    }
    w->entsb = entt + tailt;
    if ( verbose ) {
      getchar();
    }
  }
#endif

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
      printf("%9ld: entropy %8.4f(%g), %8.4f, %8.4f; %8.4f, %8.4f(%8.4f), m %6.0f(%5.0f)\n",
          t, ent, var, entc, entcb, ents, entsb, log(q), mest, sqrt(mvar));
      fprintf(fplog, "%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t,
          ent, sqrt(var), log(q),
          entc, sqrt(varc), entcb, sqrt(varcb),
          ents, sqrt(vars), entsb, sqrt(varsb),
          mest, sqrt(mest));
      fflush(fplog);
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
  //mtscramble(time(NULL) + clock());
  work();
  return 0;
}

#ifndef AVE_H__
#define AVE_H__

/* simple version */
typedef struct {
  double cnt, sum, sqr;
} av_t;

__inline static void av_clear(av_t *a)
{
  a->cnt = 0;
  a->sum = 0;
  a->sqr = 0;
}

__inline static void av_add(av_t *a, double x)
{
  a->cnt += 1;
  a->sum += x;
  a->sqr += x * x;
}

__inline static double av_getave(av_t *a, double *var)
{
  double ave = (a->cnt > 0) ? a->sum / a->cnt : 0;
  if ( var != NULL )
    *var = (a->cnt > 0) ? a->sqr / a->cnt - ave * ave : 0;
  return ave;
}


/* more careful version */
typedef struct {
  double cnt, sum, sqr, ave, var;
} ave_t;

__inline static void ave_clear(ave_t *a)
{
  a->cnt = 0;
  a->sum = 0;
  a->sqr = 0;
  a->ave = 0;
  a->var = 0;
}

__inline static void ave_add(ave_t *a, double x)
{
  double dx;
  a->sum += x;
  if ( a->cnt > 0 ) { /* update the variance accumulator */
    dx = x - a->ave;
    a->sqr += dx * dx * a->cnt / (a->cnt + 1);
  }
  a->cnt += 1;
  a->ave = a->sum / a->cnt;
  a->var = a->sqr / a->cnt;
}


#endif /* AVE_H__ */

#define ZCOM_PICK
#define ZCOM_MD
#include "zcom.h"



#define INFTIME  ((real) 1e9)



typedef struct {
  int n;
  int d;
  int dof;
  real time;
  real *box; /* box dimensions, box[0] ... box[d-1] */
  real *x, *v;
  real *xij, *vij;
  real *coltime; /* collision time */
  int *partner; /* collide with whom */
} hsmd_t;



static void hsmd_close(hsmd_t *m)
{
  free(m->box);
  free(m->x);
  free(m->v);
  free(m->xij);
  free(m->vij);
  free(m->coltime);
  free(m->partner);
  free(m);
}



static int hsmd_initpos3d(hsmd_t *m)
{
  const double stretch = 1.05;
  int d = m->d;
  int *nm, ntot;
  int ix, iy, iz, i, j;

  die_if ( d != 3, "currently only working on three dimensions d %d\n", d);
  xnew(nm, d);
  ntot = 1;
  for ( j = 0; j < d; j++ ) {
    nm[j] = (int) (m->box[j] / stretch);
    ntot *= nm[j];
  }
  die_if (ntot < m->n, "cannot afford %d particles, max %d\n", m->n, ntot);

  i = 0;
  for ( ix = 0; ix < nm[0]; ix++ )
    for ( iy = 0; iy < nm[1]; iy++ )
      for ( iz = 0; iz < nm[2]; iz++ ) {
        m->x[i*d + 0] = ix * stretch;
        m->x[i*d + 1] = iy * stretch;
        m->x[i*d + 2] = iz * stretch;
        if ( ++i >= m->n ) goto END;
      }
END:
  free(nm);
  return 0;
}



/* apply the minimal image convention */
__inline static real wrap(real x, real box)
{
  return x > box*.5 ? x - box : x < -box*.5 ? x + box : x;
}



/* apply the minimal image convention */
__inline static real wrap1(real x, real box)
{
  return x - ((int) (1.0*x/box + 1000) - 1000) * box;
}



/* compute the collision time of i and j */
static real hsmd_coltimeij(hsmd_t *m, int i, int j,
    real *coltime, int *partner)
{
  int d = m->d, k;
  real bij, xij2, vij2, discr, tij = INFTIME;

  bij = 0;
  for ( k = 0; k < d; k++ ) {
    m->xij[k] = wrap(m->x[i*d + k] - m->x[j*d + k], m->box[k]);
    m->vij[k] = m->v[i*d + k] - m->v[j*d + k];
    bij += m->xij[k] * m->vij[k];
  }

  if ( bij < 0 ) {
    xij2 = vij2 = 0;
    for ( k = 0; k < d; k++ ) {
      xij2 += m->xij[k] * m->xij[k];
      vij2 += m->vij[k] * m->vij[k];
    }
    discr = bij * bij - vij2 * (xij2 - 1);
    if ( discr > 0 ) {
      tij = -(bij + (real) sqrt(discr)) / vij2;
      //printf("i %d, j %d, bij %g, discr %g, tij %g\n", i, j, bij, discr, tij);
      if ( tij < coltime[i] ) {
        coltime[i] = tij;
        partner[i] = j;
      }
      if ( tij < coltime[j] ) {
        coltime[j] = tij;
        partner[j] = i;
      }
    }
  }
  return tij;
}



/* initialize the collision list */
static void hsmd_initcoltime(hsmd_t *m)
{
  int n = m->n, i, j;

  for ( i = 0; i < n; i++ ) m->coltime[i] = INFTIME;

  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ )
      hsmd_coltimeij(m, i, j, m->coltime, m->partner);
}



/* verify the collision list */
__inline static void hsmd_verifycoltime(hsmd_t *m)
{
  int n = m->n, i, j;
  real *coltime;
  int *partner;

  xnew(coltime, n);
  xnew(partner, n);
  for ( i = 0; i < n; i++ ) coltime[i] = INFTIME;
  for ( i = 0; i < n - 1; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      hsmd_coltimeij(m, i, j, coltime, partner);
    }
  }
  for ( i = 0; i < n; i++ )
    if ( (partner[i] != m->partner[i] || fabs(coltime[i] - m->coltime[i]) > 1e-5) && coltime[i] < INFTIME*.5) {
      printf("i %d: %g %g %d %d\n", i, coltime[i], m->coltime[i], partner[i], m->partner[i]);
      getchar();
    }
  free(coltime);
  free(partner);
}



static void hsmd_vscale(hsmd_t *m, real tp)
{
  real ekin = md_ekin(m->v, m->n * m->d, m->dof, NULL);
  md_vscale(m->v, m->n * m->d, m->dof, tp, ekin, NULL, NULL);
}



static hsmd_t *hsmd_open(int n, int d, real *box)
{
  hsmd_t *m;
  int i, j;

  xnew(m, 1);
  m->n = n;
  m->d = d;
  m->dof = (d == 3) ? n * 3 - 6 : (d == 2) ? n * 2 - 3 : n * d - d;
  m->time = 0;

  xnew(m->box, d);
  for ( j = 0; j < d; j++ ) m->box[j] = box[j];

  xnew(m->xij, d);
  xnew(m->vij, d);

  xnew(m->x, n*d);
  hsmd_initpos3d(m);

  xnew(m->v, n*d);
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < d; j++ )
      m->v[i*d+j] = rnd0() - 0.5;

  md_shiftcom(m->v, m->n, m->d);
  md_shiftang(m->x, m->v, m->n, m->d);
  hsmd_vscale(m, 1); /* scale the velocity to unit temperature */

  xnew(m->coltime, n);
  xnew(m->partner, n);

  hsmd_initcoltime(m); /* initialize the collision time */
  return m;
}



/* find out the next collision pair */
static void hsmd_getij(const hsmd_t *m, int *i, int *j)
{
  int k;

  *i = 0;
  for ( k = 1; k < m->n; k++ )
    if ( m->coltime[k] < m->coltime[*i] )
      *i = k;
  *j = m->partner[*i];
}



/* move along the velocity by tij */
static void hsmd_move(hsmd_t *m, real tij)
{
  int n = m->n, d = m->d;
  int l, k;

  /* move forward for tij */
  for ( l = 0; l < n; l++ ) {
    m->coltime[l] -= tij;
    for ( k = 0; k < d; k++ ) {
      m->x[l*d + k] = wrap1(m->x[l*d + k] + m->v[l*d + k] * tij, m->box[k]);
    }
  }
}



/* change the velocity by collision */
static real hsmd_collide(hsmd_t *m, int i, int j)
{
  int d = m->d, k;
  real bij = 0, xij2 = 0;

  for ( k = 0; k < d; k++ ) {
    m->xij[k] = wrap(m->x[i*d + k] - m->x[j*d + k], m->box[k]);
    m->vij[k] = m->v[i*d + k] - m->v[j*d + k];
    bij += m->xij[k] * m->vij[k];
    xij2 += m->xij[k] * m->xij[k];
  }
  for ( k = 0; k < d; k++ ) {
    m->v[i*d + k] -= bij * m->xij[k];
    m->v[j*d + k] += bij * m->xij[k];
  }
  //printf("xij2 %g\n", xij2);
  return xij2;
}



/* update the collision time list
 * after a collision of i0 and j0 */
static void hsmd_update(hsmd_t *m, int i0, int j0)
{
  int n = m->n, i, j;

  for ( i = 0; i < n; i++ ) {
    /* skip a particle that is unaffected by the collision */
    if ( i != i0 && m->partner[i] != i0
      && i != j0 && m->partner[i] != j0 )
      continue;

    /* update the collision time of i */
    m->coltime[i] = INFTIME;
    for ( j = 0; j < n; j++ )
      if ( j != i )
        hsmd_coltimeij(m, i, j, m->coltime, m->partner);
  }
}



static real hsmd_step(hsmd_t *m)
{
  int i, j;

  hsmd_getij(m, &i, &j);
  hsmd_move(m, m->coltime[i]);
  hsmd_collide(m, i, j);
  hsmd_update(m, i, j);
  m->time += m->coltime[i];
  //hsmd_verifycoltime(m);
  //printf("step i %d, %g, j %d, %g, %g\n", i, m->coltime[i], j, m->coltime[j], m->time);
  return 0;
}





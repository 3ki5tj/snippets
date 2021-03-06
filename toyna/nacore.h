#ifndef NACORE_H__
#define NACORE_H__



#ifndef D
#define D 3
#endif



#include "mtrand.h"
#include "nautil.h"
#include "napot.h"
#include <ctype.h>



#define NA_UIS 1
#define NA_BIS 2
#define NA_TIS 3

#define TIS_IP 0
#define TIS_IS 1
#define TIS_IB 2

/* bond length parameters */
const double R0_PS = 4.660;
const double K0_PS = 23.0;
const double R0_SP = 3.766;
const double K0_SP = 64.0;
const double R0_SB[4] = {4.700, 4.157, 4.811, 4.163};
const double K0_SB = 10.0;

/* bond angle parameters */
const double A0_PSB[4] = {D2R( 94.06), D2R( 86.52), D2R( 98.76), D2R( 86.31)};
const double A0_BSP[4] = {D2R(106.69), D2R(108.27), D2R(106.17), D2R(108.28)};
const double A0_PSP = D2R(83.53);
const double KA_PSB = 5.0;
const double KA_BSP = 5.0;
const double KA_PSP = 20.0;

/* WCA parameters */
#define WCA_SIG 3.2
const double WCA_SIG2 = WCA_SIG * WCA_SIG;
const double WCA_EPS = 1.0;

/* stacking parameters */
const double ST_R0[4][4] = {
/*         A       C       G       U     */
/* A */  {4.164,  3.832,  4.450,  3.822},
/* C */  {4.667,  4.241,  4.992,  4.230},
/* G */  {3.971,  3.661,  4.236,  3.651},
/* U */  {4.675,  4.250,  5.000,  4.237}
/* Note U-A, A-U, C-U are unavailable, and deduced from
 *      C-A, A-C, C-C, respectively */
};
const double ST_PHI10 = D2R(-148.16);
const double ST_PHI20 = D2R( 175.97);
const double ST_KR = 1.4;
const double ST_KPHI = 4.0;

/* from Table 1. Denesyuk, 2013 */
const double ST_TM[4][4] = {
/*         A       C       G       U     */
/* A */  { 26.0,   26.0,   68.0,   26.0},
/* C */  { 26.0,   13.0,   42.0,   13.0},
/* G */  { 68.0,   70.0,   93.0,   65.0},
/* U */  { 26.0,   13.0,   65.0,  -21.0}
};
/* from Table 2. Denesyuk, 2013 */
const double ST_H[4][4] = {
/*         A       C       G       U     */
/* A */  { 4.35,   4.31,   5.12,   4.31},
/* C */  { 4.29,   4.01,   4.60,   3.99},
/* G */  { 5.08,   5.07,   5.56,   4.98},
/* U */  { 4.29,   3.99,   5.03,   3.37}
};
/* from Table 2. Denesyuk, 2013 */
const double ST_S[4][4] = {
/*         A       C       G       U     */
/* A */  {-0.32,  -0.32,   5.30,  -0.32},
/* C */  {-0.32,  -1.57,   0.77,  -1.57},
/* G */  { 5.30,   4.37,   7.35,   2.92},
/* U */  {-0.32,  -1.57,   2.92,  -3.56}
};

/* hydrogen-bond parameters */
const double HB_R0[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  5.801},
/* C */  { 0.000,  0.000,  5.548,  0.000},
/* G */  { 0.000,  5.548,  0.000,  0.000},
/* U */  { 5.801,  0.000,  0.000,  0.000}
};
const double HB_TH10[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  2.678},
/* C */  { 0.000,  0.000,  2.416,  0.000},
/* G */  { 0.000,  2.808,  0.000,  0.000},
/* U */  { 2.461,  0.000,  0.000,  0.000}
};
const double HB_TH20[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  2.461},
/* C */  { 0.000,  0.000,  2.808,  0.000},
/* G */  { 0.000,  2.416,  0.000,  0.000},
/* U */  { 2.678,  0.000,  0.000,  0.000}
};
const double HB_PSI0[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  0.878},
/* C */  { 0.000,  0.000,  0.969,  0.000},
/* G */  { 0.000,  0.969,  0.000,  0.000},
/* U */  { 0.878,  0.000,  0.000,  0.000}
};
const double HB_PSI10[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  2.670},
/* C */  { 0.000,  0.000,  2.803,  0.000},
/* G */  { 0.000,  2.546,  0.000,  0.000},
/* U */  { 2.760,  0.000,  0.000,  0.000}
};
const double HB_PSI20[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  2.670},
/* C */  { 0.000,  0.000,  2.803,  0.000},
/* G */  { 0.000,  2.546,  0.000,  0.000},
/* U */  { 2.760,  0.000,  0.000,  0.000}
};
const double HB_MUL[4][4] = {
/*         A       C       G       U     */
/* A */  { 0.000,  0.000,  0.000,  2.000},
/* C */  { 0.000,  0.000,  3.000,  0.000},
/* G */  { 0.000,  3.000,  0.000,  2.000},
/* U */  { 2.000,  0.000,  2.000,  0.000}
};
const double HB_KR = 5;
const double HB_KTH = 1.5;
const double HB_KPSI = 0.15;



typedef struct {
  int nr; /* number of residues */
  char *seq; /* sequence */
  int *iseq; /* A: 0, C: 1, G: 2, U: 3 */

  int apr; /* number of atoms per residue */
  int na; /* number of atoms */
  int dof; /* degrees of freedom */

  double tp; /* temperature */
  double debyel; /* Debye screen length */
  double uhb0;

  double *m; /* mass */
  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double ekin;
  double epot;
} na_t;



/* initialize an RNA chain for the two-interaction-site model */
static void na_initchain2(na_t *na)
{
  int i, ib, is, nr = na->nr;
  double rb = 10.0, rs = 7.0, ang = 2*PI/10, dh = 3.4, th, c, s;

  for ( i = 0; i < nr; i++ ) {
    th = ang * i;
    c = cos(th);
    s = sin(th);
    ib = i*2;
    is = i*2 + 1;
    na->m[ib] = 1.0; // TODO
    na->x[ib][0] = rb * c;
    na->x[ib][1] = rb * s;
    na->x[ib][2] = dh * i;
    na->m[is] = 1.0; // TODO
    na->x[is][0] = rs * c;
    na->x[is][1] = rs * s;
    na->x[is][2] = dh * i;
  }
}



/* initialize an RNA chain for the three-interaction-site model */
static void na_initchain3(na_t *na)
{
  int i, ic, nr = na->nr;
  double ang = D2R(32.7), dh = 2.81;
  const double rp = 8.710, thp = D2R(-70.502), dhp = 3.750;
  const double rs = 9.214, ths = D2R(-41.097), dhs = 2.864;
  const double rbarr[4] = {5.458, 5.643, 5.631, 5.633};
  const double thbarr[4] = {D2R(-25.976), D2R(-33.975), D2R(-22.124), D2R(-34.102)};
  const double dhbarr[4] = {0.742, 0.934, 0.704, 0.932};
  double rb, thb, dhb, th;
  const double bmass[4] = {134.0, 110.0, 150.0, 111.0};

  for ( i = 0; i < nr; i++ ) {
    th = ang * i;

    na->m[i*3 + TIS_IP] = 95.0;
    na->x[i*3 + TIS_IP][0] = rp * cos(th + thp);
    na->x[i*3 + TIS_IP][1] = rp * sin(th + thp);
    na->x[i*3 + TIS_IP][2] = dh * i + dhp;

    na->m[i*3 + TIS_IS] = 99.0;
    na->x[i*3 + TIS_IS][0] = rs * cos(th + ths);
    na->x[i*3 + TIS_IS][1] = rs * sin(th + ths);
    na->x[i*3 + TIS_IS][2] = dh * i + dhs;

    ic = na->iseq[ i ];
    rb = rbarr[ ic ];
    thb = thbarr[ ic ];
    dhb = dhbarr[ ic ];
    na->m[i*3 + TIS_IB] = bmass[ ic ];
    na->x[i*3 + TIS_IB][0] = rb * cos(th + thb);
    na->x[i*3 + TIS_IB][1] = rb * sin(th + thb);
    na->x[i*3 + TIS_IB][2] = dh * i + dhb;
  }
}



/* open an LJ system */
static na_t *na_open(const char *seq, int model,
    double tp, double debyel, double uhb0)
{
  na_t *na;
  int i, d, n, nr;

  xnew(na, 1);

  na->nr = nr = strlen(seq);
  xnew(na->seq, nr);
  xnew(na->iseq, nr);
  for ( i = 0; i < nr; i++ ) {
    na->seq[i] = toupper( seq[i] );
    if ( na->seq[i] == 'A' ) {
      na->iseq[i] = 0;
    } else if ( na->seq[i] == 'C' ) {
      na->iseq[i] = 1;
    } else if ( na->seq[i] == 'G' ) {
      na->iseq[i] = 2;
    } else if ( na->seq[i] == 'U' || na->seq[i] == 'T' ) {
      na->iseq[i] = 3;
    } else {
      fprintf(stderr, "unknown residue type %c at %d\n", na->seq[i], i);
      free(na->seq);
      free(na->iseq);
      free(na);
      return NULL;
    }
  }

  na->apr = model;
  na->na = n = nr * na->apr;
  na->dof = n * D - D * (D+1)/2;

  na->tp = tp;
  na->debyel = debyel;
  na->uhb0 = uhb0;

  xnew(na->m, n);
  xnew(na->x, n);
  xnew(na->v, n);
  xnew(na->f, n);

  if ( model == 2 ) {
    na_initchain2(na);
  } else {
    na_initchain3(na);
  }

  /* initalize random velocities */
  for (i = 0; i < n; i++) {
    double amp = sqrt( BOLTZK * tp / na->m[i] );
    for ( d = 0; d < D; d++ ) {
      na->v[i][d] = amp * randgaus();
    }
  }

  rmcom(na->v, na->m, n);
  shiftang(na->x, na->v, na->m, n);

  return na;
}



static void na_close(na_t *na)
{
  free(na->m);
  free(na->x);
  free(na->v);
  free(na->f);
  free(na);
}



/* extract the sequence from the PDB file */
__inline static char *na_getseqpdb(const char *fn, int *resmin)
{
  char *seq, buf[128];
  int nr, ir, irmin, irmax;
  FILE *fp;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return NULL;
  }

  /* try to determine the number of residues */
  irmin =  10000000;
  irmax = -10000000;
  while ( fgets(buf, sizeof buf, fp) ) {
    /* we will skip residues without the P atom */
    buf[26] = '\0';
    if ( strncmp(buf, "ATOM  ", 6) != 0
      || strncmp(buf + 12, " P  ", 4) != 0 ) {
      continue;
    }
    buf[26] = '\0';
    ir = atoi(buf+22);
    if ( ir > irmax ) {
      irmax = ir;
    }
    if ( ir < irmin ) {
      irmin = ir;
    }
  }
  if ( resmin != NULL ) {
    *resmin = irmin;
  }
  nr = irmax - irmin + 1;

  /* allocate space for the sequence */
  if ( (seq = malloc(nr + 1)) == NULL ) {
    fclose(fp);
    return NULL;
  }

  /* load the sequence */
  rewind(fp);
  while ( fgets(buf, sizeof buf, fp) ) {
    if ( strncmp(buf, "ATOM  ", 6) != 0
      || strncmp(buf + 12, " P  ", 4) != 0 ) {
      continue;
    }
    buf[26] = '\0';
    ir = atoi(buf+22) - irmin;
    seq[ir] = buf[19];
  }
  seq[nr] = '\0';
  fclose(fp);
  return seq;
}



/* load coordinates from a pdb */
__inline static int na_loadpdb(na_t *na, const char *fn, int irmin)
{
  char buf[128];
  int ir, id;
  FILE *fp;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    /* we will skip residues without the P atom */
    if ( strncmp(buf, "ATOM  ", 6) != 0 ) {
      continue;
    }

    /* copy the residue id */
    buf[26] = '\0';
    ir = atoi(buf+22) - irmin;
    if ( ir < 0 || ir >= na->nr ) {
      continue;
    }

    if ( buf[13] == 'P' ) {
      id = ir*3 + TIS_IP;
    } else if ( buf[13] == 'S' ) {
      id = ir*3 + TIS_IS;
    } else if ( buf[13] == 'B' ) {
      id = ir*3 + TIS_IB;
    } else {
      continue;
    }

    /* copy the coordinates */
    buf[54] = '\0';
    na->x[id][2] = atof( buf+46 );
    buf[46] = '\0';
    na->x[id][1] = atof( buf+38 );
    buf[38] = '\0';
    na->x[id][0] = atof( buf+30 );
    //printf("ir %d, id %d, %g %g %g\n", ir, id, na->x[id][0], na->x[id][1], na->x[id][2]);
  }
  return 0;
}



/* write positions (and possibly velocities) */
__inline static int na_writepos(na_t *na,
    double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d, n = na->na;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %d\n", D, na->nr, na->na, (v != NULL));
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ )
      fprintf(fp, "%.14e ", x[i][d]);
    fprintf(fp, "%d ", i % 2);
    if ( v != NULL )
      for ( d = 0; d < D; d++ )
        fprintf(fp, "%.14e ", v[i][d]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



#define na_energy(na) \
  na->epot = na_energyTIS_low(na, na->x, na->tp, na->debyel)

/* compute force, return energy */
__inline static double na_energyTIS_low(na_t *na, double (*x)[D],
    double tp, double debyel)
{
  int i, j, ic, jc, nr = na->nr, n = na->na;
  double ep = 0;
  double eps, Q, QQ;

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS],
                   NULL, NULL);
    ep += ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   NULL, NULL);
    if ( i < nr - 1 ) {
      ep += ebondlen(R0_SP, K0_SP,
                     x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL);
    }
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    ep += ebondang(A0_PSB[ic], KA_PSB,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   NULL, NULL, NULL);
    if ( i < nr - 1 ) {
      ep += ebondang(A0_BSP[ic], KA_BSP,
                     x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL, NULL);
      ep += ebondang(A0_PSP, KA_PSP,
                     x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     NULL, NULL, NULL);
    }
  }

  // excluded volume
  for ( i = 0; i < n - 1; i++ ) {
    for (j = i + 1; j < n; j++) {
      ep += ewca(WCA_SIG2, WCA_EPS, x[i], x[j], NULL, NULL);
    }
  }

  // stack energy
  for ( i = 0; i < nr - 1; i++ ) {
    double ust0 = 0, *xp3 = NULL;
    ic = na->iseq[i];
    jc = na->iseq[i + 1];
    ust0 = -ST_H[ic][jc] + ST_S[ic][jc] * (tp - (T0 + ST_TM[ic][jc]));
    if ( i < nr - 2 ) {
      xp3 = x[(i + 2)*3 + TIS_IP];
    }
    ep += estack(ST_R0[ic][jc], ST_PHI10, ST_PHI20, ST_KR, ST_KPHI, ust0,
                 x[i*3       + TIS_IP], x[i*3       + TIS_IS], x[i*3       + TIS_IB],
                 x[(i + 1)*3 + TIS_IP], x[(i + 1)*3 + TIS_IS], x[(i + 1)*3 + TIS_IB], xp3,
                 NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  }

  // hydrogen-bond energy
  for ( i = 0; i < nr - 1; i++ ) {
    ic = na->iseq[i];
    for ( j = i + 2; j < nr; j++ ) {
      jc = na->iseq[j];
      if ( ic + jc != 3 ) {
        continue;
      }
      ep += ehbond(HB_R0[ic][jc], HB_TH10[ic][jc], HB_TH20[ic][jc],
                   HB_PSI0[ic][jc], HB_PSI10[ic][jc], HB_PSI20[ic][jc],
                   HB_KR, HB_KTH, HB_KPSI, na->uhb0 * HB_MUL[ic][jc],
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   x[j*3 + TIS_IB], x[j*3 + TIS_IS], x[j*3 + TIS_IP],
                   NULL, NULL, NULL, NULL, NULL, NULL);
    }
  }

  // electrostatic interaction
  Q = getchargeQ(tp, &eps);
  QQ = Q*Q*KE2/eps;
  for ( i = 0; i < nr - 1; i++ ) {
    for ( j = i + 1; j < nr; j++ ) {
      ep += echargeDH(QQ, debyel,
                      x[i*3 + TIS_IP], x[j*3 + TIS_IP],
                      NULL, NULL);
    }
  }
  return ep;
}



#define na_force(na) \
  na->epot = na_forceTIS_low(na, na->x, na->f, na->tp, na->debyel)

/* compute force, return energy */
__inline static double na_forceTIS_low(na_t *na, double (*x)[D], double (*f)[D],
    double tp, double debyel)
{
  int i, j, ic, jc, nr = na->nr, n = na->na;
  double v[3], ep = 0, Q, QQ, eps;

  for (i = 0; i < n; i++) vzero(f[i]);

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    v[0] = ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS],
                   f[i*3 + TIS_IP], f[i*3 + TIS_IS]);
    v[1] = ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      v[2] = ebondlen(R0_SP, K0_SP,
                      x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                      f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    } else {
      v[2] = 0;
    }
    ep += v[0] + v[1] + v[2];
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = na->iseq[i];
    v[0] = ebondang(A0_PSB[ic], KA_PSB,
                    x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                    f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      v[1] = ebondang(A0_BSP[ic], KA_BSP,
                      x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                      f[i*3 + TIS_IB], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
      v[2] = ebondang(A0_PSP, KA_PSP,
                      x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                      f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    } else {
      v[1] = v[2] = 0;
    }
    ep += v[0] + v[1] + v[2];
  }

  // excluded volume
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      v[0] = ewca(WCA_SIG2, WCA_EPS, x[i], x[j], f[i], f[j]);
      ep += v[0];
    }
  }

  // stack energy
  for ( i = 0; i < nr - 1; i++ ) {
    double ust0 = 0, *xp3 = NULL, *fp3 = NULL;
    ic = na->iseq[i];
    jc = na->iseq[i + 1];
    ust0 = -ST_H[ic][jc] + ST_S[ic][jc] * (tp - (T0 + ST_TM[ic][jc]));
    if ( i < nr - 2 ) {
      xp3 = x[(i + 2)*3 + TIS_IP];
      fp3 = f[(i + 2)*3 + TIS_IP];
    }
    v[0] = estack(ST_R0[ic][jc], ST_PHI10, ST_PHI20, ST_KR, ST_KPHI, ust0,
                  x[i*3       + TIS_IP], x[i*3       + TIS_IS], x[i*3       + TIS_IB],
                  x[(i + 1)*3 + TIS_IP], x[(i + 1)*3 + TIS_IS], x[(i + 1)*3 + TIS_IB], xp3,
                  f[i*3       + TIS_IP], f[i*3       + TIS_IS], f[i*3       + TIS_IB],
                  f[(i + 1)*3 + TIS_IP], f[(i + 1)*3 + TIS_IS], f[(i + 1)*3 + TIS_IB], fp3);
    ep += v[0];
  }

  // hydrogen-bond energy
  for ( i = 0; i < nr - 1; i++ ) {
    ic = na->iseq[i];
    for ( j = i + 2; j < nr; j++ ) {
      jc = na->iseq[j];
      if ( ic + jc != 3 ) {
        continue;
      }
      v[0] = ehbond(HB_R0[ic][jc], HB_TH10[ic][jc], HB_TH20[ic][jc],
                    HB_PSI0[ic][jc], HB_PSI10[ic][jc], HB_PSI20[ic][jc],
                    HB_KR, HB_KTH, HB_KPSI, na->uhb0 * HB_MUL[ic][jc],
                    x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                    x[j*3 + TIS_IB], x[j*3 + TIS_IS], x[j*3 + TIS_IP],
                    f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[i*3 + TIS_IB],
                    f[j*3 + TIS_IB], f[j*3 + TIS_IS], f[j*3 + TIS_IP]);
      //fprintf(stderr, "%4d(%c)-%4d(%c): %g\n", i, na->seq[i], j, na->seq[j], v[0]);
      ep += v[0];
    }
  }

  // electrostatic interaction
  Q = getchargeQ(tp, &eps);
  QQ = Q*Q*KE2/eps;
  for ( i = 0; i < nr - 1; i++ ) {
    for ( j = i + 1; j < nr; j++ ) {
      ep += echargeDH(QQ, debyel,
                      x[i*3 + TIS_IP], x[j*3 + TIS_IP],
                      f[i*3 + TIS_IP], f[j*3 + TIS_IP]);
    }
  }
  return ep;
}



/* velocity-verlet */
__inline static void na_vv(na_t *na, double dt)
{
  int i, n = na->na;
  double dth = dt * 0.5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(na->v[i], na->f[i], dth/na->m[i]);
    vsinc(na->x[i], na->v[i], dt);
  }
  na_force(na);
  for (i = 0; i < n; i++) { /* VV part 2 */
    vsinc(na->v[i], na->f[i], dth/na->m[i]);
  }
}



/* compute the kinetic energy */
#define na_ekin(na, v) md_ekin(v, na->m, na->na)



/* velocity rescaling thermostat */
#define na_vrescale(na, tp, dt) \
  md_vrescale(na->v, na->m, na->na, na->dof, BOLTZK * tp, dt)



#if 0



/* displace a random particle i, return i */
static int na_randmv(na_t *na, double *xi, double amp)
{
  int i, d;

  i = (int) (rand01() * na->n);
  for ( d = 0; d < D; d++ )
    xi[d] = na->x[i][d] + (rand01() * 2 - 1) * amp;
  return i;
}



/* compute pair energy */
static int na_pair(double *xi, double *xj, double l, double invl,
    double rc2, double *u)
{
  double dx[D], dr2, invr2, invr6;

  dr2 = na_dist2(dx, xi, xj);
  if (dr2 < rc2) {
    invr2 = 1 / dr2;
    invr6 = invr2 * invr2 * invr2;
    *u  = 4.f * invr6 * (invr6 - 1);
    return 1;
  } else {
    *u = 0;
    return 0;
  }
}



/* return the energy change from displacing x[i] to xi */
__inline static double na_depot(na_t *na, int i, double *xi)
{
  int j, n = na->n;
  double l = na->l, invl = 1/l, rc2 = na->rc2, u, du;

  u = 0;
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) continue;
    if ( na_pair(na->x[i], na->x[j], rc2, &du) ) {
      u -= du;
    }
    if ( na_pair(xi, na->x[j], rc2, &du) ) {
      u += du;
    }
  }
  return u;
}



/* commit a particle displacement */
__inline static void na_commit(na_t *na, int i,
    const double *xi, double du)
{
  vcopy(na->x[i], xi);
  na->ep0 += du;
  na->epot += du;
}



/* Metropolis algorithm */
__inline static int na_metro(na_t *na, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0.;

  i = na_randmv(na, xi, amp);
  du = na_depot(na, i, xi);
  if ( du < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -bet * du ) );
  }
  if ( acc ) {
    na_commit(na, i, xi, du);
    return 1;
  }
  return 0;
}
#endif


#endif /* NACORE_H__ */

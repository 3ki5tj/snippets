#ifndef CAGOCORE_H__
#define CAGOCORE_H__



/* alpha-carbon based Go-model */



#include "cagoutil.h"



typedef struct {
  int n; /* number of residues */
  int dof; /* degree of freedom */
  unsigned flags; /* input flags */

  double kb; /* .5 kb (b - b0)^2, ~ 200 */
  double ka; /* .5 ka (a - a0)^2, ~ 40 */
  double kd1, kd3; /* kd1 (1 - cos(d - d0)) + kd3 (1 - cos(3*(d-d0))), ~ 1 & 0.5 */
  double nbe, nbc; /* nbc ~ 4 A */

  double (*xref)[3];
  double epotref; /* energy of the reference structure */
  int *iaa; /* integer amino-acid type */
  int *ires; /* residue */

  double *bref; /* bonds */
  double *aref; /* angles */
  double *dref; /* dihedrals */
  double *r2ref; /* pair distances */

  int ncont; /* number of defined contacts */
  int *iscont; /* if the pair i-j is a contact */

  /* variables for MD simulations */
  double *m; /* masses */
  double (*x)[3];
  double (*v)[3];
  double (*f)[3];
  double (*x1)[3];
  double ekin, epot;
} cago_t;



/* convenient macro for computing RMSD from the reference structure */
#define cago_rmsd(go, x, xf) \
  vrmsd(x, xf, go->xref, go->m, go->n, 0, NULL, NULL)



/* compute the reference bond lengths, angles, dihedrals and pair distances */
__inline static int cago_refgeo(cago_t *go)
{
  int i, j, n = go->n;
  double dx[D];

  /* calculate reference bond lengths, angles, dihedrals */
  xnew(go->bref, n - 1); /* bonds */
  for (i = 0; i < n - 1; i++) {
    go->bref[i] = vdistx(dx, go->xref[i], go->xref[i + 1]);
  }

  xnew(go->aref, n - 2); /* angles */
  for (i = 1; i < n - 1; i++) {
    go->aref[i - 1]  = vang(go->xref[i - 1], go->xref[i], go->xref[i + 1],
      NULL, NULL, NULL);
  }

  xnew(go->dref, n - 3); /* dihedrals */
  for (i = 0; i < n - 3; i++) {
    go->dref[i] = vdih(go->xref[i], go->xref[i + 1],
        go->xref[i + 2], go->xref[i + 3], NULL, NULL, NULL, NULL);
  }

  /* reference pair distances */
  xnew(go->r2ref, n*n);
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      go->r2ref[j*n + i] = go->r2ref[i*n + j]
          = vsqr( vdiff(dx, go->xref[i], go->xref[j]) );
    }
  }
  return 0;
}



enum { PDB_CONTACT_CA, PDB_CONTACT_HEAVY, PDB_CONTACT_ALL }; /* ways of searching contacts */



/* read C-alpha coordinates, and amino acid types from a PDB file */
__inline static int cago_loadbasic(cago_t *go, const char *fn)
{
  int nres = 0, ncap;
  const int blksz = 32;
  char s[128];
  FILE *fp;

  ncap = blksz;
  xnew(go->xref, ncap);
  xnew(go->iaa, ncap);
  xnew(go->ires, ncap);

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  while ( fgets(s, sizeof s, fp) ) {
    /* we only load a single model */
    if ( strncmp(s, "TER", 3) == 0
      || strncmp(s, "END", 3) == 0
      || strncmp(s, "ENDMDL", 6) == 0 ) {
      break;
    }

    /* we are only interested in C-alpha atom */
    if ( strncmp(s, "ATOM  ", 6) != 0
      || strncmp(s + 13, "CA", 2) != 0 ) {
      continue;
    }

    /* discard alternative positions */
    if ( s[16] != ' ' && s[16] != 'A' ) {
      continue;
    }

    /* reallocate space if needed */
    if ( nres + 1 >= ncap ) {
      ncap += blksz;
      xrenew(go->xref, ncap);
      xrenew(go->iaa, ncap);
      xrenew(go->ires, ncap);
    }

    /* load the coordinates */
    if ( 3 != sscanf(s + 30, "%lf%lf%lf",
                     &go->xref[nres][0],
                     &go->xref[nres][1],
                     &go->xref[nres][2]) ) {
      fprintf(stderr, "%s: cannot load coordinates for residue %d\n", fn, nres);
      break;
    }

    s[20] = '\0';
    go->iaa[nres] = res2iaa(s + 17);

    s[26] = '\0';
    go->ires[nres] = atoi(s + 22);

    nres += 1;
  }

  fclose(fp);
  return go->n = nres;
}



/* make the contact map */
__inline static int cago_mkcont(cago_t *go, const char *fn,
    double rc, int ctype, int nsexcl)
{
#define RAMAX 32 /* maximal number of atoms in the residues */
  int ires, ir, jr, rai, n = go->n;
  int *ra; /* number of atoms in residue i */
  double (*x)[RAMAX][3];
  double rc2 = rc * rc, dx[D];
  char s[128];
  FILE *fp;

  xnew(go->iscont, n * n);
  xnew(ra, n);
  xnew(x, n);

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  while ( fgets(s, sizeof s, fp) ) {
    /* we only load a single model */
    if ( strncmp(s, "TER", 3) == 0
      || strncmp(s, "END", 3) == 0
      || strncmp(s, "ENDMDL", 6) == 0 ) {
      break;
    }

    /* we are only interested in C-alpha atom */
    if ( strncmp(s, "ATOM  ", 6) != 0 ) {
      continue;
    }

    /* discard alternative positions */
    if ( s[16] != ' ' && s[16] != 'A' ) {
      continue;
    }

    if ( ctype == PDB_CONTACT_CA
      && strncmp(s + 13, "CA", 2) != 0 ) {
      continue;
    }

    if ( ctype == PDB_CONTACT_HEAVY
      && s[13] == 'H' ) {
      continue;
    }

    s[26] = '\0';
    ires = atoi(s + 22);

    /* find the actual residue index */
    for ( ir = 0; ir < go->n; ir++ ) {
      if ( go->ires[ir] == ires ) {
        break;
      }
    }
    if ( ir >= go->n ) {
      fprintf(stderr, "%s: unknown residue index %d\n", fn, ires);
      continue;
    }

    if ( ra[ir] >= RAMAX ) {
      fprintf(stderr, "%s: too many atoms for residue %d/%d\n", fn, ir, ires);
      continue;
    }
    rai = ra[ir];
    ra[ir] += 1;

    /* load the coordinates */
    if ( 3 != sscanf(s + 30, "%lf%lf%lf", &x[ir][rai][0], &x[ir][rai][1], &x[ir][rai][2]) ) {
      fprintf(stderr, "%s: cannot load coordinates for residue %d\n", fn, ires);
      break;
    }
  }

  fclose(fp);

  for ( ir = 0; ir < n; ir++ ) {
    for ( jr = ir + nsexcl; jr < n; jr++ ) {
      int ia, ja, isc = 0;

      /* scan over atoms in the two residues */
      for ( ia = 0; ia < ra[ir] && !isc; ia++ ) {
        for ( ja = 0; ja < ra[jr] && !isc; ja++ ) {
          if ( vsqr( vdiff(dx, x[ir][ia], x[jr][ja]) ) < rc2 )
            isc = 1;
        }
      }

      go->iscont[ir*n + jr] = go->iscont[jr*n + ir] = isc;
    }
  }

  free(ra);
  free(x);
  return 0;
}




/* return cago_t from pdb file fnpdb
 * `rcc' is the cutoff radius for defining contacts
 * `ctype' is one of PDB_CONTACT_CA, _HEAVY, _ALL
 * `nsexcl' is the number of successive residues to be excluded as contacts
 * e.g., nsexcl = 4 means `a' and `d' in -a-b-c-d- are excluded */
__inline static cago_t *cago_open0(const char *fnpdb,
    double rcc, int ctype, int nsexcl)
{
  cago_t *go;
  int i, j;

  xnew(go, 1);
  /* read coordinates, and amino acid types */
  if ( cago_loadbasic(go, fnpdb) < 0 ) {
    fprintf(stderr, "cannot load %s\n", fnpdb);
    free(go);
    return NULL;
  }

  go->dof = go->n*3 - 6;

  cago_mkcont(go, fnpdb, rcc, ctype, nsexcl);

  xnew(go->m, go->n);
  for ( i = 0; i < go->n; i++ ) {
    go->m[i] = 1.0;
  }

  /* compute the reference bond length, angles, etc. */
  cago_refgeo(go);

  /* count the number of contacts */
  for (go->ncont = 0, i = 0; i < go->n - 1; i++)
    for (j = i + 1; j < go->n; j++)
      go->ncont += go->iscont[ i*go->n + j ];
  return go;
}



/* set force parameters */
__inline static void cago_setfparam(cago_t *go,
    double kb, double ka, double kd1, double kd3,
    double nbe, double nbc)
{
  go->kb = kb;
  go->ka = ka;
  go->kd1 = kd1;
  go->kd3 = kd3;
  go->nbe = nbe;
  go->nbc = nbc;
}



/* return a pointer to cago_t from PDB file `fnpdb'
 * cago_open1() + bond parameters */
__inline static cago_t *cago_open(const char *fnpdb,
    double kb, double ka, double kd1, double kd3, double nbe, double nbc,
    double rcc, int ctype, int nsexcl)
{
  cago_t *go = cago_open0(fnpdb, rcc, ctype, nsexcl);
  if (go == NULL) return NULL;
  cago_setfparam(go, kb, ka, kd1, kd3, nbe, nbc);
  return go;
}



/* destroy cago_t */
__inline static void cago_close(cago_t *go)
{
  free(go->m);
  free(go->x);
  free(go->v);
  free(go->f);
  free(go->x1);
  free(go->iscont);
  free(go->bref);
  free(go->aref);
  free(go->dref);
  free(go->r2ref);
  free(go->xref);
  free(go->iaa);
  free(go->ires);
  free(go);
}



/* bond energy 1/2 k (r - r0)^2 */
__inline static double potbond(vct a, vct b, double r0, double k,
    vct fa, vct fb)
{
  double dx[3], r, dr;

  r = vnorm( vdiff(dx, a, b) );
  dr = r - r0;
  if (fa != NULL) {
    double amp = k * dr / r;
    vsinc(fa, dx, -amp);
    vsinc(fb, dx,  amp);
  }
  return 0.5 * k * dr * dr;
}



/* harmonic angle 1/2 k (ang - ang0)^2 */
__inline static double potang(vct a, vct b, vct c, double ang0, double k,
    vct fa, vct fb, vct fc)
{
  double dang, amp, ga[3], gb[3], gc[3];

  if (fa) { /* compute gradient */
    dang = vang(a, b, c, ga, gb, gc) - ang0;
    amp = -k * dang;
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
  } else {
    dang = vang(a, b, c, NULL, NULL, NULL) - ang0;
  }
  return .5f * k * dang * dang;
}



/* 1-3 dihedral: k1 * (1 - cos(dang)) + k3 * (1 - cos(3*dang)) */
__inline static double potdih13(vct a, vct b, vct c, vct d, double ang0,
    double k1, double k3, vct fa, vct fb, vct fc, vct fd)
{
  double dang, amp, ga[3], gb[3], gc[3], gd[3];

  if (fa) {
    dang = vdih(a, b, c, d, ga, gb, gc, gd) - ang0;
    amp  = (double)( -k1 * sin(dang) - 3 * k3 * sin(3*dang) );
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
    vsinc(fd, gd, amp);
  } else {
    dang = vdih(a, b, c, d, NULL, NULL, NULL, NULL) - ang0;
  }
  return (double)( k1 * (1 - cos(dang)) + k3 * (1 - cos(3 * dang)) );
}



/* 12-10 potential: u = 5(rc/r)^12 - 6(rc/r)^10,
 * the minimum is at r = rc, and u = -1 */
__inline static double pot1210(vct a, vct b, double rc2, double eps, vct fa, vct fb)
{
  double dx[3], dr2, invr2, invr4, invr6, invr10, amp;

  dr2 = vsqr( vdiff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr4 = invr2 * invr2;
  invr6 = invr4 * invr2;
  invr10 = invr4 * invr6;
  if (fa) {
    amp = 60 * eps * (invr2 - 1) * invr10 * (1/dr2);
    vsinc(fa, dx,  amp);
    vsinc(fb, dx, -amp);
  }
  return eps * (5 * invr2 - 6) * invr10;
}



/* repulsive potential: (rc/r)^12 */
__inline static double potr12(vct a, vct b, double rc2, double eps, vct fa, vct fb)
{
  double dx[3], dr2, invr2, invr6, u, amp;

  dr2 = vsqr( vdiff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;
  u = eps * invr6 * invr6;
  if (fa) {
    amp = 12 * u / dr2;
    vsinc(fa, dx,  amp);
    vsinc(fb, dx, -amp);
  }
  return u;
}



/* force field from C. Clementi, H. Nymeyer, J. N. Onuchic,
 * J. Mol. Biol, Vol. 298 (2000) 937-953 */
__inline static double cago_force(cago_t *go, vct *x, vct *f)
{
  int i, j, id, n = go->n;
  double ene = 0, kb = go->kb, ka = go->ka, kd1 = go->kd1, kd3 = go->kd3;
  double nbe = go->nbe, nbc2 = go->nbc * go->nbc;

  if (f != NULL) {
    for (i = 0; i < n; i++) {
      vzero(f[i]);
    }
  }

  /* bonds */
  for (i = 0; i < n - 1; i++)
    ene += potbond(x[i], x[i + 1], go->bref[i], kb, f[i], f[i + 1]);

  /* angles */
  for (i = 0; i < n - 2; i++)
    ene += potang(x[i], x[i + 1], x[i + 2], go->aref[i],
              ka, f[i], f[i + 1], f[i + 2]);

  /* dihedrals */
  for (i = 0; i < n - 3; i++)
    ene += potdih13(x[i], x[i + 1], x[i + 2], x[i + 3], go->dref[i],
          kd1, kd3, f[i], f[i + 1], f[i + 2], f[i + 3]);

  /* non-bonded */
  for (i = 0; i < n - 4; i++)
    for (j = i + 4; j < n; j++) {
      id = i*n + j;
      if ( go->iscont[id] ) { /* contact pair */
        ene += pot1210(x[i], x[j], go->r2ref[id], nbe, f[i], f[j]);
      } else {
        ene += potr12(x[i], x[j], nbc2, nbe, f[i], f[j]);
      }
    }
  return ene;
}



/* remove center of mass motion, linear and angular */
__inline static void cago_rmcom(cago_t *go, double (*x)[D], double (*v)[D])
{
  rmcom(v, go->m, go->n);
  shiftang(x, v, go->m, go->n);
}



/* compute the kinetic energy */
__inline static int cago_ekin(cago_t *go, double (*v)[D])
{
  int i, n = go->n;
  double ek = 0.0;

  for ( i = 0; i < n; i++ ) {
    ek += go->m[i] * vsqr( v[i] );
  }
  return 0.5 * ek;
}



/* initialize molecular dynamics
 *  o create an initial structure
 *    if `open', start from a nearly-straight chain,
 *      with a disturbance of `rndamp' in the x, y directions
 *    otherwise start from the reference structure,
 *      with a random disturbance of `rndamp'
 *  o initialize the velocity with the center of mass motion removed
 *  o compute the initial force and energy */
__inline static int cago_initmd(cago_t *go,
    int open, double rndamp, double T0)
{
  int i, n = go->n;
  double s, dx[3];

  xnew(go->f, n);
  xnew(go->v, n);
  xnew(go->x, n);
  xnew(go->x1, n);

  /* initialize position */
  if (open) { /* open chain */
    for (i = 0; i < n - 1; i++) {
      dx[0] = 1.0;
      dx[1] = rndamp * randgaus();
      dx[2] = rndamp * randgaus();
      s = sqrt( vsqr( dx ) );
      /* x_{i+1} = x_i + dx * bref[i] */
      vsadd(go->x[i + 1], go->x[i], dx, go->bref[i] / s);
    }
  } else { /* copy from xref, slightly disturb it */
    for (i = 0; i < n; i++) {
      dx[0] = rndamp * randgaus();
      dx[1] = rndamp * randgaus();
      dx[2] = rndamp * randgaus();
      vadd(go->x[i], go->xref[i], dx);
    }
  }
  rmcom(go->x, go->m, go->n);
  go->epotref = cago_force(go, go->xref, go->f);

  /* initialize velocities */
  for (i = 0; i < n; i++) {
    s = sqrt( T0 / go->m[i] );
    go->v[i][0] = s * randgaus();
    go->v[i][1] = s * randgaus();
    go->v[i][2] = s * randgaus();
  }
  cago_rmcom(go, go->x, go->v); /* remove center of mass motion */
  go->ekin = cago_ekin(go, go->v);
  return 0;
}



/* velocity Verlet */
__inline static int cago_vv(cago_t *go, double fs, double dt)
{
  int i, n = go->n;
  double dth = .5f * dt * fs;

  for ( i = 0; i < n; i++ ) { /* VV part 1 */
    vsinc(go->v[i], go->f[i], dth / go->m[i]);
    vsinc(go->x[i], go->v[i], dt);
  }

  go->epot = cago_force(go, go->x, go->f); /* calculate force */

  for ( i = 0; i < n; i++ ) { /* VV part 2 */
    vsinc(go->v[i], go->f[i], dth / go->m[i]);
  }
  go->ekin = cago_ekin(go, go->v);
  return 0;
}



/* Exact velocity rescaling thermostat */
__inline static double cago_vrescale(cago_t *go,
    double (*v)[D], double tp, double dt)
{
  int i, n = go->n, dof = go->dof;
  double ekav = 0.5 * tp * dof, ek1, ek2, s, c = 0, r, r2;

  c = exp(-dt);
  ek1 = cago_ekin(go, v);
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = ek1 + (1 - c) * (ekav*(r2 + r*r)/dof - ek1)
      + 2 * r * sqrt(c*(1 - c) * ekav/dof*ek1);
  if ( ek2 < 0 ) {
    ek2 = 0;
  }
  s = sqrt(ek2 / ek1);
  for ( i = 0; i < n; i++ ) {
    vsmul(v[i], s);
  }
  return ek2;
}



/* the change of the potential energy */
__inline static double cago_depot(cago_t *go,
    double (*x)[D], int i, double *xi)
{
  int j, id, n = go->n;
  double (*xn)[D] = go->x1; /* we use x1 freely here */
  double ene = 0;
  double ka = go->ka, kb = go->kb, kd1 = go->kd1, kd3 = go->kd3;
  double nbe = go->nbe, nbc2 = go->nbc * go->nbc;

  /* copy coordinates */
  for (j = 0; j < n; j++) {
    if (j == i) vcopy(xn[i], xi);
    else vcopy(xn[j], x[j]);
  }

  /* bonds */
  for (j = i - 1; j <= i; j++) {
    if (j < 0 || j >= n - 1) continue;
    ene -= potbond(x[j], x[j + 1], go->bref[j], kb, NULL, NULL);
    ene += potbond(xn[j], xn[j + 1], go->bref[j], kb, NULL, NULL);
  }

  /* angles */
  for (j = i - 1; j <= i + 1; j++) {
    if (j < 1 || j >= n - 1) continue;
    ene -= potang(x[j - 1], x[j], x[j + 1], go->aref[j - 1], ka,
        NULL, NULL, NULL);
    ene += potang(xn[j - 1], xn[j], xn[j + 1], go->aref[j - 1], ka,
        NULL, NULL, NULL);
  }

  /* dihedrals */
  for (j = i - 3; j <= i; j++) {
    if (j < 0 || j >= n - 3) continue;
    ene -= potdih13(x[j], x[j + 1], x[j + 2], x[j + 3], go->dref[j],
        kd1, kd3, NULL, NULL, NULL, NULL);
    ene += potdih13(xn[j], xn[j + 1], xn[j + 2], xn[j + 3], go->dref[j],
        kd1, kd3, NULL, NULL, NULL, NULL);
  }

  /* non-bonded interaction */
  for (j = 0; j < n; j++) {
    if ( j > i - 4 || j < i + 4 ) continue;

    /* subtract the old energies */
    id = i*n + j;
    if ( go->iscont[id] ) { /* contact pair */
      ene -= pot1210(x[i], x[j], go->r2ref[id], nbe, NULL, NULL);
    } else { /* non-contact pair */
      ene -= potr12(x[i], x[j], nbc2, nbe, NULL, NULL);
    }

    /* add the new energies */
    if ( go->iscont[id] ) { /* contact pair */
      ene += pot1210(x[i], x[j], go->r2ref[id], nbe, NULL, NULL);
    } else { /* non-contact pair */
      ene += potr12(x[i], x[j], nbc2, nbe, NULL, NULL);
    }
  }
  return ene;
}



/* Metropolis algorithm */
__inline static int cago_metro(cago_t *go, double amp, double bet)
{
  int i;
  double du, xi[D];

  i = (int) (go->n * rand01());
  xi[0] = amp * (rand01() * 2 - 1);
  xi[1] = amp * (rand01() * 2 - 1);
  xi[2] = amp * (rand01() * 2 - 1);
  vinc(xi, go->x[i]);
  du = cago_depot(go, go->x, i, xi);
  if (du < 0 || rand01() < exp(-bet * du)) {
    vcopy(go->x[i], xi);
    go->epot += du;
    return 1;
  } else
    return 0;
}



/* count the number of native contacts that are formed in the structure `x'
 * this counting process is independent of the process of defining contacts.
 * here, given a set of defined contacts, we simple observe how many pairs
 *   are close enough to be regarded as contacts
 * a contact is formed if the pair distance is <= gam * native-distance
 * return the number of contacts
 * `*Q' is the ratio of formed contacts / the total number of contacts  */
__inline static int cago_ncontacts(cago_t *go,
    double (*x)[D], double gam, double *Q, int *mat)
{
  int i, j, id, nct = 0, n = go->n;
  double dx[D], gam2;

  if ( gam < 0 ) {
    gam = 1.2;
  }
  gam2 = gam * gam;

  if ( mat ) {
    for (id = 0; id < n * n; id++) {
      mat[id] = 0;
    }
  }

  for ( i = 0; i < n - 1; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      id = i * n + j;
      if ( !go->iscont[id] ) { /* skip a noncontact pair */
        continue;
      }
      if ( vsqr( vdiff(dx, x[i], x[j]) ) < go->r2ref[id] * gam2 ) {
        if ( mat ) mat[id] = mat[j*n + i] = 1;
        nct++;
      }
    }
  }

  if ( Q ) *Q = nct / go->ncont;
  return nct;
}



/* write position/velocity file */
__inline static int cago_writepos(cago_t *go, vct *x, vct *v, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  /* print the information line */
  fprintf(fp, "# %d %d\n", n, (v != NULL));

  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%g %g %g", x[i][0], x[i][1], x[i][2]);
    if ( v != NULL ) {
      fprintf(fp, " %g %g %g", v[i][0], v[i][1], v[i][2]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



/* read position/velocity file */
__inline static int cago_readpos(cago_t *go, vct *x, vct *v, const char *fn)
{
  char s[1024];
  FILE *fp;
  int i, n = go->n, hasv = 0, next, ret = -1;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  if ( fgets(s, sizeof s, fp) == NULL || s[0] != '#' ) {
    fprintf(stderr, "%s has no information line\n", fn);
    rewind(fp);
  } else {
    if ( 2 != sscanf(s + 1, "%d%d", &n, &hasv) || n != go->n ) {
      fprintf(stderr, "%s: first line is corrupted:\n%s", fn, s);
      fclose(fp);
      return -1;
    }
  }

  for ( i = 0; i < n; i++ ) {
    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s: cannot read line for i %d\n", fn, i);
      break;
    }

    if ( s[0] == '#' || s[0] == '\n' ) {
      continue;
    }

    if ( 3 != sscanf(s, "%lf%lf%lf%n", &x[i][0], &x[i][1], &x[i][2], &next) ) {
      fprintf(stderr, "%s: cannot read position for i %d\n", fn, i);
      break;
    }

    if ( hasv && 3 != sscanf(s + next, "%lf%lf%lf", &v[i][0], &v[i][1], &v[i][2]) ) {
      fprintf(stderr, "%s: cannot read velocity for i %d\n", fn, i);
      break;
    }
  }
  if ( i >= n ) {
    ret = 0;
  }

  fclose(fp);
  return ret;
}



/* output pdb format */
__inline static int cago_writepdb(cago_t *go, vct *x, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for (i = 0; i < n; i++)
    fprintf(fp, "ATOM  %5d  CA  %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n",
        i + 1, aanames[ go->iaa[i] ], i + 1, x[i][0], x[i][1], x[i][2]);
  fprintf(fp, "END%77s\n", " ");
  fclose(fp);
  return 0;
}



#endif /* CAGOCORE_H__ */


/* obtain helix parameters by a structural fitting */
#include "util.h"
#include "mat.h"


char *fninp = "alpha_n_22.pdb";
int itmax = 10;



/* data structure for the coordinates */



/* atoms types that we care about */
enum { N, HN, CA, HA, CB, C, O, ATOM_TYPES };

const char *atom_types[] = {"N", "HN", "CA", "HA", "CB", "C", "O"};



/* parse the input file into coordinates
 * return the number of residues
 * The array `x` should be able to hold nresmax * ATOM_TYPES atoms
 * */
static int parse(const char *fn, double (*x)[3], int nresmax)
{
  FILE *fp;
  char s[512], atnm[8], tok[16];
  int i, atp, resid, resmin = 10000, resmax = -1, nres;
  double y[3];

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  /* first pass of the scan, get the number of residues
   * and the offset */
  while ( fgets(s, sizeof s, fp) ) {
    if ( strncmp(s, "ATOM  ", 6) != 0 ) {
      continue;
    }

    /* grab the residue number */
    strncpy(tok, s + 22, 4);
    tok[4] = '\0';
    resid = atoi( strstrip(tok) );

    if ( resid < resmin ) {
      resmin = resid;
    }

    if ( resid > resmax ) {
      resmax = resid;
    }
  }

  nres = resmax - resmin + 1;

  /* clean up the coordinates */
  for ( i = 0; i < ATOM_TYPES * nresmax; i++ ) {
    vzero( x[i] );
  }

  /* go back to the beginning of the file */
  rewind(fp);

  /* second scan get the atom coordinates */
  while ( fgets(s, sizeof s, fp) ) {
    if ( strncmp(s, "ATOM  ", 6) != 0 ) {
      continue;
    }

    /* grab the atom name */
    strncpy(atnm, s + 12, 4);
    atnm[4] = '\0';
    strstrip(atnm);

    /* compare the residue name from the requested atom list */
    for ( atp = 0; atp < ATOM_TYPES; atp++ ) {
      if ( strcmp(atnm, atom_types[atp]) == 0 ) {
        break;
      }
    }
    if ( atp >= ATOM_TYPES ) { /* not a wanted atom */
      continue;
    }

    /* grab the residue number */
    strncpy(tok, s + 22, 4);
    tok[4] = '\0';
    resid = atoi( strstrip(tok) ) - resmin;

    if ( resid >= nresmax ) {
      fprintf(stderr, "too many residues %d >= %d\n", resid, nresmax);
      break;
    }

    /* grab the coordinates */
    for ( i = 0; i < 3; i++ ) {
      strncpy(tok, s + 30 + 8 * i, 8);
      tok[8] = '\0';
      sscanf(tok, "%lf", &y[i]);
      x[resid * ATOM_TYPES + atp][i] = y[i];
    }

    //fprintf(stderr, "%d, %d, %lf %lf %lf\n", resid, atp, y[0], y[1], y[2]);
  }

  fclose(fp);

  return nres;
}



static int savepdb(double (*x)[3], int nres, const char *fn)
{
  FILE *fp;
  int i, j, atp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  for ( i = 0; i < nres; i++ ) {
    for ( atp = 0; atp < ATOM_TYPES; atp++ ) {
      j = i * ATOM_TYPES + atp;
      fprintf(fp, "ATOM  %5d  %-3s ALA  %4d    %8.3f%8.3f%8.3f  1.00  1.00           %c  \n",
          j + 1, atom_types[atp], i + 1, x[j][0], x[j][1], x[j][2], atom_types[atp][0]);
    }
  }

  fclose(fp);

  return 0;
}



typedef struct {
  double r;    /* radius from the helical axis */
  double dz;   /* rise along the helical axis */
  double z0;   /* initial rise */
  double dth;  /* angle per residue */
  double th0;  /* initial angle */
} params_t;



/* build model from the parameters */
static int build_model(double (*xm)[3], int nres,
    const params_t *params)
{
  int i, atp, j;
  double th;

  for ( j = 0, i = 0; i < nres; i++ ) {
    for ( atp = 0; atp < ATOM_TYPES; atp++, j++ ) {
      th = params[atp].dth * i + params[atp].th0;
      xm[j][0] = params[atp].r * cos(th);
      xm[j][1] = params[atp].r * sin(th);
      xm[j][2] = params[atp].dz * i + params[atp].z0;
    }
  }
  return 0;
}



static int improve_model(double (*xm)[3], double (*xf)[3],
    int n, params_t *params, int updatedth)
{
  int i, atp, j;
  double r, th, c, s, m, mm;
  params_t np[ATOM_TYPES] = {{0}}; /* numerators */
  params_t dp[ATOM_TYPES] = {{0}}; /* denominators */

  for ( j = 0, i = 0; i < n; i++ ) {
    for ( atp = 0; atp < ATOM_TYPES; atp++, j++ ) {
      /* for the radius */
      r = sqrt( xf[j][0] * xf[j][0] + xf[j][1] * xf[j][1] );
      np[atp].r += r;

      /* for the z parameters */
      np[atp].z0 += xf[j][2];
      np[atp].dz += i * xf[j][2];

      /* for the theta parameters */
      th = params[atp].dth * i + params[atp].th0;
      c = cos(th);
      s = sin(th);
      m = -xf[j][0] * r * s + xf[j][1] * r * c;
      mm = r * r; /* more stable than  xf[j][0] * r * c + xf[j][1] * r * s; */
      np[atp].th0 += m;
      dp[atp].th0 += mm;

      /* the updates to dth can be aggregated */
      np[0].dth += m * i;
      dp[0].dth += mm * i;
    }
  }

  /* update the new parameters */
  if ( updatedth ) { /* the dth parameters is very vulnerable */
    params[0].dth += 0.5 * np[0].dth / dp[0].dth;
    printf("dth %g\n", np[0].dth / dp[0].dth);
  }
  for ( atp = 0; atp < ATOM_TYPES; atp++ ) {
    printf("Before: %-2s: r %8.5f, z0 %8.5f, dz %8.5f, th0 %8.5f(%8.3f), dth %8.5f(%8.3f)\n",
        atom_types[atp], params[atp].r, params[atp].z0, params[atp].dz,
        params[atp].th0, params[atp].th0 * 180 / M_PI,
        params[atp].dth, params[atp].dth * 180 / M_PI);

    params[atp].r = np[atp].r / n;

    /* np.z0 = Sum z_i
     * np.dz = Sum i z_i
     * */
    params[atp].z0 = ((4 * n - 2) * np[atp].z0 - 6 * np[atp].dz)
                   / (n*n + n);
    params[atp].dz = 12 * (np[atp].dz - 0.5 * (n - 1) * np[atp].z0)
                   / (n*n*n - n);

    params[atp].th0 += np[atp].th0 / dp[atp].th0;
    params[atp].dth = params[0].dth;

    printf("After:  %-2s: r %8.5f, z0 %8.5f, dz %8.5f, th0 %8.5f(%8.3f), dth %8.5f(%8.3f)\n",
        atom_types[atp], params[atp].r, params[atp].z0, params[atp].dz,
        params[atp].th0, params[atp].th0 * 180 / M_PI,
        params[atp].dth, params[atp].dth * 180 / M_PI);
  }

  build_model(xm, n, params);

  return 0;
}



#define NRESMAX 100

static int dofit(void)
{
  int it;
  int nres, natm;
  static double x[NRESMAX * ATOM_TYPES][3];
  static double xf[NRESMAX * ATOM_TYPES][3];
  static double xm[NRESMAX * ATOM_TYPES][3];
  double R[D][D], T[D], rmsd;
  params_t params[ATOM_TYPES] = {
    /* r, dz, z0, dth, th0 */
    {2.5, 1.5, 0.0, M_PI * 100 / 180, M_PI *  3 / 180},
    {2.5, 1.5, 0.0, M_PI * 100 / 180, M_PI *  3 / 180},
    {2.5, 1.5, 0.3, M_PI * 100 / 180, M_PI * 30 / 180},
    {2.5, 1.5, 0.3, M_PI * 100 / 180, M_PI * 30 / 180},
    {2.5, 1.5, 0.6, M_PI * 100 / 180, M_PI * 45 / 180},
    {2.5, 1.5, 0.9, M_PI * 100 / 180, M_PI * 60 / 180},
    {2.5, 1.5, 1.2, M_PI * 100 / 180, M_PI * 55 / 180}
  };

  /* parse the input PDB file */
  nres = parse(fninp, x, NRESMAX);

  /* set the initial parameters */
  build_model(xm, nres, params);

  natm = nres * ATOM_TYPES;

  rmsd = vrmsd(x, xf, xm, NULL, natm, 0, R, T);

  /* iteratively improve the parameters */
  for ( it = 0;  it < itmax; it++ ) {
    /* improve the angular parameters */
    improve_model(xm, xf, nres, params, 0);

    /* improve the global rotation matrix and translation vector
     * hence the fitting structure */
    rmsd = vrmsd(x, xf, xm, NULL, natm, 0, R, T);
    printf("t %d, rmsd %g\n", it, rmsd);
  }

  /* save the model structure */
  savepdb(xm, nres, "out.pdb");

  return 0;
}



int main(int argc, char **argv)
{
  if ( argc > 1 ) fninp = argv[1];
  dofit();
  return 0;
}



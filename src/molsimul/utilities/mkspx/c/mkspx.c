/*
  Building an extended chain (beta-strand, spiral-like) conformation

  Copyright (C) 2010-2013  Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  gcc -O3 mkspx.c -lm

  TODO:
  *  more advanced rotamer specification
  *  collision detection
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

int verbose = 0;
const char *prog = "mkspx"; /* name of the program */
const char *fnout = "out.pdb";
const char *fnin  = "in.seq";
static void help(void);

/* standard bond lengths in angstroms */
#define B_CAC     1.53
#define B_CN_PEP  1.33
#define B_NCA     1.46
#define B_CC      1.54
#define B_CN      1.48 /* peptide bond */
#define B_CO      1.22 /* carbonyl */
#define B_CS      1.81
#define BH_OHOH   2.80
#define BH_NHOH   2.90
#define BH_OHOC   2.80
#define B_CC_RING 1.40

#define D2R (M_PI/180.0)
#define R2D (180.0/M_PI)

/* standard cosine and sine values */
double c36, s36, c12, s12, c72, s72, c30, s30;

/* rotating angle in the x-y plane around z-axis (z parallel to C=O or NH) */
double rotang = 10.0*D2R;
/* swinging e angle between successive CO-NH peptide planes */
double swgang = 50.0*D2R;
/* rising angle along z-axis */
double risang = 9.0*D2R;

#define NMAX 510  /* maximal # of residues */
#define NM (NMAX+2)
int nres;
char name[NM][8] = {"ACE"};
char letterseq[NM]="";
int pgamma[NM] = {0}; /* C-gamma rotamer position */

static double *mkv(double v[], double x, double y, double z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
  return v;
}

static double *neg(double a[3])
{
  a[0] = -a[0];
  a[1] = -a[1];
  a[2] = -a[2];
  return a;
}

static double *copy(double b[3], double a[3])
{
  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  return b;
}

/* c = a*gam1+b*gam2 */
static double *lincomb2(double c[3], double a[3], double b[3], double gam1, double gam2)
{
  c[0] = gam1*a[0] + gam2*b[0];
  c[1] = gam1*a[1] + gam2*b[1];
  c[2] = gam1*a[2] + gam2*b[2];
  return c;
}

#define sadd(c, a, b, l) lincomb2(c, a, b, 1.0, l)

#define add(c, a, b) sadd(c, a, b, 1.0)

#define diff(c, a, b) sadd(c, a, b, -1.0)

static double *cross(double q[3], double u[3], double v[3])
{
  q[0] = u[1]*v[2]-u[2]*v[1];
  q[1] = u[2]*v[0]-u[0]*v[2];
  q[2] = u[0]*v[1]-u[1]*v[0];
  return q;
}

static double dot(double a[3], double b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static double norm(double v[3]) { return sqrt(dot(v, v)); }

static double *normalize(double v[3])
{
  double x = 1.0/norm(v);
  v[0] *= x; v[1] *= x; v[2] *= x;
  return v;
}

/* normalized component of a that is perpendicular to b */
static double *perpen(double p[], double a[], double b[])
{
  double d = dot(a, b);
  return normalize( sadd(p, a, b, -d/dot(b, b)) );
}


enum AATYPE {ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE,
  LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL};
char aaletter[64]="ARNDCEQGHILKMFPSTWYV";

/* return the integer index */
static int name2itype(const char *name)
{
  const char *p = name;

  /* be careful with 4-letter residues */
  if (strcmp(name, "CYS2") == 0) {
    return CYS;
  } else if (name[3] != '\0' && (name[0]=='C'||name[0]=='N'))
    p++; /* skip the leading terminal letter 'C' or 'N' */

  if      (strcmp (p, "ALA") == 0)    return ALA;
  else if (strcmp (p, "ARG") == 0)    return ARG;
  else if (strcmp (p, "ASN") == 0)    return ASN;
  else if (strcmp (p, "ASP") == 0)    return ASP;
  else if (strncmp(p, "CY", 2) == 0)  return CYS;
  else if (strcmp (p, "GLU") == 0)    return GLU;
  else if (strcmp (p, "GLN") == 0)    return GLN;
  else if (strcmp (p, "GLY") == 0)    return GLY;
  else if (strncmp(p, "HI", 2) == 0)  return HIS;
  else if (strcmp (p, "ILE") == 0)    return ILE;
  else if (strcmp (p, "LEU") == 0)    return LEU;
  else if (strncmp(p, "LY", 2) == 0)  return LYS;
  else if (strcmp (p, "MET") == 0)    return MET;
  else if (strcmp (p, "PHE") == 0)    return PHE;
  else if (strcmp (p, "PRO") == 0)    return PRO;
  else if (strcmp (p, "SER") == 0)    return SER;
  else if (strcmp (p, "THR") == 0)    return THR;
  else if (strcmp (p, "TRP") == 0)    return TRP;
  else if (strcmp (p, "TYR") == 0)    return TYR;
  else if (strcmp (p, "VAL") == 0)    return VAL;
  return -1;
}

#define c_terminal(name) check_terminal(name, 'C')
#define n_terminal(name) check_terminal(name, 'N')

static int check_terminal(char *name, char ch)
{
  if (strcmp(name, "CYS2") == 0) return 0;
  if (name[3] !='\0' && name[0] == ch) return 1;
  return 0;
}

double max_x, max_y, max_z;
double min_x, min_y, min_z;
FILE *fpout;

/* print a line for an atom to pdbfile */
static void atomput(int atomid, char *atomname, char*resname, int resid,
    double rvec[3], char *ele) {
  static char myele[4]="C";
  if (ele) {
    if (ele[0] != atomname[0]) {
      fprintf(stderr, "mismatch: atom=%s, ele=%s\n", atomname, ele);
      exit(1);
    }
    assert(ele[0] == atomname[0]);
  }
  if (ele) myele[0] = ele[0]; else myele[0] = atomname[0];
  if (rvec[0] > max_x) max_x = rvec[0]; else if (rvec[0] < min_x) min_x = rvec[0];
  if (rvec[1] > max_y) max_y = rvec[1]; else if (rvec[1] < min_y) min_y = rvec[1];
  if (rvec[2] > max_z) max_z = rvec[2]; else if (rvec[2] < min_z) min_z = rvec[2];
  fprintf(fpout, "ATOM   %4d  %-3s %-4s %4d    %8.3f%8.3f%8.3f  1.00  1.00          %s\n",
    atomid, atomname, resname, resid, rvec[0], rvec[1], rvec[2], myele);
}


/* fill trp coordinates, g: coords of CG; x and y: unit vectors of trp plane */
static void gentrp(double r[8][3], double g[3], double x[3], double y[3])
{
  double v[3];
  sadd(r[0], g,    lincomb2(v, x, y,  c36,  s36), B_CC_RING); /* CD1 */
  sadd(r[1], g,    lincomb2(v, x, y, -c36,  s36), B_CC_RING); /* CD2 */
  sadd(r[2], r[0], lincomb2(v, x, y, -c72,  s72), B_CC_RING); /* NE1 */
  sadd(r[3], r[1], lincomb2(v, x, y,  c72,  s72), B_CC_RING); /* CE2 */
  sadd(r[4], r[1], lincomb2(v, x, y, -c12, -s12), B_CC_RING); /* CE3 */
  sadd(r[5], r[3], lincomb2(v, x, y, -c72,  s72), B_CC_RING); /* CZ2 */
  sadd(r[7], r[5], lincomb2(v, x, y, -c12, -s12), B_CC_RING); /* CH2 */
  sadd(r[6], r[7], lincomb2(v, x, y, -c72, -s72), B_CC_RING); /* CZ3 */
}

/* remove decorations around the residue name */
static char *restrim(char *resname)
{
  char *p, name2[8];
  char filter[] = "\"\',:;./\\~()[]{}";

  p = resname + strlen(resname) - 1;
  while (p >= resname && strchr(filter, *p) != NULL)
    *p--='\0';
  if (resname[0] == '\0') return resname;
  for (p = resname; strchr(filter, *p); p++)
    ;
  strcpy(name2, p);
  return strcpy(resname, name2);
}

/* read the sequence file */
static int readseq(const char *fname)
{
  FILE *fp;
  static char buf[1024], token[16];
  int i, j;

  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", fname);
    help();
  }

  if (fgets(buf, sizeof buf, fp) == NULL) { /* first line */
    fprintf(stderr, "unable to get the first line of %s\n", fname);
    return -1;
  }

  if (buf[0] == '#') { /* read parameters */
    double a1, a2, a3;
    if (3 <= sscanf(buf+1, "%lf%lf%lf%s", &a1, &a2, &a3, letterseq)) {
      rotang = a1 * D2R;
      swgang = a2 * D2R;
      risang = a3 * D2R;
    } else {
      fprintf(stderr, "the parameter lines is corrupted.\n%s", buf+1);
    }
  } else {
    rewind(fp); /* back to the beginning */
  }

  nres = 0;
  /* read residue sequence */
  for (i = 1; i <= NMAX; i++) {
    memset(token, '\0', sizeof token);
    if (1 != fscanf(fp, "%7s", token)) {
      nres = i - 1;
      goto CLOSE_FILE;
    }
    if (token[0] == '#') {
      nres = i - 1;
      goto ROTAMER;
    } else if (token[0] == '\0' || !isalpha(token[0]) ||
        token[1] == '\0' || (strcmp(token, "SEQRES") == 0) ) {
      if (verbose)
        fprintf(stderr, "%-3d: ignore an illegal residue name [%s]\n", i, token);
      i--;
      continue;
    } else if (strcmp(token, "ACE") == 0) {
      i = 0;
      nres = 0;
      continue;
    } else if (strcmp(token, "NH2") == 0) {
      nres = i - 1;
      /* skip the rest of the file */
      while (1 == fscanf(fp, "%7s", token)) {
        if (token[0] == '#') goto ROTAMER;
      }
      goto CLOSE_FILE;
    }
    restrim(strcpy(name[i], token));
  }
  if (i > NMAX) {
    fprintf(fp, "too many residues!\n");
    fclose(fp);
    return 1;
  }

ROTAMER:
  for (j = 0;  fgets(buf, sizeof buf, fp);  j++) {
    int res;
    char rota[32]="";
    if (j == 0) {
      /* previous token something like #93 G3 */
      if (token[1] != '\0') {
        sscanf(token+1, "%d", &res);
        sscanf(buf, "%s", rota);
      } else {
        sscanf(buf, "%d%s", &res, rota);
      }
    } else {
      if (buf[0] != '#') continue;
      sscanf(buf+1, "%d%s", &res, rota);
    }
    if (rota[0] != 'G' || strchr("123", rota[1]) == NULL ) {
      fprintf(stderr, "Cannot understand the rotamer code [%s]!\n", rota);
      continue;
    }
    if (res <= 0 || res > nres) {
      fprintf(stderr, "the residue # %d out of range (1,%d)\n", res, nres);
      continue;
    }
    pgamma[res] = rota[1] - '0';
    printf("Gamma rotamer of residue %s%d is set to %d\n", name[res], res, pgamma[res]);
  }

CLOSE_FILE:
  fclose(fp);
  printf("I got %d residues, rotation %g, swinging %g, rising %g (degrees)\n",
      nres, rotang*R2D, swgang*R2D, risang*R2D);
  for (i = 1; i <= nres; i++) {
    printf("%-4s ", name[i]);
    if (i%10 == 0) printf("\n");
  }
  if (nres%10 != 0) printf("\n");

  /* also print out the one-letter sequence for comparison */
  for (i = 1; i <= nres; i++) {
    int type = name2itype(name[i]);
    if (letterseq[0]) {
      if (letterseq[i-1] != aaletter[type]) {
        printf("\n");
        fprintf(stderr, "The sequence is corrupted at %d: should be %c, I got %s (%c)\n",
            i, letterseq[i-1], name[i], aaletter[type]);
        return 1;
      }
    }
    printf("%c", aaletter[type]);
  }
  printf("\n");
  if (letterseq[0]) {
    fprintf(stderr, "the sequence is verified.\n");
  }
  return 0;
}

/* print a help message and die */
void help(void)
{
  static char options[]=
  "OPTIONS:\n"
  " -h: print this message\n"
  " -o: followed by output pdb file\n"
  " -i: followed by input  sequence file\n"
  " -v: verbose\n";
  static char desc[]=
  "  I read a .seq file, which specifies the 3- or 4-letter amino acid sequence,\n"
  "  such as NMET ALA LEU .... Use space to separate residues, \n"
  "  but you add some funny character like comma or quotes, such as\n"
  "    \"NMET\", [ALA], (LEU), ...\n"
  "  The entire SEQRES section of a PDB file would also work.\n"
  "  For terminal caps, use NMET instead of MET to avoid an additional N-cap,\n"
  "  the same convention goes to C-cap.\n\n"
  "  There can be an optional first line, which looks like\n"
  "  # rotang swgang risang letterseq\n"
  "  - all angles are in degrees\n"
  "  - rotang is the angle of rotation in the x-y plane\n"
  "  - swgang is the angle of swinging between two successive amino acids\n"
  "  - risang is the angle of rising along the z-axis\n"
  "  - letterseq is the one-letter amino acid sequence (no spaces),"
  "    but only for checking purpose\n\n"
  "  After the sequence, optional lines can be added to specify rotamer positions, e.g.,\n"
  "    # 93 G3\n"
  "  or\n"
  "    #93 G3\n"
  "  specifies that C-gamma of the 93rd residue will placed at gamma 3 position.\n"
  "  G1 opposes CO; G2 opposes N; G3 is the remaining direction.\n"
  "  The output is in AMBER format, use amberize -u to convert it to a normal PDB.\n";

  fprintf(stderr, "%s  Copyright (C) 2010-2013 Cheng Zhang\n", prog);
  fprintf(stderr, "USAGE:\n");
  fprintf(stderr, "%s [OPTIONS] your.seq\n\n", prog);
  fprintf(stderr, "%s\n", options);
  fprintf(stderr, "%s\n", desc);
  exit(1);
}

/* handle command-line arguments */
static int doargs(int argc, char *argv[])
{
  char *next;
  int i, j, ch;

  prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnin = argv[i];
      continue;
    }

    ch = argv[i][1];
    if (strchr("io", ch) != NULL) { /* options that require an argument */
      if (argv[i][2] == '\0') { /* read the next argument */
        if (i == argc-1 || argv[i+1][0] == '-') /* invalid argument */
          help();
        next = argv[i+1];
        i++; /* skip the next argument (filename) */
      } else {
        next = argv[i]+2;
      }
      switch (ch) {
      case 'o': fnout = next; break;
      case 'i': fnin  = next; break;
      default: fprintf(stderr, "unexpected option %s\n", argv[i]); exit(1);
      }
      continue;
    }

    for (j = 1; (ch = argv[i][j]) != '\0'; j++)
      switch (ch) {
      case 'v': verbose = 1; break;
      case 'h': default: help();
      }
  }
  return 0;
}



int main(int argc, char *argv[])
{
  int i, k, sgn;
  double xc[3], xca[3], xn[3], xo[3], xcb[3], xnp[3];
  double dir_nca[3] = {0,0,0}, dir_cac[3];
  double xog1[3], xog2[3], xog3[3], xcg1[3], xcg2[3], xcg3[3], *xg; /* GAMMA position */
  double xcd[3]; /* DELTA position */
  double u[3] = {0}, v[3] = {0}, w[3] = {0}, p[3] = {0}, q[3] = {0};
  double os[3];
  double c1, s1, c2, s2, c3, s3, cr, sr;
  double c1p, s1p, c2p, s2p, c1m, s1m, c2m, s2m;
  double thp, thm, ang, cp, sp;

  int resid, atomid;
  int itype;

  doargs(argc, argv);
  readseq(fnin);

  /* calculate constants */
  c30 = sqrt(.75); /* cos(M_PI/6) */
  s30 = 0.5; /* sin(M_PI/6) */
  c36 = cos(36.0*D2R);
  s36 = sin(36.0*D2R);
  c12 = cos(12.0*D2R);
  s12 = sin(12.0*D2R);
  c72 = cos(72.0*D2R);
  s72 = sin(72.0*D2R);

  sr = sin(risang);
  cr = cos(risang);

  /* compute the average angle between N-CA with the horizontal plane
   * The calculation uses the approximation of `rotang = 0'
   * it assumes that the angle N-CA-C is 109.28
   * the bond lengths are irrelevant. Let a = swgang, b = risang,
   * u and v = angles of N-CA and CA-C with the horizontal plane
   * so u = ang + b, v = ang - b
   * For convenience, assuming the distances |N-CA| = |CA-C| = 1
   *   (if not true, relocate N and C at the respective line segments),
   * we want |N-C| = sqrt(8/3) to make cos(N-CA-C) = -1/3. But
   * |N-C|^2 = [cos^2 u + cos^2 v + 2 cos u cos v cos a] + (sin u - sin v)^2
   *         = 2 + cos(u+v)(1 + cos a) + cos(u-v)(cos a - 1)
   * By u + v = 2 ang, and u - v = 2 b, `ang' can be solved from
   *  sin^2 ang = [sin^2 b cos a + sin^2 b - 1/3]/(1 + cos a) */
  ang = cos(swgang) * (cr*cr) + sr*sr - 1.0/3;
  ang /= 1 + cos(swgang);
  printf("rotating angle %g, swinging angle %g, rising angle %g, cos %g\n",
      rotang, swgang, risang, ang);
  if (ang < 0.0) {
    printf("turning angle is too large\n");
    return 1;
  }
  ang = asin(sqrt(ang));

  thp = ang + risang; /* the angle for even residues */
  c1p = cos(thp); /* for CA-C and N-CA */
  s1p = sin(thp);
  c2p = cos(thp - M_PI/3); /* for C-N */
  s2p = sin(thp - M_PI/3);

  thm = risang - ang; /* the angle for odd residues */
  c1m = cos(thm); /* for CA-C and N-CA */
  s1m = sin(thm);
  c2m = cos(thm + M_PI/3); /* for C-N */
  s2m = sin(thm + M_PI/3);
  printf("ang %g, %g, %g (degrees)\n", ang*R2D, thp*R2D, thm*R2D);

  ang = 0.5*acos(-1.0/3);
  c3 = cos(ang);
  s3 = sin(ang);

  if ((fpout = fopen(fnout, "w")) == NULL) {
    fprintf(stderr, "cannot open output.\n");
    return 1;
  }

  resid = atomid = 1;
  os[0] = os[1] = os[2] = 0;
  ang = 0;

  /* loop over residues one by one
   * we always start from i=0,
   * even if we don't want the N-cap ACE, we want to
   * calculate the position of the first N */
  for (i = 0; i <= nres; i++) {
    sgn = (i % 2) ? -1 : 1;
    if (sgn > 0) { /* `+' sign for even-index residues */
      c1 = c1p;
      s1 = s1p;
      c2 = c2p;
      s2 = s2p;
    } else { /* `-' sign for odd-index residues */
      c1 = c1m;
      s1 = s1m;
      c2 = c2m;
      s2 = s2m;
    }
    ang += rotang;
    cp = cos(ang); /* planar rotation */
    sp = sin(ang);

    /* OS is the current backbone atom */
    copy(xca, os); /* CA */
    mkv(dir_cac, c1 * cp, c1 * sp, s1); /* direction from CA to C, */
    /* s1 is the vertical component, c1 is the horizontal component,
     * which is subject to the planar angle cp and sp */

    sadd(xc, xca, dir_cac, B_CAC); /* C */
    sadd(xo, xc, mkv(v, -sr * cp, -sr * sp, cr), B_CO * sgn); /* O */

    /* Warning CB and CG atoms depends on dir_nca, which is unavailable
     * for the 0th atoms. But 0th atom is ACE, so it is OK */
    /* calculate the position of CB */
    normalize( diff(p, dir_nca, dir_cac) ); /* p is the direction from C and N midpoint to CA */
    normalize( cross(q, dir_cac, dir_nca) );
    normalize( lincomb2(w, p, q, c3, s3) ); /* w is the direction from CA to CB */
    sadd(xcb, xca, w, B_CC);

    /* calculate G1, G2 and G3 positions
     * definition:
     * view along the line of CB-CA,
     * G1 is at the opposite of C(O), i.e., CG-CB is along
     * the same direction of CA-C(O)
     * G2 is at the opposite of N (the one bonded to CA)
     * G3 is the rest position
     * if CB is chiral, and the priority in the R-S system is
     *   CA > CG1 > CG2 > CG3,
     * then CB is left-handed (S or L),
     * note: the conclusion is invalid for theonine, in which case
     * the G oxygen, if on G1, has the highest priority! */
    /* CG1-CB should be in the same direction as CA-C(O) */
    normalize( diff(u, xca, xc) );
    sadd(xog1, xcb, u, B_CO);
    sadd(xcg1, xcb, u, B_CC);
    /* OG2-CB should be along the same direction as CA-N
       G2 is opposite to N, viewing along CA-CB */
    normalize( diff(v, xca, xn) );
    sadd(xog2, xcb, v, B_CO);
    sadd(xcg2, xcb, v, B_CC);
    /* G3: opposing H in CA, view along CA-CB */
    perpen(p, u, w); /* w is CA-->CB */
    perpen(q, v, w);
    add(u, p, q); /* add the two perpendicular components */
    normalize( lincomb2(v, u, w, -1.0, 1.0/3) );
    sadd(xog3, xcb, v, B_CO);
    sadd(xcg3, xcb, v, B_CC);

    copy(xnp, xn);
    sadd(xn, xc, mkv(p, c2*cp, c2*sp, s2), B_CN_PEP); /* N */
    /* dir_nca is the direction from N to the next CA
       it is along the same direction from this CA to C */
    copy(dir_nca, dir_cac);
    /* os is the position of the next O */
    sadd(os, xn, dir_nca, B_NCA);

    ang += swgang * sgn; /* so it is flipping around */
    cp = cos(ang);
    sp = sin(ang);

    if (i == 0) {
      if (n_terminal(name[1])) { /* we don't need ACE */
        continue;
      } else { /* we need ACE */
        itype = GLY; /* pretend to be GLY; */
      }
    } else {
      itype = name2itype(name[i]);
      assert(itype >= 0);
    }

    /* start putting atom coordinates */
    if (i > 0) {
      atomput(atomid++, "N", name[i], resid, xnp, "N");
      atomput(atomid++, "CA", name[i], resid, xca, "C");
    } else {
      atomput(atomid++, "CH3", name[i], resid, xca, "C");
    }

    atomput(atomid++, "C", name[i], resid, xc, "C");
    if ( i == nres  &&  c_terminal(name[i]) ) { /* no need for NH2 */
      atomput(atomid++, "OC1", name[i], resid, xo, "O");
      normalize(diff(u, xc, xo));
      normalize(diff(v, xc, xca));
      sadd(p, xc, add(w, u, v), B_CO);
      atomput(atomid++, "OC2", name[i], resid, p, "O");
    } else {
      atomput(atomid++, "O", name[i], resid, xo, "O");
    }

    if (itype == GLY) goto NEXT;

#define GCHOOSE(x, g, g1, g2, g3) { \
    if (pgamma[i] == 1) x = g1; \
    else if (pgamma[i] == 2) x = g2; \
    else if (pgamma[i] == 3) x = g3; \
    else x = g; }

    /* building side chains */
    atomput(atomid++, "CB", name[i], resid, xcb, "C");
    if (itype == ALA) {
      ; /* already done, don't do anything, */
    } else if (itype == SER) {
      GCHOOSE(xg, xog2, xog1, xog2, xog3);
      atomput(atomid++, "OG", name[i], resid, xg, "O");
    } else if (itype == CYS) {
      GCHOOSE(xg, xog2, xog1, xog2, xog3);
      atomput(atomid++, "SG", name[i], resid, xg, "S");
    } else if (itype == VAL) {
      atomput(atomid++, "CG1", name[i], resid, xcg1, "C");
      atomput(atomid++, "CG2", name[i], resid, xcg2, "C");
    } else if (itype == LEU) {
      GCHOOSE(xg, xcg1, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      add(p, xg, diff(q, xcb, xca));
      atomput(atomid++, "CD1", name[i], resid, p, "C");
      diff(u, xg, xcb); /* z axis; */
      normalize(cross(w, q, u)); /* x axis; */
      sadd(v, q,  u, -1.0/3); /* y axis; */
      sadd(p, xg, u,  1.0/3); /* (xcd1 + xcd2 + xcd3)/3 */
      sadd(p, p,  v, -1.0/2); /* (xcd2 + xcd3)/2 */
      sadd(p, p,  w, sqrt(2./3) * B_CC);
      atomput(atomid++, "CD2", name[i], resid, p, "C");
    } else if (itype == ILE) {
      xg  = (pgamma[i] == 3) ? xcg3 : xcg1;
      atomput(atomid++, "CG1", name[i], resid, xg, "C");
      add(p, xg, diff(q, xcb, xca));
      atomput(atomid++, "CD", name[i], resid, p, "C");
      xg  = (pgamma[i] == 3) ? xcg1 : xcg2;
      atomput(atomid++, "CG2", name[i], resid, xg, "C");
    } else if (itype == TRP) {
      static double xtrp[8][3];
      static char trpatomnm[8][4] = {"CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};
      GCHOOSE(xg, xcg1, xcg1, xcg2, xcg3);
      /* calculate two othogonal directions for the trp plane */
      normalize( diff(v, xg, xcb) ); /* v the y direction */
      /* the x direction is define from CA-CB cross CB-CG */
      normalize( cross(u, v, diff(q, xca, xcb) ) ); /* u is the x direction */
      if (xg == xcg2) neg(u); /* flip sign */

      gentrp(xtrp, xg, u, v);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      for (k = 0; k < 8; k++) {
        atomput(atomid++, trpatomnm[k], name[i], resid, xtrp[k], NULL);
      }
    } else if (itype == HIS) {
      GCHOOSE(xg, xcg2, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      normalize(diff(v, xg, xcb)); /* y axis */
      normalize(cross(w, diff(p, xcb, xca), v)); /* x axis */
      lincomb2(q, v, w, s36*B_CC_RING, c36*B_CC_RING);
      add(p, xg, q);
      atomput(atomid++, "ND1", name[i], resid, p, "N");
      lincomb2(q, v, w, s72*B_CC_RING, -c72*B_CC_RING);
      add(u, p, q);
      atomput(atomid++, "CE1", name[i], resid, u, "C");
      lincomb2(q, v, w, s36*B_CC_RING, -c36*B_CC_RING);
      add(p, xg, q);
      atomput(atomid++, "CD2", name[i], resid, p, "C");
      lincomb2(q, v, w, s72*B_CC_RING, c72*B_CC_RING);
      add(u, p, q);
      atomput(atomid++, "NE2", name[i], resid, u, "N");
    } else if (itype == PHE || itype == TYR) {
      GCHOOSE(xg, xcg2, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      normalize(diff(v, xg, xcb)); /* y axis */
      normalize(cross(w, diff(p, xcb, xca), v)); /* x axis */
      lincomb2(q, v, w, s30*B_CC_RING, c30*B_CC_RING);
      add(p, xg, q);
      atomput(atomid++, "CD1", name[i], resid, p, "C");
      sadd(q, p, v, B_CC_RING);
      atomput(atomid++, "CE1", name[i], resid, q, "C");
      sadd(q, xg, v, 2*B_CC_RING);
      atomput(atomid++, "CZ", name[i], resid, q, "C");
      if (itype == TYR) {
        sadd(u, q, v, B_CO);
        atomput(atomid++, "OH", name[i], resid, u, "O");
      }
      lincomb2(q, v, w, s30*B_CC_RING, -c30*B_CC_RING);
      add(p, xg, q);
      atomput(atomid++, "CD2", name[i], resid, p, "C");
      sadd(q, p, v, B_CC_RING);
      atomput(atomid++, "CE2", name[i], resid, q, "C");

    } else if (itype == THR) {
      /* make sure CB is right-handed */
      atomput(atomid++, "OG1", name[i], resid, xog1, "O");
      atomput(atomid++, "CG2", name[i], resid, xcg2, "C");

    } else if (itype == GLU || itype == GLN) {
      GCHOOSE(xg, xcg2, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      add(xcd, xg, diff(v, xcb, xca) );
      atomput(atomid++, "CD", name[i], resid, xcd, "C");
      normalize( v ); /* v is y */
      normalize( cross(u, v, diff(q, xcb, xg) ) ); /* u is x */
      /* NOTE I'm not sure which one should be called as OE1 */
      sadd(p, xcd, lincomb2(q, u, v, c30, s30), B_CO);
      atomput(atomid++, "OE1", name[i], resid, p, "O");
      sadd(p, xcd, lincomb2(q, u, v, -c30, s30), B_CO);
      if (itype == GLU) {
        atomput(atomid++, "OE2", name[i], resid, p, "O");
      } else {
        atomput(atomid++, "NE2", name[i], resid, p, "N");
      }
    } else if (itype == ASP || itype == ASN) {
      GCHOOSE(xg, xcg1, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      normalize( diff(v, xg, xcb) ); /* v is y axis */
      perpen(u, diff(q, xca, xcb), v); /* u is x axis */
      sadd(p, xg, lincomb2(q, u, v, c30, s30), B_CO);
      atomput(atomid++, "OD1", name[i], resid, p, "O");
      sadd(p, xg, lincomb2(q, u, v, -c30, s30), B_CO);
      if (itype == ASP) {
        atomput(atomid++, "OD2", name[i], resid, p, "O");
      } else {
        atomput(atomid++, "ND2", name[i], resid, p, "N");
      }
    } else if (itype == LYS) {
      GCHOOSE(xg, xcg2, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      add(p, xg, diff(q, xcb, xca));
      atomput(atomid++, "CD", name[i], resid, p, "C");
      add(u, p, diff(q, xg, xcb));
      atomput(atomid++, "CE", name[i], resid, u, "C");
      sadd(v, u, normalize(diff(q, xcb, xca)), B_CN);
      atomput(atomid++, "NZ", name[i], resid, v, "N");
    } else if (itype == ARG) {
      GCHOOSE(xg, xcg2, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      add(p, xg, diff(q, xcb, xca));
      atomput(atomid++, "CD", name[i], resid, p, "C");
      sadd(u, p, normalize(diff(v, xg, xcb)), B_CN);
      atomput(atomid++, "NE", name[i], resid, u, "N");
      normalize(cross(w, v, q)); /* w is x axis, v is y axis */
      add(u, u, lincomb2(p, w, v, c30*B_CN, s30*B_CN));
      atomput(atomid++, "CZ", name[i], resid, u, "C");
      add(q, u, lincomb2(p, w, v, c30*B_CN, -s30*B_CN));
      atomput(atomid++, "NH1", name[i], resid, q, "N");
      sadd(p, u, v, B_CN);
      atomput(atomid++, "NH2", name[i], resid, p, "N");
    } else if (itype == MET) {
      GCHOOSE(xg, xcg1, xcg1, xcg2, xcg3);
      atomput(atomid++, "CG", name[i], resid, xg, "C");
      /* the position is not right,
       * SN bond is significantly longer than CN bond,
       * so the angle is wrong, but we don't take that into account */
      sadd(p, xg, normalize(diff(q, xcb, xca)), B_CS);
      atomput(atomid++, "SD", name[i], resid, p, "S");
      sadd(u, p,  normalize(diff(q, xg,  xcb)), B_CS);
      atomput(atomid++, "CE", name[i], resid, u, "C");
    } else if (itype == PRO) {
      normalize(diff(u, xcb, xnp)); /* x axis */
      diff(p, xcb, xca);
      diff(q, xnp, xca);
      normalize(add(v, p, q)); /* y axis */
      lincomb2(p, u, v, -c72*B_CC, s72*B_CC);
      add(q, p, xcb);
      atomput(atomid++, "CD", name[i], resid, q, "C");
      lincomb2(p, u, v, c72*B_CC, s72*B_CC);
      add(q, p, xnp);
      atomput(atomid++, "CG", name[i], resid, q, "C");
    } else {
      fprintf(stderr, "%d: %s not supported.\n", resid, name[i]);
    }

NEXT:
    resid++;
  } /* end of the loop over residues */

  if ( !c_terminal(name[nres]) ) {
    nres++;
    strcpy(name[nres], "NH2");
    atomput(atomid++, "N", name[nres], resid++, xn, "N");
  }

  fprintf(fpout, "TER    %4d      %-4s %4d%56s\n",
    atomid, name[nres], resid - 1, " ");

  fclose(fpout);

  xc[0] = max_x-min_x;
  xc[1] = max_y-min_y;
  xc[2] = max_z-min_z;
  fprintf(stderr,
      "x: min %g, max %g, delta=%g\n"
      "y: min %g, max %g, delta=%g\n"
      "z: min %g, max %g, delta=%g\n",
      min_x, max_x, xc[0],
      min_y, max_y, xc[1],
      min_z, max_z, xc[2]);
  if (xc[1] > xc[0]) xc[0] = xc[1];
  if (xc[2] > xc[0]) xc[0] = xc[2];
  fprintf(stderr, "max dimension is %g\n", xc[0]);
  fprintf(stderr, "output saved to %s\n", fnout);
  return 0;
}

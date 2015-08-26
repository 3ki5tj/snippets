double msolvezero_lasty;
double msolvezero_lasttol;

/* Solve matrix equation a x = 0 by Gaussian elimination (full-pivot)
 * The matrix 'a' is destroyed, solutions are saved as *row* vectors in 'x'
 * return the number of solutions */
__inline static int msolvezero(double a[D][D], double (*x)[D], double tol)
{
  double y;
  int i, j, k, cmap[D], sgn = 1;

  for ( i = 0; i < D; i++ ) {
    cmap[i] = i;
  }

  for ( i = 0; i < D; i++ ) {
    /* find the pivot, the largest element in the matrix */
    y = mpivotf_(a, i, cmap, &sgn);
    msolvezero_lasty = y;
    msolvezero_lasttol = tol;
    if ( y <= tol ) { /* we have D - i solutions */
      break;
    }

    /* normalize the row i */
    for ( y = a[i][i], k = i; k < D; k++ ) a[i][k] /= y;

    /* use the pivot to simplify the matrix */
    for ( j = 0; j < D; j++ ) { /* for rows j */
      if ( j != i ) {
        for ( y = a[j][i], k = i; k < D; k++ ) { /* for columns k >= i*/
          a[j][k] -= y * a[i][k];
        }
      }
    }
  }

  /* solve the D - i solutions */
  for ( j = 0; j < D - i; j++ ) {
    vzero( x[j] );
    for ( k = 0; k < i; k++ ) {
      x[j][ cmap[k] ] = -a[k][i + j];
    }
    x[j][ cmap[i + j] ] = 1.0;
    vnormalize( x[j] );
  }

  return D - i;
}



__inline static int meigvecs_low(double (*vecs)[D],
    double mat[D][D], double val, double tol)
{
  double m[D][D];
  int d;

  mcopy(m, mat); /* make a matrix */
  for ( d = 0; d < D; d++ ) {
    m[d][d] -= val;
  }
  return msolvezero(m, vecs, tol);
}



/* maximal acceptable relative tolerance for solving eigenvectors */
double meig_reltol = 1e-6;

#define meigvecs(vecs, mat, val) meigvecs_(vecs, mat, val, meig_reltol)

/* given an eigenvalue, return the corresponding eigenvectors
 * Note: there might be multiple eigenvectors for the eigenvalue */
__inline static int meigvecs_(double (*vecs)[D], double mat[D][D],
    double val, double reltol)
{
  double rtol, max = 0;
  int i = 0, j;

  for ( i = 0; i < D; i++ )
    for ( j = 0; j < D; j++ )
      if ( fabs(mat[i][j]) > max ) max = mat[i][j];

  /* increase the tolerance, until a solution is found */
  for ( i = 0, rtol = DBL_EPSILON * 10; rtol < reltol; rtol *= 10 ) {
    if ( (i = meigvecs_low(vecs, mat, val, max * rtol)) > 0 )
      break;
  }
  return i;
}



/* compute eigenvalues of a 3x3 matrix
 * by solving a cubic equation */
__inline static double *meigval(double v[3], double a[3][3])
{
  double m, p, q, pr, pr3, a00, a11, a22;

  m = (a[0][0] + a[1][1] + a[2][2])/3;
  a00 = a[0][0] - m;
  a11 = a[1][1] - m;
  a22 = a[2][2] - m;
  q = ( a00 * (a11*a22 - a[1][2]*a[2][1])
      + a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a22)
      + a[0][2] * (a[1][0]*a[2][1] - a11*a[2][0]) ) / 2.0;
  p = (a00*a00 + a11*a11 + a22*a22) / 6.0
    + (a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]) / 3.0;
  /* solve x^3 - 3 p x  - 2 q = 0 */
  pr = sqrt(p);
  pr3 = p * pr;
  if ( pr3 <= fabs(q) ) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = m + pr;
      v[2] = m - 2.0 * pr;
    } else { /* phi = 0 */
      v[0] = m + 2.0 * pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    double phi = acos(q/pr3)/3.0; /* 0 < phi < pi/3 */

    v[0] = m + 2.0 * pr * cos(phi);  /* largest */
    v[1] = m + 2.0 * pr * cos(phi - 2*M_PI/3); /* second largest */
    v[2] = m + 2.0 * pr * cos(phi + 2*M_PI/3); /* smallest */
  }
  return v;
}



#define meigsys(v, vecs, mat, nt) meigsys_(v, vecs, mat, nt, meig_reltol)

/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs'
 * ideally, eigenvalues are sorted in descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
__inline static int meigsys_(double v[3], double vecs[3][3], double mat[3][3],
    int nt, double reltol)
{
  double vs[5][3] = {{0}}; /* for safety, vs needs 5 rows */
  int n = 0, nn, i = 0;

  /* eigenvalues are sorted in descending order */
  meigval(v, mat);

  for ( nn = i = 0; i < 3; i++ ) {
    n = meigvecs_(vs + nn, mat, v[nn], reltol);
    if ( n == 0 ) {
      fprintf(stderr, "meigsys failed: try to increase msolvezero_reltol, i %d, nn %d, %g > %g\n",
          i, nn, msolvezero_lasty, msolvezero_lasttol);
      return -1;
    }
    /* if we get multiple 'n' eigenvectors for the same eigenvalue
     * we have to advance the index 'nn' for eigenvalues */
    if ( (nn += n) >= 3 ) break;
  }

  /* maybe we can obtain the last eigenvector by cross-product? */

  mcopy(vecs, vs);
  msort2(v, vecs, NULL);

  if ( !nt ) {
    mtrans(vecs);
  }
  return 0;
}





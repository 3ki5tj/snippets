


"use strict";



var PI = 3.141592653589793;



/* remove the center of mass motion */
function rmcom(x, m, n)
{
  var i;
  var xc = [0,0,0];
  var mtot = 0;

  for ( i = 0; i < n; i++ ) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
function shiftang(x, v, m, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  var mtot = 0;

  for (i = 0; i < n; i++) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);

  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vsinc(am, ang, m[i]);

    xx += m[i] * xi[0] * xi[0];
    yy += m[i] * xi[1] * xi[1];
    zz += m[i] * xi[2] * xi[2];
    xy += m[i] * xi[0] * xi[1];
    yz += m[i] * xi[1] * xi[2];
    zx += m[i] * xi[2] * xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  var inv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  minv(inv, mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  // ang is the solution of M^(-1) * L
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(dv, ang, xi);
    vinc(v[i], dv);
  }
}



var aanames = [
  "GLY", "ALA", "VAL", "LEU", "ILE",
  "PRO", "SER", "THR", "CYS", "MET",
  "ASN", "GLN", "ASP", "GLU",
  "LYS", "ARG", "HIS",
  "PHE", "TYR", "TRP"];

var aaletters = [
  "G", "A", "V", "L", "I",
  "P", "S", "T", "C", "M",
  "N", "Q", "D", "E",
  "K", "R", "H",
  "F", "Y", "W"];

var aacolors = [
  "#808080", "#c0c0c0", "#a0a0c0", "#a0c0c0", "#a0c0a0",
  "#c0a0a0", "#c0a0c0", "#c0c0a0", "#c0c000", "#a0a000",
  "#a040a0", "#e080e0", "#a00000", "#e00000",
  "#0000a0", "#0000e0", "#8080e0",
  "#00a0a0", "#00a000", "#804000"];

var aaradii = [
  0.90, 1.00, 1.10, 1.20, 1.20,
  1.10, 1.20, 1.20, 1.20, 1.40,
  1.20, 1.30, 1.25, 1.35,
  1.35, 1.40, 1.40,
  1.40, 1.40, 1.50];

/* residue name to integer of amino acid */
function res2iaa(res)
{
  var i;

  for ( i = 0; i < aanames.length; i++ ) {
    if ( res === aanames[i] ) {
      return i;
    }
  }
  return -1;
}




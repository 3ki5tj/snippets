/* Three-dimensional Lennard-Jones fluid */



"use strict";



/* initialize a fcc lattice */
function lj_initfcc3d(lj)
{
  var i, j, k, id, n = lj.n;

  var n1 = Math.floor(Math.pow(2*n, 1.0/3) + 0.999999); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++) {
    for (j = 0; j < n1 && id < n; j++) {
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 === 0) {
          /* add some noise to prevent two atoms happened to
           * be separated by precisely some special cutoff distance,
           * which might be half of the box */
          lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][2] = (k + 0.5) * a + noise * (2*rand01() - 1);
          id++;
        }
      }
    }
  }
}



/* get the tail correction */
function lj_gettail3d(lj, rho, n)
{
  var irc, irc3, irc6, utail, ptail;

  irc = 1 / lj.rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = 8 * Math.PI * rho * n / 9 * (irc6 - 3) * irc3;
  ptail = 32 * Math.PI * rho * rho / 9 * (irc6 - 1.5) * irc3;
  return [utail, ptail];
}



/* apply the view matrix */
function transform(x, l)
{
  var n = x.length;
  var xyz = newarr2d(n, 3), xc = [l * 0.5, l * 0.5, l * 0.5], xi = [0, 0, 0];

  for ( var i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vmxv(xyz[i], viewmat, xi);
    vinc(xyz[i], xc);
  }
  return xyz;
}





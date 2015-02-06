/* Two-dimensional Lennard-Jones fluid */



"use strict";



/* initialize a fcc lattice */
function lj_initfcc2d(lj)
{
  var i, j, id = 0, n = lj.n;
  var n1 = Math.floor( Math.sqrt(2*n) + 0.999999 ); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;

  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      if ( (i+j) % 2 === 0 ) {
        // add some noise to prevent two atoms happened to
        // be separated precisely by the cutoff distance,
        // which might be half of the box
        lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
        lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
        id++;
      }
    }
  }
}



/* get the tail correction */
function lj_gettail2d(lj, rho, n)
{
  var irc, irc2, irc6, utail, ptail;

  irc = 1 / lj.rc;
  irc2 = irc * irc;
  irc6 = irc2 * irc2 * irc2;
  utail = Math.PI * rho * n * (0.4*irc6 - 1) * irc2 * irc2;
  ptail = Math.PI * rho * rho * (2.4*irc6 - 3) * irc2 * irc2;
  return [utail, ptail];
}




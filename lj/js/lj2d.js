/* Two-dimensional Lennard-Jones fluid */



"use strict";



/* initialize a fcc lattice */
function lj_initfcc2d(lj)
{
  var i, j, id = 0, n = lj.n;
  var n1 = Math.floor( Math.sqrt(2*n) + .999999 ); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;

  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      if ( (i+j) % 2 == 0 ) {
        // add some noise to prevent two atoms happened to
        // be separated precisely by the cutoff distance,
        // which might be half of the box
        lj.x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
        lj.x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
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



/* annihilate the total angular momentum */
function lj_shiftang2d(x, v, n)
{
  var i;
  var am, r2, xc = [0, 0], xi = [0, 0];

  for (i = 0; i < n; i++) vinc(xc, x[i]);
  vsmul(xc, 1.0 / n);
  for (am = r2 = 0.0, i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    am += vcross2d(xi, v[i]);
    r2 += vsqr(x[i]);
  }
  am = -am / r2;
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am*xi[1];
    v[i][1] +=  am*xi[0];
  }
}



// draw all atoms in the box
function ljdraw2d(lj, target)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // draw the background
  ctx.fillStyle = "#d0d0d0";
  ctx.fillRect(0, 0, width, height);

  // the system dimension is L + two radii
  var scale = Math.min(width, height) / (lj.l + 1);
  var radius = 0.5 * scale;
  var margin = 0.5 * scale;

  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect(margin, margin, width - 2*margin, height - 2*margin);

  // draw each particle
  for (var i = 0; i < lj.n; i++) {
    var x = lj.x[i][0] * scale;
    var y = lj.x[i][1] * scale;
    var spotcolor = "#a0a0e0"
    var color = "#2040a0";
    drawBall(ctx, x + margin, y + margin, radius, color, spotcolor);
  }
}




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



var mousedown = 0;
var mousemoved = 0;
var mousex = -1;
var mousey = -1;
var viewmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];



function ljmousedown(e)
{
  e = e || window.event;
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown = 1;
  //console.log("mousedown", e.clientX, e.clientY, m2str(viewmat));
}



function ljmouseup(e)
{
  e = e || window.event;
  mousex = -1;
  mousey = -1;
  mousemoved = mousedown - 1;
  mousedown = 0;
  //console.log("mouseup", e.clientX, e.clientY, m2str(viewmat));
}



function ljmousemove(e)
{
  if ( !mousedown ) {
    return;
  }
  e = e || window.event;
  if ( mousex >= 0 && mousey >= 0 ) {
    var target = e.target ? e.target : e.srcElement;
    viewmat = mxrot3d(viewmat, 180.0 * (e.clientY - mousey) / target.height);
    viewmat = myrot3d(viewmat, 180.0 * (e.clientX - mousex) / target.width);
    paint(); // defined in ljmain.js
  }
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown += 1;
  //console.log("mousemove", e.clientX, e.clientY, m2str(viewmat));
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





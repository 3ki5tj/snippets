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



/* annihilate the total angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
function lj_shiftang3d(x, v, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) {
    vinc(xc, x[i]);
  }
  vsmul(xc, 1.0/n);
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vinc(am, ang);
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  var inv = rm3_inv(mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  // ang is the solution of M^(-1) * I
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(dv, ang, xi);
    vinc(v[i], dv);
  }
}



var mousedown = false;
var mousex = -1;
var mousey = -1;
var viewmat = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]];



function ljmousedown(e)
{
  e = e || window.event;
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown = true;
  //console.log("mousedown", e.clientX, e.clientY, m2str(viewmat));
}



function ljmouseup(e)
{
  e = e || window.event;
  mousex = -1;
  mousey = -1;
  mousedown = false;
  //console.log("mouseup", e.clientX, e.clientY, m2str(viewmat));
}



function ljmousemove(e)
{
  if ( !mousedown ) return;
  e = e || window.event;
  if ( mousex >= 0 && mousey >= 0 ) {
    var target = e.target ? e.target : e.srcElement;
    viewmat = mxrot3d(viewmat, 180.0 * (e.clientY - mousey) / target.height);
    viewmat = myrot3d(viewmat, 180.0 * (e.clientX - mousex) / target.width);
    paint(); // defined in ljmain.js
  }
  mousex = e.clientX;
  mousey = e.clientY;
  //console.log("mousemove", e.clientX, e.clientY, m2str(viewmat));
}



/* for mouse wheel event */
function ljwheel(e){
  var delta = 0; // positive for scrolling up
  e = e || window.event;
  if ( e.wheelDelta ) { // IE/Opera
    delta = e.wheelDelta / 120;
  } else if ( e.detail ) { // Firefox
    delta = -e.detail / 3;
  }
  if ( delta > 0 ) {
    userscale *= 1.05;
  } else if ( delta < 0 ) {
    userscale *= 0.95;
  }
  //console.log("wheel", delta);
  if ( e.preventDefault ) {
    e.preventDefault();
  }
  e.returnValue = false;
}



function transform(x, l)
{
  var n = x.length;
  var xyz = newarr(n), xc = [l * 0.5, l * 0.5, l * 0.5], xi = [0, 0, 0];

  for ( var i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    xyz[i] = mmulv(viewmat, xi);
    //console.log(x[i], xi, xc, xyz[i]);
    vinc(xyz[i], xc);
  }
  return xyz;
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = [0, 0, 0];
  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }
  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    for ( j = i + 1; j < n; j++ ) {
      if ( xyz[j][2] < xyz[k][2] ) {
        l = k;
        k = j;
        j = l;
      }
    }
    if ( k != i ) {
      vcopy(rt, xyz[k]);
      vcopy(xyz[k], xyz[i]);
      vcopy(xyz[i], rt);
    }
  }
  return xyz;
}



// draw all atoms in the box
function ljdraw3d(lj, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (lj.l + 1.0);
  var radius = 0.5 * scale;

  var xyz = transform(lj.x, lj.l); // apply the rotation matrix
  xyz = sortbyz(xyz); // sort particles by the z order

  // draw each particle
  var zmax = xyz[lj.n - 1][2], zmin = xyz[0][2];
  for (var i = 0; i < lj.n; i++) {
    var x = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scale + width  * 0.5 );
    var y = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scale + height * 0.5 );
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    var spotcolor = rgb2str(100 + 100 * zf, 100 + 100 * zf, 120 + 100 * zf);
    var color = rgb2str(20, 32, 80 + 160 * zf);
    // make closer particles larger
    var rz = Math.floor( radius * (0.7 + 0.3 * zf) );
    drawBall(ctx, x, y, rz, color, spotcolor);
  }
}



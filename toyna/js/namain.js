/* Handle web interface */



"use strict";



var na = null;
var nr = 20;
var tp = 1.5;
var rc = 1000.0;

var timer_interval = 100; // in milliseconds
var natimer = null;
var simulmethod = "MD";

var mddt = 0.002;
var thdt = 0.02;
var nstepspsmd = 200; // number of steps per second for MD
var nstepspfmd = 20;  // number of steps per frame for MD

var mcamp = 0.2;
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var mctot = 0.0;
var mcacc = 0.0;

var sum1 = 1e-30;
var sumU = 0.0;

var userscale = 1.0;



function getparams()
{
  nr = get_int("nr", 55);
  tp = get_float("temperature", 1.5);
  rc = get_float("rcutoff", 1000.0);

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepspsmd = get_int("nstepspersecmd", 1000);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  userscale = get_float("nascale");
}



function changescale()
{
  userscale = get_float("nascale");
  paint();
}



var mousedown = false;
var mousex = -1;
var mousey = -1;
var viewmat = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]];



function namousedown(e)
{
  e = e || window.event;
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown = true;
  //console.log("mousedown", e.clientX, e.clientY, m2str(viewmat));
}



function namouseup(e)
{
  e = e || window.event;
  mousex = -1;
  mousey = -1;
  mousedown = false;
  //console.log("mouseup", e.clientX, e.clientY, m2str(viewmat));
}



function namousemove(e)
{
  if ( !mousedown ) return;
  e = e || window.event;
  if ( mousex >= 0 && mousey >= 0 ) {
    var target = e.target ? e.target : e.srcElement;
    viewmat = mxrot3d(viewmat, 180.0 * (e.clientY - mousey) / target.height);
    viewmat = myrot3d(viewmat, 180.0 * (e.clientX - mousex) / target.width);
    paint();
  }
  mousex = e.clientX;
  mousey = e.clientY;
  //console.log("mousemove", e.clientX, e.clientY, m2str(viewmat));
}



/* for the mouse wheel event */
function nawheel(e)
{
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
  grab("nascale").value = userscale;
  //console.log("wheel", delta);
  if ( e.preventDefault ) {
    e.preventDefault();
  }
  e.returnValue = false;
  paint();
}



/* install the mouse wheel event */
function installwheel(target, handler)
{
  if ( target.addEventListener ) {
    // for IE9+, Chrome, Safari, Opera
    target.addEventListener('mousewheel', handler, false);
    // for Firefox
    target.addEventListener('DOMMouseScroll', handler, false);
  } else { // for IE 6/7/8
    target.attachEvent("onmousewheel", handler);
  }
}



function installmouse()
{
  var target = grab("nabox");
  target.onmousedown = namousedown;
  target.onmouseup = namouseup;
  target.onmousemove = namousemove;
  installwheel(target, nawheel);
}



function domd()
{
  var istep, sinfo = "";

  //for ( istep = 0; istep < nstepspfmd; istep++ ) {
  //  na.vv(mddt);
  //  na.vrescale(tp, thdt);
  //  sum1 += 1.0;
  //  sumU += na.epot;
  //}
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  return sinfo;
}



function domc()
{
  var istep, sinfo = "";

  //for ( istep = 0; istep < nstepspfmc; istep++ ) {
  //  mctot += 1.0;
  //  mcacc += na.metro(mcamp, 1.0 / tp);
  //  sum1 += 1.0;
  //  sumU += na.epot;
  //}
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%, ";
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  return sinfo;
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
function nadraw(na, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (na.l + 1.0);

  var xyz = transform(na.x, na.l); // apply the rotation matrix
  xyz = sortbyz(xyz); // sort particles by the z order

  // draw each particle
  var zmax = xyz[na.n - 1][2], zmin = xyz[0][2];
  for (var i = 0; i < na.n; i++) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var scl = scale * (0.7 + 0.3 * zf);
    var x = Math.floor(  (xyz[i][0] - na.l * 0.5) * scl + width  * 0.5 );
    var y = Math.floor( -(xyz[i][1] - na.l * 0.5) * scl + height * 0.5 );
    var spotcolor = rgb2str(100 + 100 * zf, 100 + 100 * zf, 120 + 100 * zf);
    var color = rgb2str(20, 32, 80 + 160 * zf);
    var rz = Math.floor( 0.5 * scl );
    drawBall(ctx, x, y, rz, color, spotcolor);
  }
}


function paint()
{
  if ( !na ) return;
  nadraw(na, "nabox", userscale);
}



function pulse()
{
  var sinfo;

  if ( simulmethod === "MD" ) {
    sinfo = domd();
  } else if ( simulmethod === "MC" ) {
    sinfo = domc();
  }
  grab("sinfo").innerHTML = sinfo;

  paint();
}



function stopsimul()
{
  if ( natimer !== null ) {
    clearInterval(natimer);
    natimer = null;
  }
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
  munit(viewmat);
}



function pausesimul()
{
  if ( !na ) return;
  if ( natimer !== null ) {
    clearInterval(natimer);
    natimer = null;
    grab("pause").value = "Resume";
  } else {
    natimer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").value = "Pause";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  na = new NA(nr, rc);
  na.force();
  installmouse();
  natimer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( natimer !== null ) {
    startsimul();
  }
}


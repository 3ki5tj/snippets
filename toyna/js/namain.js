/* Handle web interface */



"use strict";



var na = null;
var nr = 20;
var tp = 300.0;
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
  tp = get_float("temperature", 300.0);
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
  if ( !mousedown ) {
    return;
  }
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

  for ( istep = 0; istep < nstepspfmd; istep++ ) {
    na.vv(mddt);
    na.vrescale(tp, thdt);
    sum1 += 1.0;
    sumU += na.epot;
  }
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



function transform(x)
{
  var i, d, n = x.length, l = 0;
  var xyz = newarr(n), xc = [0, 0, 0], xi = [0, 0, 0];

  // compute the center of mass
  for ( i = 0; i < n; i++ ) {
    vinc(xc, x[i]);
  }
  vsmul(xc, 1.0/n);

  // rotate the coordinates of each particle
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    xyz[i] = mmulv(viewmat, xi);
    //console.log(x[i], xi, xc, xyz[i]);
    //vinc(xyz[i], xc);
    for ( d = 0; d < D; d++ ) {
      l = Math.max( Math.abs( xi[d] ), l );
    }
  }
  return [xyz, xc, l];
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = newarr(D);
  var idmap = newarr(n);
  var invmap = newarr(n);

  for ( i = 0; i < n; i++ ) {
    idmap[i] = i;
    // i:         index of the output array `xyz`
    // idmap[i]:  index of the input array `x`
    // so xyz[i] --> x[ idmap[i] ];
    invmap[i] = i;
  }

  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }

  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    var zmin = x[ idmap[i] ][2];
    for ( j = i + 1; j < n; j++ ) {
      if ( x[ idmap[j] ][2] < zmin ) {
        k = j;
        zmin = x[ idmap[j] ][2];
      }
    }
    if ( k != i ) {
      // before
      //  xyz[i] --> x[ idmap[i] ]
      //  xyz[k] --> x[ idmap[k] ]
      // after
      //  xyz[i] --> x[ idmap[k] ]
      //  xyz[k] --> x[ idmap[i] ]
      l = idmap[i];
      idmap[i] = idmap[k];
      idmap[k] = l;
    }
  }

  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[ idmap[i] ]);
  }
  // compute the inverse map
  for ( i = 0; i < n; i++ ) {
    invmap[ idmap[i] ] = i;
  }
  return [xyz, idmap, invmap];
}



// draw a metallic line
function drawLineFancy(ctx, xi, yi, xj, yj)
{
  drawLine(ctx, xi, yi, xj, yj, '#aaaaaa', 4);
  drawLine(ctx, xi, yi, xj, yj, '#bbbbbb', 2);
  drawLine(ctx, xi, yi, xj, yj, '#cccccc', 1);
}



// get the scaling factor due to z
function getzscale(r, zmin, zmax, ortho)
{
  if ( ortho ) {
    return 0.7;
  } else {
    var zf = (r[2] - zmin) / (zmax - zmin);
    return 0.7 + 0.3 * zf;
  }
}



// draw all atoms in the box
function nadraw(na, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var i, j, jb, k, ir, ic, ret;
  var apr = na.apr;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ret = transform(na.x); // apply the rotation matrix
  var xt = ret[0];
  na.l = ret[2];
  ret = sortbyz(xt); // sort particles by the z order
  var xyz = ret[0];
  var idmap = ret[1], invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]

  var ortho = grab("orthographic").checked;
  var scale = userscale * Math.min(width, height) / (na.l * 2.0);

  // draw each particle
  var zmax = xyz[na.n - 1][2], zmin = xyz[0][2];

  // draw lines that were used to group clusters
  {
    ctx.lineWidth = 2;
    ctx.strokeStyle = '#808080';
    for ( ir = 0; ir < na.nr; ir++ ) {
      i = invmap[ ir*apr ];
      j = invmap[ ir*apr + 1 ];
      if ( apr == 3 ) {
        jb = invmap[ ir*apr + 2 ];
      }

      var scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
      var xi = Math.floor(  xyz[i][0] * scli + width  * 0.5 );
      var yi = Math.floor( -xyz[i][1] * scli + height * 0.5 );

      var sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
      var xj = Math.floor(  xyz[j][0] * sclj + width  * 0.5 );
      var yj = Math.floor( -xyz[j][1] * sclj + height * 0.5 );

      var scljb = scale * getzscale(xyz[jb], zmin, zmax, ortho);
      var xjb = Math.floor(  xyz[jb][0] * scljb + width  * 0.5 );
      var yjb = Math.floor( -xyz[jb][1] * scljb + height * 0.5 );

      if ( apr === 2 ) {
        drawLineFancy(ctx, xi, yi, xj, yj);
      } else {
        drawLineFancy(ctx, xi, yi, xj, yj);
        drawLineFancy(ctx, xj, yj, xjb, yjb);
      }

      if ( ir < na.nr - 1 ) {
        k = invmap[ (ir + 1) * apr ];
        var sclk = scale * getzscale(xyz[k], zmin, zmax, ortho);
        var xk = Math.floor(  xyz[k][0] * sclk + width  * 0.5 );
        var yk = Math.floor( -xyz[k][1] * sclk + height * 0.5 );
        if ( apr === 2 ) {
          drawLineFancy(ctx, xi, yi, xk, yk);
        } else {
          drawLineFancy(ctx, xj, yj, xk, yk);
        }
      }
    }
  }


  for (i = 0; i < na.n; i++) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var scl = scale * getzscale(xyz[i], zmin, zmax, ortho);
    var x = Math.floor(  xyz[i][0] * scl + width  * 0.5 );
    var y = Math.floor( -xyz[i][1] * scl + height * 0.5 );
    var spotcolor = rgb2str(100 + 100 * zf, 100 + 100 * zf, 120 + 100 * zf);
    var color, rad;
    var i0 = idmap[ i ];
    var tp = i0 % apr;
    if ( tp === 0 ) {
      color = rgb2str(20, 32, 80 + 160 * zf);
      rad = 1.5;
    } else if ( tp === apr - 1 ) {
      color = rgb2str(80 + 160 * zf, 32, 20);
      rad = 1.0;
    } else {
      color = rgb2str(20, 80 + 60 * zf, 20);
      rad = 1.2;
    }
    var rz = Math.floor( rad * scl );
    paintBall(ctx, x, y, rz, color, spotcolor);
  }
}



function paint()
{
  if ( na ) {
    nadraw(na, "nabox", userscale);
  }
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
  if ( !na ) {
    return;
  }
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
  na = new NA("ACGGUUCAGCU", rc);
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


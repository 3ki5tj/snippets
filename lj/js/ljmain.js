/* Handle web interface */



"use strict";



var lj = null;
var n = 55;
var rho = 0.7;
var tp = 1.5;
var rcdef = 1000.0;

var timer_interval = 100; // in milliseconds
var ljtimer = null;
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
var sumP = 0.0;

var userscale = 1.0;



function getparams()
{
  n = get_int("n", 55);
  var dim = get_int("dimension", 2);
  if ( dim === 2 || dim === 3 ) {
    D = dim;
  }
  rho = get_float("density", 0.7);
  tp = get_float("temperature", 1.5);
  rcdef = get_float("rcutoff", 1000.0);

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepspsmd = get_int("nstepspersecmd", 1000);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  userscale = get_float("ljscale");
}



function changescale()
{
  userscale = get_float("ljscale");
  paint();
}



/* for the mouse wheel event */
function ljwheel(e)
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
  grab("ljscale").value = userscale;
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
  var target = grab("ljbox");
  target.onmousedown = ljmousedown;
  target.onmouseup = ljmouseup;
  target.onmousemove = ljmousemove;
  installwheel(target, ljwheel);
}



function domd()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmd; istep++ ) {
    lj.vv(mddt);
    lj.vrescale(tp, thdt);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
  }
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  sinfo += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3) + ".";
  return sinfo;
}



function domc()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    mctot += 1.0;
    mcacc += lj.metro(mcamp, 1.0 / tp);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
  }
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%, ";
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  sinfo += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3) + ".";
  return sinfo;
}



function paint()
{
  if ( !lj ) return;
  if ( lj.dim === 2 ) {
    ljdraw2d(lj, "ljbox", userscale);
  } else if ( lj.dim === 3 ) {
    ljdraw3d(lj, "ljbox", userscale);
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
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
  sumP = 0.0;
  munit(viewmat);
}



function pausesimul()
{
  if ( !lj ) return;
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
    grab("pause").value = "Resume";
  } else {
    ljtimer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").value = "Pause";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  lj = new LJ(n, D, rho, rcdef);
  lj.force();
  installmouse();
  ljtimer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( ljtimer !== null ) {
    startsimul();
  }
}


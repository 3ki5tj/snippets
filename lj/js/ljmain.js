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
var thtype = "vrescale"; // thermostat type
var thdamp = 5;  // thermostat damping factor
var nstepspsmd = 200; // number of steps per second for MD
var nstepspfmd = 20;  // number of steps per frame for MD

var mcamp = 0.2;
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var mctot = 0.0;
var mcacc = 0.0;

// reference values from the equation of state
var Uref;
var Pref;

var sum1 = 1e-30;
var sumU = 0.0;
var sumP = 0.0;
var sumK = 0.0;




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

  if ( dim === 3 ) {
    var ref = lj_eos3dPVEhBH(rho, tp);
    Uref = ref[0];
    Pref = ref[1];
  } else {
    Uref = null;
    Pref = null;
  }

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thtype = grab("thermostattype").value;
  thdamp = get_float("thdamp", 5);
  nstepspsmd = get_int("nstepspersecmd", 1000);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  mousescale = get_float("ljscale");
}



function changescale()
{
  mousescale = get_float("ljscale");
  paint();
}



function thermostat(dt)
{
  if ( thtype === "vrescale" ) {
    return lj.vrescale(tp, dt * thdamp * 2);
  } else {
    return lj.vlang(tp, dt * thdamp);
  }
}



function reportUP()
{
  var s;
  s = '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3);
  if ( Uref ) s += " (" + roundto(Uref, 3) + ")";
  s += ", ";
  s += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3);
  if ( Pref ) s += " (" + roundto(Pref, 3) + ")";
  s += ", " + sum1 + " samples.";
  return s;
}



function domd()
{
  var dof = ( thtype === "Langevin" ) ? lj.dim * lj.n : lj.dof;

  for ( var istep = 0; istep < nstepspfmd; istep++ ) {
    thermostat(mddt * 0.5);
    lj.vv(mddt);
    var ek = thermostat(mddt * 0.5);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
    sumK += ek / dof;
  }
  var sinfo = '<span class="math"><i>E<sub>K</sub></i>/<i>N<sub>f</sub></i></span>: '
        + roundto(sumK/sum1, 3) + "(" + roundto(tp/2, 3) + "), ";
  return sinfo + reportUP();
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
  return sinfo + reportUP();
}



function paint()
{
  if ( !lj ) {
    return;
  }
  if ( lj.dim === 2 ) {
    ljdraw2d(lj, "ljbox", mousescale);
  } else if ( lj.dim === 3 ) {
    ljdraw3d(lj, "ljbox", mousescale);
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



function resetdata()
{
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
  sumP = 0.0;
  sumK = 0.0;
}



function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  resetdata();
  munit(viewmat);
}



function pausesimul()
{
  if ( !lj ) {
    return;
  }
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
  installmouse("ljbox", "ljscale");
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


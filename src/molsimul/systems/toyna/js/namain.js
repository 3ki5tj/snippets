/* Handle web interface */



"use strict";



var na = null;
var seq = ""
var nr = 20;
var tp = 300.0;
var conc = 1.0; // salt concentration
var debyel = 4.36; // Debye screening length
var uhb0 = 2.43; // magnitude of the hydrogen bond

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

var mousescale = 1.0;



function getparams()
{
  seq = grab("sequence").value.trim();
  tp = get_float("temperature", 300.0);

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepspsmd = get_int("nstepspersecmd", 1000);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  mousescale = get_float("nascale");

  var conc = get_float("saltconc");
  debyel = getDebyel([1.0], [conc], 1, tp);
  uhb0 = get_float("uhb0");
}



function changescale()
{
  mousescale = get_float("nascale");
  paint();
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
  sinfo += '<span class="math"><i>U</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
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



function paint()
{
  if ( na ) {
    nadraw(na, "nabox", mousescale);
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
  na = new NA(seq, tp, debyel, uhb0);
  na.force();
  installmouse("nabox", "nascale");
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


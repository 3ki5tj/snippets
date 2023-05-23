/* Handle web interface */



"use strict";



var potts = null;
var l = 16;
var q = 10;
var tp = 0.7;

var timer_interval = 100; // in milliseconds
var potts_timer = null;
var mc_algorithm = "Metropolis";

var nstepspsmc = 1000; // number of steps per second for MC
var nstepspfmc = 100;  // number of steps per frame for MC
var mctot = 0.0;
var mcacc = 0.0;

var sum1 = 1e-300;
var sumU = 0;




function getparams()
{
  l = get_int("L", 32);
  q = get_int("q", 10);
  tp = get_float("temperature", 0.7);

  mc_algorithm = grab("mc_algorithm").value;

  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  mousescale = get_float("pottsscale");
}



function changescale()
{
  mousescale = get_float("pottsscale");
  paint();
}



function dometropolis()
{
  var istep, sinfo = "";
  var i, id, n = potts.n;

  potts.setproba(1.0/tp);
  //potts.energy();
  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    for ( i = 0; i < n; i++ ) {
      mctot += 1.0;
      id = potts.pick();
      if ( potts.h <= 0 || rand01() < potts.proba[potts.h] ) {
        mcacc += 1;
        potts.flip(id, potts.sn, potts.h);
      }
      sum1 += 1.0;
      sumU += potts.E;
    }
  }
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%, ";
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: '
         + roundto(sumU/sum1/n, 3) + ".";
  return sinfo;
}



function dowolff()
{
  var istep, sinfo = "";
  var i, id, n = potts.n;
  var padd;

  padd = 1 - Math.exp(-1/tp);
  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    potts.wolff(padd);
    sum1 += 1.0;
    sumU += potts.E;
  }
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: '
         + roundto(sumU/sum1/n, 3) + ".";
  return sinfo;
}



function paint()
{
  if ( !potts ) {
    return;
  }

  /* draw all atoms in the box */
  var c = grab("pottsbox");
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + 1
  var dx = 1.0 * width / (potts.l + 1) * mousescale;
  var dy = 1.0 * height / (potts.l + 1) * mousescale;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // draw each spin
  var l = potts.l, id = 0, i, j, s;
  var q = potts.q, colors = new Array(q), spotcolors = new Array(q);
  for ( s = 0; s < q; s++ ) {
    colors[s] = getHueColor(1.0*s/q);
    spotcolors[s] = lightenColor(colors[s], 0.3);
  }
  for ( i = 0; i < l; i++ ) {
    for ( j = 0; j < l; j++, id++ ) {
      var x = (i + 0.5 - 0.5*l) * dx + width * 0.5;
      var y = (j + 0.5 - 0.5*l) * dy + height * 0.5;
      var radius = dx * 0.5;
      s = potts.s[id];
      var color = colors[s];
      var spotcolor = spotcolors[s];
      paintBall(ctx, x, y, radius, color, spotcolor);
    }
  }
}



function pulse()
{
  var sinfo;

  if ( mc_algorithm === "Metropolis" ) {
    sinfo = dometropolis();
  } else if ( mc_algorithm === "Wolff" ) {
    sinfo = dowolff();
  }
  grab("sinfo").innerHTML = sinfo;

  paint();
}



function stopsimul()
{
  if ( potts_timer !== null ) {
    clearInterval(potts_timer);
    potts_timer = null;
  }
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
}



function pausesimul()
{
  if ( !potts ) {
    return;
  }
  if ( potts_timer !== null ) {
    clearInterval(potts_timer);
    potts_timer = null;
    grab("pause").value = "Resume";
  } else {
    potts_timer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").value = "Pause";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  potts = new Potts(l, q);
  installmouse("pottsbox", "pottsscale");
  potts_timer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( potts_timer !== null ) {
    startsimul();
  }
}


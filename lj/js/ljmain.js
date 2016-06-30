/* Handle web interface */



"use strict";



var lj = null;
var n = 55;
var rho = 0.7;
var tp = 1.5;
var rcdef = 1000.0;
var dof;

var timer_interval = 100; // in milliseconds
var ljtimer = null;
var simulmethod = "MD";

var mddt = 0.002;
var thtype = "v-rescale"; // thermostat type
var vresdamp = 20.0;  // velocity-rescaling damping factor
var langdamp = 1.0;  // Langevin-dynamics damping factor
var zeta = null;
var zmass = null;

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

var dke = 1.0;
var kehistn = 0;
var kehist = null;
var histplot = null;

var keseq = null;
var keseqmax = 50000;

var vsize = 5;
var vseq = null
var vseqmax = 50000;

var corrplot = null;

var paintcnt = 0;


function getparams()
{
  stopsimul();

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
  vresdamp = get_float("vresdamp", 20.0);
  langdamp = get_float("langdamp", 1.0);
  nstepspsmd = get_int("nstepspersecmd", 1000);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  var nhclen = get_int("nhclen", 5);
  var nhcmass1 = get_float("nhcmass1", 1.0);
  var nhcmass2 = get_float("nhcmass2", 1.0);
  var i;
  zeta = newarr(nhclen);
  zmass = newarr(nhclen);
  for ( i = 0; i < nhclen; i++ ) {
    zeta[i] = 0.0;
    zmass[i] = ( i === 0 ) ? nhcmass1 : nhcmass2;
  }

  vsize = get_int("vactsize", 1);

  lj = new LJ(n, D, rho, rcdef);
  lj.force();
  dof = ( thtype === "Langevin" ) ? lj.dim * lj.n : lj.dof;

  kehistn = Math.floor(tp * dof / dke) * 2;
  kehist = new Array(kehistn);
  for ( i = 0; i < kehistn; i++ ) kehist[i] = 0;

  keseq = [];
  vseq = [];

  mousescale = get_float("ljscale");

  paintcnt = 0;
}



function changescale()
{
  mousescale = get_float("ljscale");
  paint();
}



function thermostat(dt)
{
  if ( thtype === "v-rescale" ) {
    return lj.vrescale(tp, dt * vresdamp);
  } else if ( thtype === "Nose-Hoover" ) {
    return lj.nhchain(tp, dt, zeta, zmass);
  } else if (thtype === "Langevin" ) {
    return lj.langevin(tp, dt * langdamp);
  }
}



function reportUP()
{
  var s;
  s = '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3);
  if ( Uref ) s += " (ref: " + roundto(Uref, 3) + ")";
  s += ". ";
  s += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3);
  if ( Pref ) s += " (ref: " + roundto(Pref, 3) + ")";
  s += ".<br>" + sum1 + " samples.";
  return s;
}



function domd()
{
  for ( var istep = 0; istep < nstepspfmd; istep++ ) {
    thermostat(mddt * 0.5);
    lj.vv(mddt);
    var ekin = thermostat(mddt * 0.5);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
    sumK += ekin / dof;
    var ike = Math.floor(ekin / dke);
    if ( ike < kehistn ) kehist[ike] += 1.0;
    keseq.push( 2 * ekin / (dof * tp) - 1 );
    for ( var j = 0; j < vsize; j++ ) {
      vseq.push( lj.v[j][0] );
      vseq.push( lj.v[j][1] );
      vseq.push( lj.v[j][2] );
    }
  }
  if ( keseq.length > keseqmax )
    keseq = keseq.slice( keseq.length - keseqmax );
  if ( vseq.length > vseqmax * vsize * D )
    vseq = vseq.slice( vseq.length - vseqmax * vsize * D );
  var sinfo = '<span class="math"><i>E<sub>K</sub></i>/<i>N<sub>f</sub></i></span>: '
        + roundto(sumK/sum1, 3) + " (" + roundto(tp/2, 3) + ").<br>";
  return sinfo + reportUP();
}



function domc()
{
  for ( var istep = 0; istep < nstepspfmc; istep++ ) {
    mctot += 1.0;
    mcacc += lj.metro(mcamp, 1.0 / tp);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
  }
  var sinfo = "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%.<br>";
  return sinfo + reportUP();
}



// return ln Gamma(dof/2)
function getnorm(dof)
{
  var s = (dof % 2) ? 0.5 * Math.log( Math.PI ) : 0;
  for ( var i = 2 - dof % 2; i < dof; i += 2 )
    s += Math.log(i*0.5);
  return s;
}

/* update the histogram plot */
function updatehistplot()
{
  if ( simulmethod === "MC" ) {
    histplot = null;
    return;
  }

  var i, htot = 0, norm = getnorm(dof);
  var dat = "Kinetic energy,Histogram,Reference\n";
  for ( i = 0; i < kehistn; i++ )
    htot += kehist[i];
  var kemin = Math.max(0, (dof - 5 * Math.sqrt(dof)) * tp * 0.5);
  var kemax = (dof + 5 * Math.sqrt(dof)) * tp * 0.5;
  for ( i = 0; i < kehistn; i++ ) {
    var ke = (i + 0.5) * dke;
    if ( ke < kemin ) continue;
    if ( ke > kemax ) break;
    var hval = kehist[i] / (htot * dke);
    var ref = Math.exp(Math.log(ke/tp) * (dof*0.5-1) -ke/tp - norm) / tp;
    dat += "" + ke + "," + hval + "," + ref + "\n";
  }

  if ( histplot === null ) {
    var options = {
      xlabel: '<small>Kinetic energy</small>',
      ylabel: '<small>Histogram',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      width: 360,
      height: 240,
      xRangePad: 2,
    };
    histplot = new Dygraph(document.getElementById("histplot"), dat, options);
  } else {
    histplot.updateOptions({ file: dat });
  }
}



/* update the correlation plot */
function updatecorrplot()
{
  if ( simulmethod === "MC" ) {
    corrplot = null;
    return;
  }

  var dat = "Time,Kinetic energy,Velocity\n";
  var i, j, jmax, iblk, l, t;
  var kelen = keseq.length, kecorr0, kecorr;
  var vblk = vsize * D, vlen = vseq.length / vblk, vcorr0, vcorr;
  //console.log(kelen, vlen);
  for ( i = 0; i <= 1000; i++ ) {
    if ( i >= kelen && i >= vlen ) break;
    
    jmax = kelen - i;
    kecorr = 0;
    for ( j = 0; j < jmax; j++ )
      kecorr += keseq[j] * keseq[j+i];
    kecorr /= jmax;

    jmax = (vlen - i) * vblk;
    iblk = i * vblk;
    vcorr = 0;
    for ( j = 0; j < jmax; j++ )
      vcorr += vseq[j] * vseq[j + iblk];
    vcorr /= jmax;

    if ( i === 0 ) {
      kecorr0 = kecorr;
      vcorr0  = vcorr;
    }
    dat += "" + roundto(i*mddt, 4) + ","
        + roundto(kecorr/kecorr0, 4) + ","
        + roundto(vcorr/vcorr0, 4) + "\n";
  }

  if ( corrplot === null ) {
    var options = {
      xlabel: '<small>Time</small>',
      ylabel: '<small>Normalized autocorrelation function</small>',
      includeZero: true,
      axisLabelFontSize: 10,
      width: 360,
      height: 240,
      xRangePad: 2,
    };
    corrplot = new Dygraph(document.getElementById("corrplot"), dat, options);
  } else {
    corrplot.updateOptions({ file: dat });
  }
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
  ++paintcnt;
  if ( paintcnt % 10 == 0 ) updatehistplot();
  if ( paintcnt % 100 == 0 ) updatecorrplot();
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
  for ( var i = 0; i < kehistn; i++ ) kehist[i] = 0;
  keseq = [];
  vseq = []; 
  histplot = null;
  corrplot = null;
  paintcnt = 0;
}



function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  resetdata();
  munit(viewmat);
  grab("pause").value = "Pause";
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
  getparams();
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


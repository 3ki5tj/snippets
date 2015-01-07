


"use strict";



var N = 3;


var hess, hessmat;
var mass;
var x;
var v;
var f;
var temp = 1; // temperature
var dt = 0.005; // time step for molecular dynamics
var thdt = 0.01; // thermostat time step
var thermostat = "Langevin"; // thermostat type

var quartic_a = 0.0;
var quartic_b = 0.0;

var fsmallworld = 0.0;
var ksmallworld = 1.0;
var smallworld;

var nstepsperframe = 1000; // number of steps per frame
var timer_interval = 100; // interval of animation, in milliseconds

var nstequiv = 10000; // number of steps for equilibration
var nstpca = 1000000; // number of steps of doing PCA

var nstsamp = 10; // frequency of sampling
var nsamps; // number of samples collected so far
var qsum;
var qqsum; // correlations

// arrays used in pca()
var avq;
var varqq;
var invvarqq;
var eigval;
var eigvec;

var mdstep = 0;
var smcnt = 0, smT = 0, smU = 0;



/* allocate array spaces */
function getparams()
{
  N = get_int("npt", 3);
  hess = newarr(N);
  hessmat = newarr(N*N);
  mass = newarr(N);
  smallworld = newarr(N*N);
  x = newarr(N);
  v = newarr(N);
  f = newarr(N);
  qsum = newarr(N);
  qqsum = newarr(N*N);
  avq = newarr(N);
  varqq = newarr(N*N);
  invvarqq = newarr(N*N);
  eigval = newarr(N);
  eigvec = newarr(N*N);
  nsamps = 0;
  nstequiv = get_int("nstequiv", 10000);
  mddt = get_float("mddt", 0.005);
  temp = get_float("temp", 1.0);
  thdt = get_float("thermostatdt", 0.01);
  thermostat = grab("thermostat").value;
  nstsamp = get_int("nstsamp", 10);
  quartic_a = get_float("quartic_a", 0.0);
  quartic_b = get_float("quartic_b", 0.0);
}



/* initialize the force field */
function initff()
{
  var i, j;

  for ( i = 0; i < N; i++ ) {
    hess[i] = 1;
    if ( i > 0 ) {
      j = i - 1;
      hessmat[i*N + i] += 1;
      hessmat[j*N + j] += 1;
      hessmat[i*N + j] -= 1;
      hessmat[j*N + i] -= 1;
    }
  }

  // Warm up the random number generator
  for ( i = 0; i < 1000; i++ ) rand01();

  fsmallworld = get_float("fsmallworld", 0.0);
  ksmallworld = get_float("ksmallworld", 1.0);
  // add small-world springs
  for ( i = 0; i < N; i++ )
    for ( j = i + 2; j < N; j++ ) {
      smallworld[i*N+j] = ( rand01() < fsmallworld );
      if ( smallworld[i*N+j] ) {
        hessmat[i*N+j] -= ksmallworld;
        hessmat[j*N+i] -= ksmallworld;
        hessmat[i*N+i] += ksmallworld;
        hessmat[j*N+j] += ksmallworld;
        console.log("add spring between " + (i+1) + " and " + (j+1), "K", ksmallworld);
      }
    }

  var massdistr = grab("massdistr").value;
  var massslope = get_float("massslope", 1);
  for ( i = 0; i < N; i++ ) {
    if ( massdistr == "triangular" ) {
      mass[i] = -Math.abs(2.*i/(N-1) - 1)*massslope + 1;
      if ( massslope > 0 ) mass[i] += massslope;
    } else {
      mass[i] = 1;
    }
  }
  console.log("mass: ", mass);
}



function getpairpot(x, k, a, b)
{
  var dx2 = x*x - b*b;
  return [0.5 * k * x * x + 0.25 * a * dx2 * dx2,
          k * x + a * dx2 * x];
}



/* compute the force, return the energy */
function force()
{
  var i, j;
  var U = 0, tmp;

  for ( i = 0; i < N; i++ ) f[i] = 0;

  if ( N == 1 ) {
    tmp = getpairpot(-x[0], hess[0], quartic_a, quartic_b);
    U = tmp[0];
    f[0] = tmp[1];
  } else {
    for ( i = 0; i < N - 1; i++ ) { // loop over springs
      tmp = getpairpot(x[i+1]-x[i], hess[i], quartic_a, quartic_b);
      U       += tmp[0];
      f[i]    += tmp[1];
      f[i+1]  -= tmp[1];
    }

    // compute the force from small-world springs
    if ( fsmallworld > 0 ) {
      for ( i = 0; i < N; i++ ) {
        for ( j = i + 2; j < N; j++ ) {
          if ( smallworld[i*N + j] ) {
            tmp = getpairpot(x[j]-x[i], ksmallworld, quartic_a, quartic_b);
            U       += tmp[0];
            f[i]    += tmp[1];
            f[j]    -= tmp[1];
          }
        }
      }
    }
  }
  return U;
}



/* remove the center of mass motion */
function rmcom(arr)
{
  var i;
  var p = 0, mt = 0;

  if ( N <= 1 ) return;
  for ( i = 0; i < N; i++ ) {
    p += arr[i] * mass[i];
    mt += mass[i];
  }
  p /= mt;
  //console.log("com ", p);
  for ( i = 0; i < N; i++ )
    arr[i] -= p;
}



function getEk(vel)
{
  var ek = 0;

  for ( var i = 0; i < N; i++ )
    ek += .5 * mass[i] * vel[i] * vel[i];
  return ek;
}



/* Andersen thermostat
 * Note, DOF should be N in this case */
function andersen()
{
  var i;
  var ek;

  if ( rand01() < 0.1 ) {
    i = Math.floor(N * rand01());
    v[i] = Math.sqrt(temp/mass[i]) * randgaus();
  }
  return getEk(v);
}



/* Langevin-like thermostat */
function langevin(thdt)
{
  for ( var i = 0; i < N; i++ )
    v[i] += (-thdt * v[i] + Math.sqrt(2*temp*thdt) * randgaus())/mass[i];
  rmcom(v);
  return getEk(v);
}



/* velocity rescaling thermostat */
function vrescale(thdt, dof)
{
  var i, j;
  var ek1, ek2, s, c, r, r2;

  ek1 = getEk(v);
  c = (thdt < 700) ? Math.exp(-thdt) : 0;
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = c * ek1 + (1 - c) * (r2 + r*r) * .5 * temp
      + r * Math.sqrt(c * (1 - c) * 2 * temp * ek1);
  if (ek2 < 1e-30) ek2 = 1e-30;
  s = Math.sqrt(ek2/ek1);
  for (i = 0; i < N; i++)
    v[i] *= s;

  return ek2;
}



/* collect position correlations */
function sample()
{
  var i, j;

  for ( i = 0; i < N; i++ )
    qsum[i] += x[i];
  for ( i = 0; i < N; i++ )
    for ( j = i; j < N; j++ )
      qqsum[i*N+j] += x[i] * x[j];
  nsamps += 1;
}



function initmd()
{
  /* initialize the random velocities */
  for ( var i = 0; i < N; i++ ) {
    x[i] = 0.001 * randgaus();
    v[i] = Math.sqrt(temp/mass[i]) * randgaus();
  }
  rmcom(x);
  rmcom(v);
}



function domd()
{
  var Ep, Ek;
  var i, dof, istep;

  dof = ( thermostat == "Andersen" || N == 1 ) ? N : N - 1;
  for ( Ek = 0, i = 0; i < N; i++ )
    Ek += .5 * mass[i] * v[i] * v[i];
  Ep = force();
  //console.log("step", mdstep, Ek, Ep, dof, smT, smU, smcnt);
  for ( istep = 0; istep < nstepsperframe; istep++ ) {
    mdstep++;

    /* integrate Newton's equation of motion */
    for ( i = 0; i < N; i++ ) {
      v[i] += f[i]/mass[i] * .5 * dt;
      x[i] += v[i] * dt;
    }

    Ep = force();

    for ( i = 0; i < N; i++ ) {
      v[i] += f[i]/mass[i] * .5 * dt;
    }

    /* apply the thermostat */
    if ( thermostat == "Andersen" || N == 1 ) {
      Ek = andersen();
    } else if ( thermostat == "Langevin" ) {
      Ek = langevin(thdt);
    } else if ( thermostat == "vrescale" ) {
      Ek = vrescale(thdt, dof);
    } else {
      throw new Error("Unknown thermostat " + thermostat);
    }

    if ( mdstep <= nstequiv ) continue;

    smT += 2*Ek/dof;
    smU += Ep/N;
    smcnt += 1;

    if ( mdstep % nstsamp == 0 ) {
      if ( N > 1  && thermostat != "Andersen" ) rmcom(v);
      if ( N > 1 ) rmcom(x);
      sample();
    }
  }
  if ( N > 1 ) rmcom(x);
  if ( thermostat != "Andersen" && N > 1 ) rmcom(v);
  return "" + mdstep + " steps, " + nsamps + " samples, "
      + "Ek " + roundto(Ek, 2) + ", Ep " + roundto(Ep, 2)
      + ", Ek + Ep " + roundto(Ek + Ep, 2)
      + ", T " + roundto(smT/smcnt, 3)
      + ", U/N " + roundto(smU/smcnt, 3);
}



function printvec(v, n, name)
{
  var i, j;

  var s = name + ":\n";
  for ( i = 0; i < n; i++ )
    s += v[i] + " ";
  s += "\n";
  return s;
}

function prvec(m, nm)
{
  console.log( printvec(m, N, nm) );
}



function printmat(m, n, name)
{
  var i, j;

  var s = name + "\n";
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      s += m[i*n+j] + " ";
    s += "\n";
  }
  s += "\n";
  return s;
}

function prmat(m, nm)
{
  console.log( printmat(m, N, nm) );
}



var eigval_min = 1e-8;

/* add the zero modes */
function add_zero_modes(c, eigval, eigvec)
{
  var k, i, j;
  var big = N*100;
  var v = new Array(N);

  for ( k = 0; k < N; k++ ) {
    if ( eigval[k] > eigval_min ) continue;
    for ( i = 0; i < N; i++ )
      v[i] = eigvec[i*N + k];
    for ( i = 0; i < N; i++ )
      for ( j = 0; j < N; j++ )
        c[i*N + j] += v[i] * v[j] * big;
  }
}



/* principal component analysis */
function pca()
{
  var i, j;
  var tiny = 1e-6;

  if ( nsamps <= 0 ) return null;
  for ( i = 0; i < N; i++ )
    avq[i] = qsum[i]/nsamps;
  for ( i = 0; i < N; i++ )
    for ( j = i; j < N; j++ ) {
      var y = qqsum[i*N+j]/nsamps - avq[i] * avq[j];
      varqq[j*N+i] = varqq[i*N+j] =  y/temp*Math.sqrt(mass[i]*mass[j]);
      //console.log(i, j, y, avq[i], avq[j], qqsum[i*N+j], temp, mass[i], mass[j]);
    }
  //prmat(varqq, "q-q correlation");

  eigsym(varqq, eigval, eigvec, N);
  //prvec(eigval, "eigenvalues (omega^(-2))");
  //prmat(eigvec, "eigenvectors");

  // compute the frequencies
  var omgs = []
  var modes = []
  for ( i = 0; i < N; i++ ) {
    var omg = eigval[i];
    var mode = new Array(N);
    var mmax = 0;

    // copy the vector for mode i
    for ( j = 0; j < N; j++ ) {
      mode[j] = eigvec[j*N+i];
      if ( Math.abs(mode[j]) > Math.abs(mmax) )
        mmax = mode[j];
    }

    // normalize the mode such that the maximal value is 1
    if ( Math.abs(mmax) > 0 ) {
      for ( j = 0; j < N; j++ )
        mode[j] /= mmax;
    }

    if ( omg <= tiny ) { // a zero-frequency mode
      omg = 0;
      omgs.unshift( omg );
      modes.unshift( mode ); // add the mode at the beginning
    } else {
      omg = 1./Math.sqrt(omg);
      omgs.push( omg )
      modes.push( mode ); // add the mode at the end
    }
  }

  // compute the inverse of the correlation matrix
  // the zero modes must be added to avoid a singular matrix
  add_zero_modes(varqq, eigval, eigvec);
  //console.log(varqq, eigval, eigvec);
  luinv(varqq, invvarqq, N, 1e-10);
  for ( i = 0; i < N; i++ )
    for ( j = 0; j < N; j++ )
      invvarqq[i*N + j] *= Math.sqrt( mass[i] * mass[j] );

  return [omgs, modes, invvarqq];
}



var mdtimer = null;

var posplot = null;
var posymax = 0;
var omgplot = null;
var modesplot = null;



/* update the position plot */
function updateposplot()
{
  var i;
  var dat = "<i>i</i>,<i>x<sub>i</sub></i>\n";
  for ( i = 0; i < N; i++ ) {
    posymax = Math.max(Math.abs(x[i]), posymax);
    dat += "" + (i+1) + "," + x[i] + "\n";
  }
  if ( posplot == null ) {
    var options = {
      title: 'Positions of oscillators',
      xlabel: 'Particle index, <i>i</i>',
      ylabel: 'Position, <i>x</i><sub><i>i</i></sub>',
      includeZero: true,
      drawPoints: true,
      pointSize: 3,
      xRangePad: 3,
      //stepPlot: true,
      //yRangePad: 1,
      width: 750,
    };
    posplot = new Dygraph(grab("posplot"), dat, options);
  } else {
    posplot.updateOptions({
      file: dat,
      valueRange: [-posymax, posymax],
    });
  }
}



/* update the frequency plot */
function updateomgplot(omgs)
{
  var i;
  var dat = "<i>k</i>,<i>&omega;<sub>k</sub></i>\n";
  for ( i = 0; i < N; i++ )
    dat += "" + (i+1) + "," + omgs[i] + "\n";
  if ( omgplot == null ) {
    var options = {
      title: 'Frequencies of modes',
      xlabel: 'Mode, <i>k</i>',
      ylabel: 'Angular frequency, <i>&omega;</i><sub><i>k</i></sub>',
      includeZero: true,
      drawPoints: true,
      pointSize: 2,
      xRangePad: 2,
      width: 750,
    };
    omgplot = new Dygraph(grab("omgplot"), dat, options);
  } else {
    omgplot.updateOptions({ file: dat, });
  }
}



/* update the modes plot */
function updatemodesplot(omgs, modes)
{
  var i, j;

  var nmodes = get_int("nmodesplot", 5);
  nmodes = Math.min(N, nmodes);

  var dat = "<i>x</i>,";
  for ( j = 0; j < nmodes; j++ )
    dat += "Mode " + (j+1) + "(<i>&omega;</i> = " + roundto(omgs[j], 4) + "),";
  dat += "Zero\n";
  for ( i = 0; i < N; i++ ) {
    dat += "" + (i+1) + ",";
    for ( j = 0; j < nmodes; j++ )
      dat += "" + modes[j][i] + ",";
    dat += "0\n";
  }
  //console.log(dat);

  if ( modesplot == null ) {
    var options = {
      title: 'Modes',
      xlabel: 'Particle index, <i>i</i>',
      ylabel: '<i>x</i><sub><i>i</i></sub>',
      includeZero: true,
      drawPoints: true,
      pointSize: 1,
      xRangePad: 1,
      width: 750,
    };
    modesplot = new Dygraph(grab("modesplot"), dat, options);
  } else {
    modesplot.updateOptions({ file: dat, });
  }
}



/* draw the hessian matrix */
function drawhess(hmat, target)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var w = c.width;
  var h = c.height;
  var b = Math.min(w, h) / N;

  // find the maximal
  var mx = 0;
  for ( var i = 0; i < N; i++ )
    for ( var j = 0; j < N; j++ )
      mx = Math.max(Math.abs(hmat[i*N+j]), mx);

  ctx.clearRect(0, 0, w, h);
  for ( var i = 0; i < N; i++ ) {
    for ( var j = 0; j < N; j++ ) {
      var x = Math.floor( 255 * hmat[i*N+j] / mx );
      var red = 0, blue = 0, green = 0;
      if ( x >= 0 ) {
        blue = 255;
        red = 255 - x;
        green = 255 - x;
      } else {
        x = -x;
        red = 255;
        green = 255 - x;
        blue = 255 - x;
      }
      ctx.fillStyle = "rgb(" + red + ", " + green + ", " + blue + ")";
      ctx.fillRect(i*b, j*b, b, b);
    }
  }
}



/* upon a timer event */
function pulse()
{
  var sinfo = domd();
  grab("info").innerHTML = sinfo;
  updateposplot();
  nstpca = get_int("nstpca", 500000);
  if ( mdstep % nstpca == 0 ) {
    var ret = pca();
    if ( ret == null ) return;
    var omgs = ret[0];
    var modes = ret[1];
    var hessmatb = ret[2];
    updateomgplot(omgs);
    updatemodesplot(omgs, modes);
    drawhess(hessmat, "hessplot");
    drawhess(hessmatb, "hessbplot");
  }
}



function stopmd()
{
  if ( mdtimer != null ) {
    clearInterval(mdtimer);
    mdtimer = null;
    posymax = 0;
    posplot = null;
    omgplot = null;
    modesplot = null;
  }
}



function startmd()
{
  // stop the previous timer, if any
  stopmd();

  getparams();
  initff();
  initmd();

  nstepsperframe = get_int("nstepspersec", 1000)*timer_interval/1000;
  mdstep = 0;
  smcnt = 1e-30;
  smT = smU = 0;
  mdtimer = setInterval(function(){ pulse() }, timer_interval);
}



function changeparams()
{
  if ( mdtimer != null ) startmd();
}


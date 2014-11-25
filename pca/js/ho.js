var N = 3;


var hess;
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
var eval;
var evec;

var mdstep = 0;
var smcnt = 0, smT = 0, smU = 0;



/* allocate array spaces */
function getparams()
{
  N = get_int("npt", 3);
  hess = newnumarr(N);
  mass = newnumarr(N);
  smallworld = newnumarr(N*N);
  x = newnumarr(N);
  v = newnumarr(N);
  f = newnumarr(N);
  qsum = newnumarr(N);
  qqsum = newnumarr(N*N);
  avq = newnumarr(N);
  varqq = newnumarr(N*N);
  invvarqq = newnumarr(N*N);
  eval = newnumarr(N);
  evec = newnumarr(N*N);
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
  var i;

  for ( i = 0; i < N; i++ ) {
    hess[i] = 1;
  }

  /* Warm up the generator */ 
  for ( i = 0; i < 1000; i++ ) rand01();

  fsmallworld = get_float("fsmallworld", 0.0);
  ksmallworld = get_float("ksmallworld", 1.0);
  /* add small-world springs */
  for ( i = 0; i < N; i++ )
    for ( j = i + 2; j < N; j++ ) {
      smallworld[i*N+j] = ( rand01() < fsmallworld );
      if ( smallworld[i*N+j] )
        console.log("add spring between " + (i+1) + " and " + (j+1), fsmallworld);
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
            f[i+1]  -= tmp[1];
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



/* Langevin-like thermostat */
function langevin(thdt)
{
  for ( var i = 0; i < N; i++ )
    v[i] += (-thdt * v[i] + Math.sqrt(2*temp*thdt) * gaussrand())/mass[i];
  rmcom(v);
  return getEk(v);
}



/* randomly collide two particles */
function random_collision_MC()
{
  var m, frac, vi, vj, dv, r, amp;

  var i = Math.floor(N * rand01());
  var j = (Math.floor((N - 1) * rand01()) + i + 1) % N;

  frac = mass[i] / (mass[i] + mass[j]);
  m = mass[j] * frac;

  amp = Math.sqrt(2 * temp / m);
  dv = (rand01()*2 - 1) * amp;
  /* distribute dv to i and j such that
   * mass[i] * v[i] + mass[j] * v[j] is conserved */
  vi = v[i] + dv * (1 - frac);
  vj = v[j] - dv * frac;
  r = .5 * mass[i] * (vi*vi - v[i]*v[i])
    + .5 * mass[j] * (vj*vj - v[j]*v[j]);
  if ( r < 0 || rand01() < Math.exp(-r/temp) ) {
    v[i] = vi;
    v[j] = vj;
  }
  return getEk(v);
}




/* randomly collide two particles */
function random_collision_langevin(thdt)
{
  var m, frac, vij, dv;

  var i = Math.floor(N * rand01());
  var j = (Math.floor((N - 1) * rand01()) + i + 1) % N;

  frac = mass[i] / (mass[i] + mass[j]);
  m = mass[j] * frac;
  vij = v[i] - v[j];

  /* do a step of Langevin equation */
  dv = (-thdt * vij + Math.sqrt(2*temp*thdt) * gaussrand()) / m;
  /* distribute dv to i and j such that
   * mass[i] * v[i] + mass[j] * v[j] is conserved */
  v[i] += dv * (1 - frac);
  v[j] -= dv * frac;
  //console.log(i, j, dv, frac, m, temp, thdt, vij);
  return getEk(v);
}



/* velocity rescaling thermostat
 * do not use, not ergodic! */
function vrescale(thdt, dof)
{
  var i, j;
  var ekav = .5*temp*dof, ek1, ek2, s, amp;

  // hack: randomly swap two velocities to increase the randomness
  // for a complex system, we shouldn't need this
  // it is not ergodic, even with this hack
  if ( rand01() < 0.1 ) {
    i = Math.floor(N * rand01());
    j = (Math.floor((N - 1) * rand01()) + i + 1) % N;
    var tmp = mass[i] * v[i];
    v[i] = v[j] * mass[j] / mass[i];
    v[j] = tmp / mass[j];
  }

  // normal velocity rescaling
  for ( ek1 = 0, i = 0; i < N; i++ )
    ek1 += .5 * mass[i] * v[i] * v[i];
  amp = 2 * Math.sqrt(ek1*ekav*thdt/dof);
  ek2 = ek1 + (ekav - ek1)*thdt + amp*gaussrand();
  if (ek2 < 1e-6) ek2 = 1e-6;
  s = Math.sqrt(ek2/ek1);
  for (i = 0; i < N; i++)
    v[i] *= s;

  return ek2;
}



/* Andersen thermostat
 * Note, DOF should be N in this case */
function andersen()
{
  var i;
  var ek;

  if ( rand01() < 0.1 ) {
    i = Math.floor(N * rand01());
    v[i] = Math.sqrt(temp/mass[i]) * gaussrand();
  }
  for (ek = 0, i = 0; i < N; i++)
    ek += .5 * mass[i] * v[i] * v[i];
  return ek;
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
  for ( i = 0; i < N; i++ ) {
    x[i] = 0.001 * gaussrand();
    v[i] = Math.sqrt(temp/mass[i]) * gaussrand();
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
    } else if ( thermostat == "random_collision_MC" ) {
      Ek = random_collision_MC();
    } else if ( thermostat == "random_collision_Langevin" ) {
      Ek = random_collision_langevin(thdt);
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



/* principle component analysis */
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

  eigsym(varqq, eval, evec, N);
  //prvec(eval, "eigenvalues (omega^(-2))");
  //prmat(evec, "eigenvectors");

  // compute the frequencies
  omgs = []
  modes = []
  for ( i = 0; i < N; i++ ) {
    var omg = eval[i];
    var mode = new Array(N);
    var mmax = 0;
    for ( j = 0; j < N; j++ ) {
      mode[j] = evec[j*N+i];
      if ( Math.abs(mode[j]) > Math.abs(mmax) )
        mmax = mode[j];
    }
    // normalize the mode such that the maximal value is 1
    if ( Math.abs(mmax) > 0 ) {
      for ( j = 0; j < N; j++ )
        mode[j] /= mmax;
    }
    if ( omg <= tiny ) {
      omg = 0;
      omgs.unshift( omg );
      modes.unshift( mode ); // add the mode at the beginning
    } else {
      omg = 1./Math.sqrt(omg);
      omgs.push( omg )
      modes.push( mode ); // add the mode at the end
    }
  }

  return [omgs, modes]
}



var mdtimer = null;

var posplot = null;
var posymax = 0;
var omgplot = null;
var modesplot = null;



/* update the position plot */
function updateposplot()
{
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

  dat = "<i>x</i>,";
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



/* upon a timer event */
function pulse()
{
  sinfo = domd();
  grab("info").innerHTML = sinfo;
  updateposplot();
  nstpca = get_int("nstpca", 500000);
  if ( mdstep % nstpca == 0 ) {
    var pair = pca();
    if ( pair == null ) return;
    omgs = pair[0];
    modes = pair[1];
    updateomgplot(omgs);
    updatemodesplot(omgs, modes);
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
  getparams();
  initff();
  initmd();
  // stop the previous timer, if any
  stopmd();

  timer_interval = 100; // in milliseconds
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


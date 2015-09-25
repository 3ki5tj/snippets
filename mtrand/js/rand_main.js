


"use strict";



var timer = null;

var ngaus = 0;
var mu = []; // means
var sig = []; // width
var w = []; // relative fraction
var wacc = []; // cumulative function of `w`

// histogram
var hist = {
  xmin: 0.0,
  xmax: 20.0,
  dx:   0.1,
  n:    200,
  arr:  [],
  ref:  []
};

var histplot = null;



/* load parameters */
function loadparams(npt)
{
  var ngaus = get_int("ngaus", 1);
  var i, j, x, y;

  mu = newarr(ngaus);
  sig = newarr(ngaus);
  w = newarr(ngaus);

  // load Gaussian parameters from the web page
  for ( i = 0; i < ngaus; i++ ) {
    var i1 = i + 1;
    mu[i] = get_float("mu_" + i1, 0.0);
    sig[i] = get_float("sig_" + i1, 1.0);
    w[i] = get_float("w_" + i1, 1.0);
  }

  // compute the cumulative function of `w`
  wacc = newarr(ngaus + 1);
  wacc[0] = 0;
  for ( i = 0; i < ngaus; i++ ) {
    wacc[i+1] = wacc[i] + w[i];
  }

  hist.xmin = get_float("hist_xmin", 0.0);
  hist.xmax = get_float("hist_xmax", 20.0);
  hist.dx = get_float("hist_dx", 0.1);
  hist.n = Math.floor( (hist.xmax - hist.xmin) / hist.dx );
  hist.xmax = hist.xmin + hist.n * hist.dx;
  hist.arr = newarr(hist.n + 1);

  // compute the reference distribution
  hist.ref = newarr(hist.n);
  for ( i = 0; i < hist.n; i++ ) {
    x = hist.xmin + hist.dx * (i + 0.5);
    y = 0;
    for ( j = 0; j < ngaus; j++ ) {
      var dx = (x - mu[j]) / sig[j];
      y += w[j] / wacc[ngaus] * Math.sqrt(0.5/Math.PI) / sig[j]
           * Math.exp(-0.5 * dx * dx);
    }
    hist.ref[i] = y;
  }

  histplot = null;
}



/* generate random points */
function mkpoints()
{
  var ngaus = w.length, k, i, r, x;

  var npt = get_float("nptps", 100);

  var usetaus = (grab("RNG").value == "taus");

  for ( k = 0; k < npt; k++ ) {
    // select a random Gaussian mode
    // according to the fraction `w`
    if ( usetaus ) {
      r = taus_rand01();
    } else {
      r = rand01();
    }
    r *= wacc[ngaus];
    // find the corresponding bracket that contains r
    for ( i = 0; i < ngaus; i++ ) {
      if ( wacc[i + 1] >= r ) {
        break;
      }
    }

    // generate a point
    if ( usetaus ) {
      r = taus_randgaus();
    } else {
      r = randgaus();
    }
    x = mu[i] + sig[i] * r;

    // add the point to the histogram
    if ( x >= hist.xmin && x < hist.xmax ) {
      i = Math.floor( (x - hist.xmin) / hist.dx );
      hist.arr[i] += 1;
    }
  }
}



/* stop the previous timer, if any */
function stop_timer()
{
  if ( timer != null ) {
    window.clearInterval( timer );
    timer = null;
  }
}



function drawhist()
{
  // normalize the histogram
  var htot = 0;
  var h = newarr(hist.n);
  var i, imin = hist.n - 1, imax = 0, x, y;

  for ( i = 0; i < hist.n; i++ ) {
    y = hist.arr[i];
    if ( y <= 0 ) continue;
    if ( i > imax ) imax = i;
    if ( i < imin ) imin = i;
    htot += y;
  }

  var dat = "";
  for ( i = imin; i <= imax; i++ ) {
    x = hist.xmin + hist.dx * (i + 0.5);
    y = hist.arr[i] / (htot * hist.dx);
    dat += "" + x + "," + y + "," + hist.ref[i] + "\n";
  }

  if ( histplot == null ) {
    histplot = new Dygraph( grab("histplot"), dat, {
      xlabel: "<i>x</i>",
      ylabel: "Histogram",
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      width: 800,
      height: 480
    } );
  } else {
    histplot.updateOptions({ file: dat });
  }
}



function pulse()
{
  mkpoints();
  drawhist();
}



function start()
{
  stop_timer();
  loadparams();
  pulse();
  timer = setInterval( pulse, 1000 );
}



function stop()
{
  stop_timer();
}



function reset()
{
  stop_timer();
  loadparams();
}



function change_params()
{
  stop_timer();
  loadparams();
}



/* when a changes, change b as well */
function mapchange(a, b)
{
  grab(b).value = grab(a).value;
}



/* make a range control (i.e., slider) with a text box */
function mkrange(name, min, max, step, value, size)
{
  if ( size == null || size == undefined ) size = 8;
  // create a range control
  var s = '<input type="range" min="' + min + '" max="' + max +
    '" step="' + step + '" value="' + value + '" id="slider_' + name +
    '" class = "slider" onchange="mapchange(\'slider_' + name + '\', \'' + name + '\');change_params()">';
  // create a corresponding text box
  s += '<input type="text" value="' + value + '" id="' + name +
    '" size = "' + size + '" class="slider_num" ' +
    'onchange="mapchange(\'' + name + '\', \'slider_' + name + '\');change_params()">';
  return s;
}



function change_ngaus()
{
  stop_timer();

  var ngaus = grab("ngaus").value;
  if ( !is_int(ngaus) ) return;
  ngaus = parseInt(ngaus);
  var tab = grab("gausTable");
  var tbody = tab.lastChild;
  var rowid, row, td;

  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    rowid = "gaus_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = ngaus; i < nrows; i++ ) {
    rowid = "gaus_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create non-existing elements
  for ( var i = nrows; i < ngaus; i++ ) {
    var i1 = "" + (i+1);
    rowid = "gaus_row_" + i1;
    row = document.createElement("tr");

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // mu
    td = document.createElement("td");
    td.innerHTML = mkrange("mu_" + i1, 0, 20, 0.1, 5.0);
    row.appendChild(td);

    // sig
    td = document.createElement("td");
    td.innerHTML = mkrange("sig_"+i1, 0, 10, 0.1, 1.0);
    row.appendChild(td);

    // w
    td = document.createElement("td");
    td.innerHTML = mkrange("w_"+i1, 0, 5.0, 0.01, 1.0);
    row.appendChild(td);
  }
}



function set_val2(name, val)
{
  grab(name).value = val;
  grab("slider_" + name).value = val;
}



/* set default Gaussian parameters */
function set_gparams(i, mu, sig, w)
{
  set_val2("mu_" + i, mu);
  set_val2("sig_" + i, sig);
  set_val2("w_" + i, w);
}



/* set up the default test case
 * initialize a few Gaussians */
function init_gauss()
{
  var ng = 2;
  set_val2("ngaus", ng);
  change_ngaus();
  set_gparams(1,  5.0,  1.0, 1.0);
  set_gparams(2,  10.0, 1.5, 1.0);
  reset();
}



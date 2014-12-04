
var datarr = [];
var color = [];



function rotate(x, y, th)
{
  var c = Math.cos(th), s = Math.sin(th);
  return [ c*x - s*y, s*x + c*y ];
}



function mkpoints()
{
  var ngaus = get_int("ngaus", 1);

  datarr = [];
  color = new Array(ngaus);
  for ( var i = 0; i < ngaus; i++ ) {
    var i1 = i + 1;
    var xc = get_float("xc_" + i1, 0.0);
    var yc = get_float("yc_" + i1, 0.0);
    var a = get_float("a_" + i1, 1.0);
    var b = get_float("b_" + i1, 1.0);
    var theta = get_float("theta_" + i1, 0.0);
    var th = theta * Math.PI / 180.0;
    var npt = get_int("npt_" + i1, 100);
    for ( var j = 0; j < npt; j++ ) {
      var xy = rotate(a*gaussrand(), b*gaussrand(), th);
      datarr.push( [xc + xy[0], yc + xy[1]] );
    }
    color[i] = random_color();
  }
}



function random_color(cmin, cmax)
{
  var x = Math.random() * 6;
  var i = Math.floor( x ), r = 0, g = 0, b = 0;

  if ( cmin == undefined || cmin == null ) cmin = 0;
  if ( cmax == undefined || cmax == null ) cmax = 255;
  var cvar = cmax - cmin + 1;
  x -= i;
  if ( i < 1 ) { // red to yellow
    r = cmax;
    g = cmin + Math.floor( cvar * x );
  } else if ( i < 2 ) { // yellow to green
    r = cmin + Math.floor( cvar * (1 - x) );
    g = cmax;
  } else if ( i < 3 ) { // green to cyan
    g = cmax;
    b = cmin + Math.floor( cvar * x );
  } else if ( i < 4 ) { // cyan to blue
    g = cmin + Math.floor( cvar * (1 - x) );
    b = cmax;
  } else if ( i < 5 ) { // blue to magenta
    b = cmax;
    r = cmin + Math.floor( cvar * x );
  } else {
    b = cmin + Math.floor( cvar * (1 - x) );
    r = cmax;
  }
  return [r, g, b];
}



/* approximately draw an ellipse */
function drawEllipse(ctx, width, height, xc, yc, w, h, th, s, rgb, fill)
{
  if ( fill ) {
    ctx.fillStyle = "rgba(" + rgb[0] + ", " + rgb[1] + ", " + rgb[2] + ", 0.2)";
  } else {
    ctx.strokeStyle = "rgb(" + rgb[0] + ", " + rgb[1] + ", " + rgb[2] + ")";
  }

  ctx.save();
  ctx.translate(width/2 + xc*s, height/2 + yc*s);
  ctx.rotate(th);
  ctx.scale(w, h);
  ctx.beginPath();
  ctx.arc(0, 0, s, 0, 2 * Math.PI, false);
  ctx.restore();

  // draw oval by Bezier curves
  //xc = width/2 + xc*s;
  //yc = height/2 + yc*s;
  //w = w*s;
  //h = h*s;
  //var a = rotate(0, -h/2, th);
  //var b = rotate(w*2/3, -h/2, th);
  //var c = rotate(w*2/3, +h/2, th);
  //ctx.beginPath();
  //ctx.moveTo(xc + a[0], yc + a[1]);
  //ctx.bezierCurveTo(xc + b[0], yc + b[1], xc + c[0], yc + c[1], xc - a[0], yc - a[1]);
  //ctx.bezierCurveTo(xc - b[0], yc - b[1], xc - c[0], yc - c[1], xc + a[0], yc + a[1]);

  //ctx.closePath();

  if ( fill ) {
    ctx.fill();
  } else {
    ctx.stroke();
  }
}



function draw_gauss(ret)
{
  var c = grab("demobox");
  var ctx = c.getContext("2d");
  var w = c.width;
  var h = c.height;
  var r = Math.min(w/2, h/2) - 1;

  // draw the background
  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect(0, 0, w, h);

  var s = 60.0; // aspect ratio from real to screen coordinates
  var ngaus = get_int("ngaus");

  for ( var i = 0; i < ngaus; i++ ) {
    var i1 = i + 1;
    var xc = get_float("xc_" + i1, 0.0);
    var yc = get_float("yc_" + i1, 0.0);
    var a = get_float("a_" + i1, 1.0);
    var b = get_float("b_" + i1, 1.0);
    var theta = get_float("theta_" + i1, 0.0);
    for ( var l = 3; l > 0; l-- )
      drawEllipse(ctx, w, h, xc, yc, l*a, l*b, theta * Math.PI/180, s,
          color[i], true);
  }

  // draw random dots
  ctx.fillStyle = "black";

  for ( var i = 0; i < datarr.length; i++ ) {
    var x = datarr[i][0], y = datarr[i][1];
    // a dot === a one-by-one rectangle
    ctx.fillRect( w/2 + x * s, h/2 + y * s, 1, 1);
  }

  // draw the clustering result
  if ( ret ) {
    var l = 2.0; // 2*Math.log(2);
    for ( var i = 0; i < ret.length; i++ ) {
      var xc = ret[i].xc;
      var yc = ret[i].yc;
      var a = ret[i].a;
      var b = ret[i].b;
      var theta = ret[i].theta;
      drawEllipse(ctx, w, h, xc, yc, l*a, l*b, theta, s,
          [0,0,0], false);
    }
  }
}



/* return the two eigenvalues of the 2x2 covariance matrix m
 * | m[0]  m[1] |
 * | m[2]  m[3] | */
function get_abth(m)
{
  /* the eigenvalues satisfies
   * (x - m[0])(x - m[3]) - m[1] m[2] = 0
   * [x - (m[0]+m[3])/2]^2 = m[1] m[2] + (m[0] - m[3])^2/4 */
  var tr = m[0] + m[3];
  var d = m[1]*m[2] + .25*(m[0]-m[3])*(m[0]-m[3]);
  d = ( d > 0 ) ? Math.sqrt(d) : 0;
  var a = Math.max(.5*tr + d, 0);
  var b = Math.max(.5*tr - d, 0);
  var e1 = m[0] - a, e2 = m[1];
  var th = Math.atan2(e1, -e2);
  return [Math.sqrt(a), Math.sqrt(b), th];
}



var timer = null;
var iter = 0;
var K;
var km;



/* stop the previous timer, if any */
function stop_timer()
{
  if ( timer != null ) {
    window.clearInterval( timer );
    timer = null;
  }
}



function doclus()
{
  var ret = new Array(K);
  var niter = get_int("niter");
  var niterps = get_int("niterps");
  var it = 0;

  if ( iter >= niter ) stop_timer();

  for ( it = 0; it < niterps && iter < niter; it++, iter++ ) {
    km.estep();
    km.mstep();
  }

  for ( var k = 0; k < K; k++ ) {
    abth = get_abth(km.var[k]);
    //console.log(k, km.av[k], km.var[k]);
    ret[k] = {
      xc: km.av[k][0],
      yc: km.av[k][1],
      a: abth[0],
      b: abth[1],
      theta: abth[2]
    };
  }
  draw_gauss(ret);

  //console.log(iter, ret);
  return ret;
}



function show()
{
  stop_timer();

  // generate Gaussians
  mkpoints();

  K = get_int("nclus");
  km = new Kmeans(2, K, datarr);
  iter = 0;
  timer = setInterval( doclus, 1000 );
}



function change_params()
{
  // if the timer is running, restart it
  if ( timer != null ) show();
}



/* when a changes, change b as well */
function mapchange(a, b)
{
  grab(b).value = grab(a).value;
}



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

    // x_c
    td = document.createElement("td");
    td.innerHTML = mkrange("xc_"+i1, -5, 5, 0.1, 0.0);
    row.appendChild(td);

    // y_c
    td = document.createElement("td");
    td.innerHTML = mkrange("yc_"+i1, -5, 5, 0.1, 0.0);
    row.appendChild(td);

    // a
    td = document.createElement("td");
    td.innerHTML = mkrange("a_"+i1, 0, 5, 0.1, 1.0);
    row.appendChild(td);

    // b
    td = document.createElement("td");
    td.innerHTML = mkrange("b_"+i1, 0, 5, 0.1, 1.0);
    row.appendChild(td);

    // theta
    td = document.createElement("td");
    td.innerHTML = mkrange("theta_"+i1, -90, 90, 1, 0);
    row.appendChild(td);

    // npt
    td = document.createElement("td");
    td.innerHTML = mkrange("npt_"+i1, 1, 10000, 1, 1000);
    row.appendChild(td);
  }
}



function set_val2(name, val)
{
  grab(name).value = val;
  grab("slider_" + name).value = val;
}



/* set default Gaussian parameters */
function set_gparams(i, xc, yc, a, b, theta, npt)
{
  set_val2("xc_" + i, xc);
  set_val2("yc_" + i, yc);
  set_val2("a_" + i, a);
  set_val2("b_" + i, b);
  set_val2("theta_" + i, theta);
  set_val2("npt_" + i, npt);
}



/* set up the default test case
 * initialize a few Gaussians */
function init_gauss()
{
  var ng = 4;
  set_val2("ngaus", ng);
  change_ngaus();
  set_gparams(1,  0.0,  0.0, 1.0, 1.0,   0, 1000);
  set_gparams(2,  1.5,  0.0, 0.2, 1.5,   0, 1000);
  set_gparams(3, -2.0,  2.0, 1.0, 0.1, -45, 1000);
  set_gparams(4, -2.0, -2.0, 1.0, 0.4, +45, 1000);
}



var lj = null;
var n = 55;
var rho = 0.7;
var tp = 1.5;
var rcdef = 1000.0;
var mddt = 0.002;
var thdt = 0.02;
var nstepsps = 100; // number of steps per second
var nstepspf = 10; // number of steps per frame
var timer_interval = 100; // in milliseconds
var mdtimer = null;



function getparams()
{
  n = get_int("n", 55);
  var dim = get_int("dimension", 2);
  if ( dim === 2 || dim === 3 ) {
    D = dim;
  }
  rho = get_float("density", 0.7);
  temp = get_float("temp", 1.5);
  rcdef = get_float("rcutoff", 1000.0);
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepsps = get_int("nstepspersec", 100);
  nstepspf = nstepsps * timer_interval / 1000;
}



function domd()
{
  var istep;

  for ( istep = 0; istep < nstepspf; istep++ ) {
    lj.vv(mddt);
    lj.vrescale(tp, thdt);
  }
}



function pulse()
{
  domd();
  if ( lj.dim === 2 ) {
    ljdraw2d(lj, "ljbox");
  } else if ( lj.dim === 3 ) {
    ljdraw3d(lj, "ljbox");
  }
}



function stopmd()
{
  if ( mdtimer != null ) {
    clearInterval(mdtimer);
    mdtimer = null;
    lj = null;
  }
  munit(viewmat);
}



// for the wheel event
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



function startmd()
{
  stopmd();
  getparams();
  lj = new LJ(n, D, rho, rcdef);
  installmouse();
  mdtimer = setInterval( function(){ pulse() }, timer_interval );
}


function changeparams()
{
  if ( mdtimer !== null ) startmd();
}

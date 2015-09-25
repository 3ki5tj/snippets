/* mouse control for rotating and scaling a 3D object
 * Note: onclick is not captured */




"use strict";



var mousedown = 0;
var mousemoved = 0;
var mousex = -1;
var mousey = -1;
var viewmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

var mousescale = 1.0;
var mousescaler = null;


function mousedown3d(e)
{
  e = e || window.event;
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown = 1;
  //console.log("mousedown", e.clientX, e.clientY, m2str(viewmat));
}



function mouseup3d(e)
{
  e = e || window.event;
  mousex = -1;
  mousey = -1;
  mousemoved = mousedown - 1;
  mousedown = 0;
  //console.log("mouseup", e.clientX, e.clientY, m2str(viewmat));
}



function mousemove3d(e)
{
  if ( !mousedown ) {
    return;
  }
  e = e || window.event;
  if ( mousex >= 0 && mousey >= 0 ) {
    var target = e.target ? e.target : e.srcElement;
    viewmat = mxrot3d(viewmat, 180.0 * (e.clientY - mousey) / target.height);
    viewmat = myrot3d(viewmat, 180.0 * (e.clientX - mousex) / target.width);
    paint(); // to be defined
  }
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown += 1;
  //console.log("mousemove", e.clientX, e.clientY, m2str(viewmat));
}



/* for the mouse wheel event */
function wheel3d(e)
{
  var delta = 0; // positive for scrolling up
  e = e || window.event;
  if ( e.wheelDelta ) { // IE/Opera
    delta = e.wheelDelta / 120;
  } else if ( e.detail ) { // Firefox
    delta = -e.detail / 3;
  }
  if ( delta > 0 ) {
    mousescale *= 1.05;
  } else if ( delta < 0 ) {
    mousescale *= 0.95;
  }

  if ( mousescaler ) {
    grab( mousescaler ).value = mousescale;
  }
  //console.log("wheel", delta);
  if ( e.preventDefault ) {
    e.preventDefault();
  }
  e.returnValue = false;
  paint(); // defined later
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



function installmouse( box, scaler )
{
  var target = grab( box );
  mousescaler = scaler;
  target.onmousedown = mousedown3d;
  target.onmouseup = mouseup3d;
  target.onmousemove = mousemove3d;
  installwheel(target, wheel3d);
}




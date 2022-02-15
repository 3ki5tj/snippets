/* Common routines */



"use strict";



// check if `n' is a valid integer
function is_int(n)
{
  return !isNaN( parseInt(n, 10) ) && isFinite(n);
}



// check if `n' is a valid floating number
function is_float(n)
{
  return !isNaN(parseFloat(n)) && isFinite(n);
}



function roundto(x, decimals)
{
  var t = Math.pow(10, decimals);
  return (Math.round(x*t)/t).toFixed(decimals);
}



function newarr(n)
{
  var i, a = new Array(n);
  for ( i = 0; i < n; i++ ) {
    a[i] = 0;
  }
  return a;
}



/* create a two-dimensional numeric array */
function newarr2d(m, n)
{
  var i, j, a = new Array(m);
  for ( i = 0; i < m; i++ ) {
    a[i] = newarr(n);
  }
  return a;
}



/* copy array */
function cparr(x, y, n)
{
  var i;

  for ( i = 0; i < n; i++ ) {
    x[i] = y[i];
  }
}



/* copy two-dimensional array */
function cparr2d(x, y, m, n)
{
  var j;

  for ( j = 0; j < m; j++ ) {
    cparr(x[j], y[j], n);
  }
}



function qSel(sel)
{
  let x = document.querySelector(sel);
  if ( x === null ) {
    console.log("cannot select element ", sel);
  }
  return x;
}


function getFloat(sel, def)
{
  var x = parseFloat( qSel(sel).value );
  return !isNaN(x) && isFinite(x) ? x : def;
}



function getInt(sel, def)
{
  var x = parseInt( qSel(sel).value, 10 );
  return !isNaN(x) && isFinite(x) ? x : def;
}





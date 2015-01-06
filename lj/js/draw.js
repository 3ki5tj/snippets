/* graphics routines */



"use strict";



/* draw a ball that is centered at (x, y) with radius r
 * color is the color of the ball
 * spotcolor is the color of the spotlight
 * the format of color is "#rrggbb" */
function drawBall(ctx, x, y, r, color, spotcolor,
    spotx, spoty, spotr)
{
  if ( spotcolor === undefined || spotcolor === null ) {
    spotcolor = "#a0a0a0";
  }
  if ( spotx === undefined || spotcolor === null ) {
    spotx = r * 0.3;
  }
  if ( spoty === undefined || spoty === null ) {
    spoty = r * 0.4;
  }
  if ( spotr === undefined || spotr === null ) {
    spotr = r * 0.1;
  }
  var grd = ctx.createRadialGradient(x + spotx, y - spoty, spotr, x, y, r);
  grd.addColorStop(0, spotcolor); // spotlight color
  grd.addColorStop(1, color); // ball color
  ctx.fillStyle = grd;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.fill();
}





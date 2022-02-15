


"use strict";



/* draw all atoms in the box */
function ljdraw2d(lj, c, userscale)
{
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (lj.l + 1);
  var radius = Math.floor( 0.5 * scale );

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ctx.fillStyle = "#d0d0d0";
  ctx.fillRect( -(lj.l + 1) * 0.5 * scale + width  * 0.5,
                -(lj.l + 1) * 0.5 * scale + height * 0.5,
                scale * (lj.l + 1), scale * (lj.l + 1));

  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect( -lj.l * 0.5 * scale + width  * 0.5,
                -lj.l * 0.5 * scale + height * 0.5,
                scale * lj.l, scale * lj.l);

  // draw each particle
  for (var i = 0; i < lj.n; i++) {
    var x = Math.floor(  (lj.x[i][0] - lj.l * 0.5) * scale + width  * 0.5 );
    var y = Math.floor( -(lj.x[i][1] - lj.l * 0.5) * scale + height * 0.5 );
    var spotcolor = "#abd";
    var color = "#27c";
    paintBall(ctx, x, y, radius, color, spotcolor);
  }
}



/* apply the view matrix */
function transform(x, l, viewmat)
{
  var n = x.length;
  var xyz = newarr2d(n, 3), xc = [l * 0.5, l * 0.5, l * 0.5], xi = [0, 0, 0];

  for ( var i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vmxv(xyz[i], viewmat, xi);
    vinc(xyz[i], xc);
  }
  return xyz;
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = [0, 0, 0];
  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }
  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    for ( j = i + 1; j < n; j++ ) {
      if ( xyz[j][2] < xyz[k][2] ) {
        l = k;
        k = j;
        j = l;
      }
    }
    if ( k != i ) {
      vcopy(rt, xyz[k]);
      vcopy(xyz[k], xyz[i]);
      vcopy(xyz[i], rt);
    }
  }
  return xyz;
}



/* draw all atoms in the box */
function ljdraw3d(lj, c, userscale, viewmat)
{
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (lj.l + 1.0);

  var xyz = transform(lj.x, lj.l, viewmat); // apply the rotation matrix
  xyz = sortbyz(xyz); // sort particles by the z order

  // draw each particle
  var zmax = xyz[lj.n - 1][2], zmin = xyz[0][2];
  for (var i = 0; i < lj.n; i++) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var scl = scale * (0.7 + 0.3 * zf);
    var x = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scl + width  * 0.5 );
    var y = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scl + height * 0.5 );
    var spotcolor = rgb2str(100 + 100 * zf, 120 + 100 * zf, 150 + 105 * zf);
    var color = rgb2str(10, 60 + 60*zf, 90 + 120 * zf);
    var rz = Math.floor( 0.5 * scl );
    paintBall(ctx, x, y, rz, color, spotcolor);
  }
}




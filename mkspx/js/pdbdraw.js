


"use strict";



// NOTE: assuming the coordinates have been centered
function transform3d(x)
{
  var i, n = x.length, xyz = [];

  // rotate the coordinates of each particle
  for ( i = 0; i < n; i++ )
    xyz.push( [vdot(viewmat[0], x[i]), vdot(viewmat[1], x[i]), vdot(viewmat[2], x[i])] );
  return xyz;
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = new Array(n), rt = new Array(D);
  var idmap = new Array(n);
  var invmap = new Array(n);

  for ( i = 0; i < n; i++ ) {
    xyz[i] = new Array(3);
    idmap[i] = i;
    // i:         index of the output array `xyz`
    // idmap[i]:  index of the input array `x`
    // so xyz[i] --> x[ idmap[i] ];
    invmap[i] = i;
  }

  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }

  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    var zmin = x[ idmap[i] ][2];
    for ( j = i + 1; j < n; j++ ) {
      if ( x[ idmap[j] ][2] < zmin ) {
        k = j;
        zmin = x[ idmap[j] ][2];
      }
    }
    if ( k != i ) {
      // before
      //  xyz[i] --> x[ idmap[i] ]
      //  xyz[k] --> x[ idmap[k] ]
      // after
      //  xyz[i] --> x[ idmap[k] ]
      //  xyz[k] --> x[ idmap[i] ]
      l = idmap[i];
      idmap[i] = idmap[k];
      idmap[k] = l;
    }
  }

  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[ idmap[i] ]);
  }
  // compute the inverse map
  for ( i = 0; i < n; i++ ) {
    invmap[ idmap[i] ] = i;
  }
  return [xyz, idmap, invmap];
}



// get the scaling factor due to z
function getzscale(r, zmin, zmax, ortho)
{
  if ( ortho ) {
    return 0.9;
  } else {
    var zf = (r[2] - zmin) / (zmax - zmin);
    return 0.7 + 0.3 * zf;
  }
}



function getContactPoint(xi, xj, radius)
{
  var rji, xji = new Array(D);

  rji = vnorm( vdiff(xji, xj, xi) );
  vsmul(xji,  radius / rji);
  return vinc(xji, xi);
}

var atom_dict = {
  "H": [0.31, "#cccccc"],
  "N": [0.71, "#000080"],
  "C": [0.75, "#208080"],
  "O": [0.66, "#cc2020"],
  "S": [1.05, "#808020"]
};

// draw all atoms in the box
function pdbdraw(seq, atomls, l, boxsize,
    target, userscale, ballscale,
    overwrite, grey)
{
  var c = document.getElementById(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var n = atomls.length;
  var i, j, jb, k, ir, ic, ret;

  if ( !ballscale ) {
    ballscale = 1.0;
  }

  if ( !overwrite ) {
    // draw the background
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);
  }

  var x = [];
  for ( i = 0; i < n; i++ ) x.push( atomls[i][2] );
  var xt = transform3d(x); // apply the rotation matrix
  ret = sortbyz(xt); // sort particles by the z order
  var xyz = ret[0];
  var idmap = ret[1], invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]

  var ortho = document.getElementById("orthographic").checked;
  var scale = userscale * Math.min(width, height) / (l * 3.0);
  var zmax = xyz[n - 1][2], zmin = xyz[0][2];

  // draw the box
  if ( boxsize > 0 ) {
    var bh = 0.5 * boxsize;
    var box = transform3d([
        [-bh, -bh, -bh], [-bh, -bh, bh], [-bh,  bh, -bh], [-bh,  bh, bh],
        [ bh, -bh, -bh], [ bh, -bh, bh], [ bh,  bh, -bh], [ bh,  bh, bh] ]);
    ctx.strokeStyle = "#cccccc";
    ctx.lineWidth = "1px";
    for ( i = 0; i < 8; i++ ) {
      for ( j = i + 1; j < 8; j++ ) {
        if ( vdist(box[i], box[j]) > bh*2.1 ) continue;
        var scli = scale * getzscale(box[i], zmin, zmax, ortho);
        var xi = Math.floor(  box[i][0] * scli + width  * 0.5 );
        var yi = Math.floor( -box[i][1] * scli + height * 0.5 );
        var sclj = scale * getzscale(box[j], zmin, zmax, ortho);
        var xj = Math.floor(  box[j][0] * sclj + width  * 0.5 );
        var yj = Math.floor( -box[j][1] * sclj + height * 0.5 );
        drawLine(ctx, xi, yi, xj, yj);
      }
    }
  }

  // draw the legend
  var xo = [80, 80, 0], len = 60;
  var dx = [0, 0, 0], dy = [0, 0, 0], dz = [0, 0, 0];
  var xx = [0, 0, 0], xy = [0, 0, 0], xz = [0, 0, 0];
  var tx = [0, 0, 0], ty = [0, 0, 0], tz = [0, 0, 0];
  var colorx = "#cc2000", colory = "#20cc00", colorz = "#0020cc";
  vmxv(dx, viewmat, [1, 0, 0]);
  vmxv(dy, viewmat, [0, 1, 0]);
  vmxv(dz, viewmat, [0, 0, 1]);
  vsadd(xx, xo, dx, len);
  vsadd(xy, xo, dy, len);
  vsadd(xz, xo, dz, len);
  drawLine(ctx, xo[0], height - xo[1], xx[0], height - xx[1], colorx, "2px");
  drawLine(ctx, xo[0], height - xo[1], xy[0], height - xy[1], colory, "2px");
  drawLine(ctx, xo[0], height - xo[1], xz[0], height - xz[1], colorz, "2px");

  // draw each particle
  for (i = 0; i < n; i++) {

    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    var i0 = idmap[ i ];
    //color = rgb2str(160 + 50 * zf, 160 + 50 * zf, 160 + 50 * zf);
    var atom = atomls[i0];
    var elem = atom[0].slice(0, 1);
    //console.log(i0, elem, atom);
    var fullColor = atom_dict[elem][1];
    var color = darkenColor( fullColor, 0.8 + 0.2 * zf );
    if ( grey ) {
      color = lightenColor( greyColor( color ), 0.5 );
    }
    var spotcolor = lightenColor( color, 0.3 );
    // make closer particles larger
    var scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
    var xi = Math.floor(  xyz[i][0] * scli + width  * 0.5 );
    var yi = Math.floor( -xyz[i][1] * scli + height * 0.5 );
    var rad = atom_dict[elem][0] * ballscale;
    var rz = Math.floor( rad * scli );
    paintBall(ctx, xi, yi, rz, color, spotcolor);

    // collect bonded neighbors
    var j0, j1, jarr = [], jj;
    j0 = i0 - 20; if ( j0 < 0 ) j0 = 0;
    j1 = i0 + 20; if ( j1 > n ) j1 = n;
    var resi = atomls[i0][1];
    for ( jj = j0; jj < j1; jj++ ) {
      if ( jj === i0 ) continue;
      var resj = atomls[jj][1];
      var elemj = atomls[jj][0].slice(0, 1);
      if ( (resj === resi)
        || (resj === resi + 1 && elem === "C" && elemj === "N")
        || (resj === resi - 1 && elem === "N" && elemj === "C") ) {
        var dis = vdist(x[jj], x[i0]);
        if ( dis < 1.9 ) jarr.push( jj );
      }
    }

    // draw bonded neighbors
    for ( jj = 0; jj < jarr.length; jj++ ) {
      j0 = jarr[jj];
      j = invmap[ j0 ];
      if ( j > i ) {
        var ri = getContactPoint(xyz[i], xyz[j], rad);
        xi = Math.floor(  ri[0] * scli + width  * 0.5 );
        yi = Math.floor( -ri[1] * scli + height * 0.5 );

        var sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
        var xj = Math.floor(  xyz[j][0] * sclj + width  * 0.5 );
        var yj = Math.floor( -xyz[j][1] * sclj + height * 0.5 );

        drawLineGradient(ctx, xi, yi, xj, yj);
      }
    }
  }
  var lent = len + 20;
  vsadd(tx, xo, dx, lent);
  vsadd(ty, xo, dy, lent);
  vsadd(tz, xo, dz, lent);
  ctx.font = "italic 20px Times New Roman";
  ctx.fillStyle = colorx;
  ctx.fillText("x", tx[0], height - tx[1]);
  ctx.fillStyle = colory;
  ctx.fillText("y", ty[0], height - ty[1]);
  ctx.fillStyle = colorz;
  ctx.fillText("z", tz[0], height - tz[1]);
}


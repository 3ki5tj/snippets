


"use strict";



function transform(x)
{
  var i, d, n = x.length, l = 0;
  var xyz = newarr2d(n, 3), xc = [0, 0, 0], xi = [0, 0, 0];

  // compute the center of mass
  for ( i = 0; i < n; i++ ) {
    vinc(xc, x[i]);
  }
  vsmul(xc, 1.0/n);

  // rotate the coordinates of each particle
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vmxv(xyz[i], viewmat, xi);
    //console.log(x[i], xi, xc, xyz[i]);
    //vinc(xyz[i], xc);
    for ( d = 0; d < D; d++ ) {
      l = Math.max( Math.abs( xi[d] ), l );
    }
  }
  return [xyz, xc, l];
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = newarr(D);
  var idmap = newarr(n);
  var invmap = newarr(n);

  for ( i = 0; i < n; i++ ) {
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
    return 0.8 + 0.2 * zf;
  }
}



// draw all atoms in the box
function nadraw(na, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var i, j, jb, k, ir, ic, ret;
  var apr = na.apr;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ret = transform(na.x); // apply the rotation matrix
  var xt = ret[0];
  //na.l = ret[2];
  na.l = 8.0 * Math.pow( na.nr, 1.0/3 );
  ret = sortbyz(xt); // sort particles by the z order
  var xyz = ret[0];
  var idmap = ret[1], invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]

  var ortho = grab("orthographic").checked;
  var scale = userscale * Math.min(width, height) / (na.l * 2.0);

  // draw each particle
  var zmax = xyz[na.n - 1][2], zmin = xyz[0][2];

  // draw lines that were used to group clusters
  {
    ctx.lineWidth = 2;
    ctx.strokeStyle = '#808080';
    for ( ir = 0; ir < na.nr; ir++ ) {
      i = invmap[ ir*apr ];
      j = invmap[ ir*apr + 1 ];
      if ( apr == 3 ) {
        jb = invmap[ ir*apr + 2 ];
      }

      var scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
      var xi = Math.floor(  xyz[i][0] * scli + width  * 0.5 );
      var yi = Math.floor( -xyz[i][1] * scli + height * 0.5 );

      var sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
      var xj = Math.floor(  xyz[j][0] * sclj + width  * 0.5 );
      var yj = Math.floor( -xyz[j][1] * sclj + height * 0.5 );

      var scljb = scale * getzscale(xyz[jb], zmin, zmax, ortho);
      var xjb = Math.floor(  xyz[jb][0] * scljb + width  * 0.5 );
      var yjb = Math.floor( -xyz[jb][1] * scljb + height * 0.5 );

      if ( apr === 2 ) {
        drawLineGradient(ctx, xi, yi, xj, yj);
      } else {
        drawLineGradient(ctx, xi, yi, xj, yj);
        drawLineGradient(ctx, xj, yj, xjb, yjb);
      }

      if ( ir < na.nr - 1 ) {
        k = invmap[ (ir + 1) * apr ];
        var sclk = scale * getzscale(xyz[k], zmin, zmax, ortho);
        var xk = Math.floor(  xyz[k][0] * sclk + width  * 0.5 );
        var yk = Math.floor( -xyz[k][1] * sclk + height * 0.5 );
        if ( apr === 2 ) {
          drawLineGradient(ctx, xi, yi, xk, yk);
        } else {
          drawLineGradient(ctx, xj, yj, xk, yk);
        }
      }
    }
  }


  for (i = 0; i < na.n; i++) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var scl = scale * getzscale(xyz[i], zmin, zmax, ortho);
    var x = Math.floor(  xyz[i][0] * scl + width  * 0.5 );
    var y = Math.floor( -xyz[i][1] * scl + height * 0.5 );
    var color, rad;
    var i0 = idmap[ i ];
    var atype = i0 % apr;
    var ir = i0 / apr;
    var ra = na.seq.substr(ir, 1);
    if ( atype === 0 ) {
      color = rgb2str(40 + 20 * zf, 40 + 20 * zf, 40 + 20 * zf);
      rad = 1.6;
    } else if ( atype === apr - 1 ) {
      if ( ra == 'A' ) {
        color = rgb2str(120 + 60 * zf, 32, 20);
      } else if ( ra == 'C' ) {
        color = rgb2str(20, 80 + 60 * zf, 20);
      } else if ( ra == 'G' ) {
        color = rgb2str(20, 32, 80 + 60 * zf);
      } else if ( ra == 'U' ) {
        color = rgb2str(80 + 60 * zf, 80 + 60 * zf, 20);
      }
      rad = 1.6;
    } else {
      color = rgb2str(160 + 50 * zf, 160 + 50 * zf, 160 + 50 * zf);
      rad = 1.6;
    }
    var spotcolor = lightenColor(color, 0.3);
    var rz = Math.floor( rad * scl );
    paintBall(ctx, x, y, rz, color, spotcolor);
  }
}





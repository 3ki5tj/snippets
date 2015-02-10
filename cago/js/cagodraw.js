


"use strict";



function transform3d(x)
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
function cagodraw(go, target, userscale)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var i, j, jb, k, ir, ic, ret;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  ret = transform3d(go.x); // apply the rotation matrix
  var xt = ret[0];
  //go.l = ret[2];
  go.l = 5.0 * Math.pow( go.n, 1.0/3 );
  ret = sortbyz(xt); // sort particles by the z order
  var xyz = ret[0];
  var idmap = ret[1], invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]

  var ortho = grab("orthographic").checked;
  var scale = userscale * Math.min(width, height) / (go.l * 2.0);

  // draw each particle
  var zmax = xyz[go.n - 1][2], zmin = xyz[0][2];

  for (i = 0; i < go.n; i++) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    var i0 = idmap[ i ];
    //color = rgb2str(160 + 50 * zf, 160 + 50 * zf, 160 + 50 * zf);
    var iaa = go.iaa[ i0 ];
    var color = darkenColor( aacolors[ iaa ], 0.8 + 0.2 * zf );
    var spotcolor = lightenColor( aacolors[ iaa ], 0.7 - 0.4 * zf );
    // make closer particles larger
    var scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
    var xi = Math.floor(  xyz[i][0] * scli + width  * 0.5 );
    var yi = Math.floor( -xyz[i][1] * scli + height * 0.5 );
    var rad = aaradii[ iaa ];
    var rz = Math.floor( rad * scli );
    var xj, yj, sclj;

    paintBall(ctx, xi, yi, rz, color, spotcolor);

    // draw bonds to the adjacent residues
    var j0, jarr = [], jj;
    if ( i0 < go.n - 1 ) { jarr.push( i0 + 1 ); }
    if ( i0 >= 1 ) { jarr.push( i0 - 1 ); }

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
}




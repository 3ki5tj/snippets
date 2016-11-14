


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
  var xyz = newarr2d(n, 3), rt = new Array(D);
  var idmap = new Array(n);
  var invmap = new Array(n);

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
function pdbdraw(seq, x, atomls, l,
    target, userscale, ballscale,
    overwrite, grey)
{
  var c = document.getElementById(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var n = x.length;
  var i, j, jb, k, ir, ic, ret;

  if ( !ballscale ) {
    ballscale = 1.0;
  }

  if ( !overwrite ) {
    // draw the background
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);
  }

  ret = transform3d(x); // apply the rotation matrix
  var xt = ret[0];
  ret = sortbyz(xt); // sort particles by the z order
  var xyz = ret[0];
  var idmap = ret[1], invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]

  var ortho = document.getElementById("orthographic").checked;
  var scale = userscale * Math.min(width, height) / (l * 2.5);

  // draw each particle
  var zmax = xyz[n - 1][2], zmin = xyz[0][2];

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
}


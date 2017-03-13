"use strict";

function erfc(z)
{
  var t = 2./(2. + z), ty = 4*t - 2, tmp, d = 0, dd = 0;
  var c = [-1.3026537197817094, 6.4196979235649026e-1,
    1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4,
    3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5,
    -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
    6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10,
    9.6467911e-11, 2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13,
    -1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17];
  for ( var j = 27; j > 0; j-- ) {
    tmp = d;
    d = ty * d - dd + c[j];
    dd = tmp;
  }
  return t * Math.exp(-z * z + 0.5 * (c[0] + ty * d) - dd);
}

function erf(z)
{
  return 1 - erfc(z);
}

var echarge = 1.6021766208e-19;
var eps0 = 8.854187817e-12;
var NA = 6.022140857e23;
var cal = 4.184;
var atm = 101325;

/* vector cross product */
function vcross(a, b)
{
  return [ a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0] ];
}

function vscale(v, s)
{
  return [ v[0] * s, v[1] * s, v[2] * s ];
}

function getvol(mat)
{
  return Math.abs(mat[0][0] * mat[1][1] * mat[2][2]);
}

function getmindim(mat)
{
  var x = Math.abs( mat[0][0] ), y = Math.abs( mat[1][1] ), z = Math.abs( mat[2][2] );
  return Math.min(x, y, z);
}

function inverse(mat)
{
  var invmat = new Array(3), invvol = 1. / getvol(mat);
  invmat[0] = vcross(mat[1], mat[2]); vsmul(invmat[0], invvol);
  invmat[1] = vcross(mat[2], mat[0]); vsmul(invmat[1], invvol);
  invmat[2] = vcross(mat[0], mat[1]); vsmul(invmat[2], invvol);
  return invmat;
}

function getdist2(x, y, z, mat)
{
  var v = [x, y, z], u = [0, 0, 0], dis2 = 0;
  for ( var i = 0; i < 3; i++ ) {
    for ( var j = 0; j < 3; j++ )
      u[i] += v[j] * mat[j][i];
    dis2 += u[i] * u[i];
  }
  return dis2;
}

// NOTE: assuming the coordinates have been centered
function transform3d(x)
{
  var i, n = x.length, xyz = [];

  // rotate the coordinates of each particle
  for ( i = 0; i < n; i++ )
    xyz.push( [vdot(viewmat[0], x[i]), vdot(viewmat[1], x[i]), vdot(viewmat[2], x[i])] );
  return xyz;
}

function drawbox(mat, target, userscale)
{
  var c = document.getElementById(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;
  var i, j, jb, k, ir, ic, ret;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  //var ortho = document.getElementById("orthographic").checked;
  var l0 = Math.abs(mat[0][0]), l1 = Math.abs(mat[1][1]), l2 = Math.abs(mat[2][2]);
  var l = Math.max(l0, l1, l2);
  var ll = Math.sqrt(l0*l0 + l1*l1 + l2*l2);
  var scale = userscale * Math.min(width, height) / (1.1*ll);

  // draw the box
  var bh = 0.5 * l;
  var x0 = [
    [-0.5, -0.5, -0.5], [+0.5, -0.5, -0.5],
    [-0.5, +0.5, -0.5], [+0.5, +0.5, -0.5],
    [-0.5, -0.5, +0.5], [+0.5, -0.5, +0.5],
    [-0.5, +0.5, +0.5], [+0.5, +0.5, +0.5]];
  var x1 = new Array(8), v0, v1, v2;
  for ( i = 0; i < 8; i++ ) {
    v0 = vscale(mat[0], x0[i][0]);
    v1 = vscale(mat[1], x0[i][1]);
    v2 = vscale(mat[2], x0[i][2]);
    x1[i] = [ v0[0] + v1[0] + v2[0],
              v0[1] + v1[1] + v2[1],
              v0[2] + v1[2] + v2[2] ];
  }

  var box = transform3d(x1);
  ctx.strokeStyle = "#cccccc";
  ctx.lineWidth = "1px";
  var pr = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6],
            [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
  var colorxyz = ["#cc2000", "#20cc00", "#0020cc"];
  for ( k = 0; k < pr.length; k++ ) {
    i = pr[k][0];
    j = pr[k][1];
    var xi = Math.floor(  box[i][0] * scale + width  * 0.5 );
    var yi = Math.floor( -box[i][1] * scale + height * 0.5 );
    var xj = Math.floor(  box[j][0] * scale + width  * 0.5 );
    var yj = Math.floor( -box[j][1] * scale + height * 0.5 );
    var color = "#cccccc", lw = 2;
    if ( k < 3 ) {
      color = colorxyz[k];
      lw = 4;
    }
    drawLine(ctx, xi, yi, xj, yj, color, lw);
  }

  var sxyz = ["x", "y", "z"];
  ctx.font = "bold 24px Times New Roman, serif";
  for ( k = 0; k < 3; k++ ) {
    var tx0 = [0, 0, 0], tx1, tx = [0, 0, 0], vv = [0, 0, 0];
    vnormalize( vcopy(vv, mat[k]) );
    vsadd(tx0, x1[1 << k], vv, 0.07 * l);
    tx1 = [vdot(viewmat[0], tx0),
           vdot(viewmat[1], tx0),
           vdot(viewmat[2], tx0)];
    tx[0] = Math.floor(  tx1[0] * scale + width  * 0.5 );
    tx[1] = Math.floor( -tx1[1] * scale + height * 0.5 );
    ctx.fillStyle = colorxyz[k];
    ctx.fillText(sxyz[k], tx[0], tx[1]);
  }
}


function ewald(kappa, tol, mat)
{
  var i, j, l, xm = [0, 0, 0], km = [0, 0, 0];
  var r, k2, eself = 0, ereal = 0, erecip = 0, ebg = 0;
  var vol = getvol(mat), dis2;

  var dim = getmindim(mat);
  var invmat = inverse(mat);
  var invdim = getmindim(invmat);

  var sqrta = (kappa > 0) ? kappa : Math.sqrt(Math.PI*invdim/dim);
  var inva = Math.PI * Math.PI / (sqrta * sqrta);

  //console.log(dim, invmat[0], invmat[1], invmat[2], vol, invdim, kappa, sqrta, inva);

  // estimate the number of neighboring cells
  for ( i = 0; i < 3; i++ ) {
    for ( xm[i] = 1; xm[i] <= 1000; xm[i]++ ) {
      r = Math.abs(mat[i][i]) * xm[i];
      if ( erfc(sqrta*r)/r < tol ) break;
    }
    //console.log("dim " + i + ", xm " + xm[i] + ", sqrta " + sqrta + ", r " + r + " err " + erfc(sqrta*r)/r);
  }

  // estimate the number of wave vectors
  for ( i = 0; i < 3; i++ ) {
    for ( km[i] = 1; km[i] <= 1000; km[i]++ ) {
      r = Math.abs(invmat[i][i]) * km[i];
      k2 = r * r;
      if ( Math.exp(-k2*inva)/k2/Math.PI/vol < tol ) break;
    }
  }

  // real-space sum
  for ( i = -xm[0]; i <= xm[0]; i++ ) {
    for ( j = -xm[1]; j <= xm[1]; j++ ) {
      for ( l = -xm[2]; l <= xm[2]; l++ ) {
        dis2 = getdist2(i, j, l, mat);
        if ( dis2 <= 0 ) continue;
        r = Math.sqrt(dis2);
        ereal += erfc(sqrta*r)/r;
      }
    }
  }
  ereal *= 0.5;

  // reciprocal-space sum
  for ( i = -km[0]; i <= km[0]; i++ ) {
    for ( j = -km[1]; j <= km[1]; j++ ) {
      for ( l = -km[2]; l <= km[2]; l++ ) {
        k2 = getdist2(i, j, l, invmat);
        if ( k2 <= 0 ) continue;
        erecip += Math.exp(-k2*inva)/k2;
      }
    }
  }
  erecip *= 0.5 / (Math.PI * vol);

  // self energy
  eself = -sqrta / Math.sqrt(Math.PI);

  // background energy
  ebg = -0.5 * Math.PI / (sqrta * sqrta * vol);
  //console.log(xm, km, ereal, erecip, eself, ebg);
  return [ereal, erecip, eself, ebg, sqrta, xm, km];
}

// format number in HTML format
function numformat(x, pres)
{
  if ( !pres ) pres = 5;
  var s = x.toPrecision(pres), i;
  if ( s.substring(0, 1) === "-" ) {
    s = "&minus;" + s.substring(1);
  }
  var p = s.indexOf("e");
  if ( p >= 0 ) { // scientific format
    var s1 = s.substring(0, p);
    // remove the trailing zeros
    for ( i = s1.length; i > 0; i-- )
      if ( s1.substring(i - 1, i) !== "0" )
        break;
    s1 = s1.substring(0, i);
    // remove the final "."
    if ( s1.substring(i - 1, i) === "." )
      s1 = s1.substring(0, i - 1);
    // the exponent
    var s2 = s.substring(p + 1);
    if ( s2.substring(0, 1) === "+" ) { // remove the leading `+` sign
      s2 = s2.substring(1);
    }
    if ( s2.substring(0, 1) === "-" ) { // convert the leading `-` sign to &minus;
      s2 = "&minus;" + s2.substring(1);
    }
    s = s1 + "&times;10<sup>" + s2 + "</sup>";
  } else {
    for ( i = s.length; i > 0; i-- )
      if ( s.substring(i - 1, i) !== "0" )
        break;
    s = s.substring(0, i);
    // remove the final "."
    if ( s.substring(i - 1, i) === "." )
      s = s.substring(0, i - 1);
  }
  return s;
}

function readPDB()
{
  var pdb = document.getElementById("inputPDB").value.trim();
  if ( pdb === "" ) return;
  var lines = pdb.split("\n"), i, line;
  var resid = -1, resi, resnm;
  var qtot = 0;
  for ( i = 0; i < lines.length; i++ ) {
    line = lines[i].trim();
    if ( line === "" || line.substring(0, 6) !== "ATOM  " )
      continue;
    resi = parseInt( line.substring(22, 26).trim() );
    if ( resi === resid ) { // old residue
      continue;
    } else {
      resid = resi;
    }
    resnm = line.substring(17, 21).trim().toUpperCase();
    if ( resnm === "SOD" || resnm === "POT" || resnm === "CES"
      || resnm === "NA"  || resnm === "K" || resnm === "CS"
      || resnm === "ARG" || resnm === "LYS" || resnm === "LYSP"
      || resnm === "HSP" || resnm === "HIP" || resnm === "HISP"
      || resnm === "PROP" ) {
      qtot += 1;
    } else if ( resnm === "CAL" || resnm === "MG" || resnm === "ZN2" ) {
      qtot += 2;
    } else if ( resnm === "CLA" || resnm === "CL"
      || resnm === "ASP" || resnm === "GLU"
      || resnm === "GUA" || resnm === "ADE" || resnm === "CYT"
      || resnm === "THY" || resnm === "URA"
      || resnm === "CTER" ) {
      qtot += -1;
    } else if ( resnm === "DISU" ) {
      qtot += -0.36;
    }
  }
  document.getElementById("qtot").value = qtot;
}

function changelunit()
{
  var lunit = document.getElementById("lunit").value;
  if ( lunit !== "" ) {
    document.getElementById("invl_1").innerHTML = " (" + lunit + "<sup>&minus;1</sup>)";
    document.getElementById("l_0").innerHTML = " (" + lunit + ")";
    document.getElementById("l_1").innerHTML = " (" + lunit + ")";
    document.getElementById("l_2").innerHTML = " (" + lunit + ")";
    document.getElementById("l_3").innerHTML = " (" + lunit + ")";
  }
  return lunit;
}

function readNAMDconf()
{
  var conf = document.getElementById("NAMDconf").value.trim();
  if ( conf === "" ) return;
  document.getElementById("eunit").selectedIndex = 1;
  document.getElementById("punit").selectedIndex = 1;
  document.getElementById("lunit").selectedIndex = 1;
  changelunit();
  var cutoff = 12, pmetolerance = 1e-6;
  var lines = conf.split("\n"), i, line;
  var x = [0, 0, 0], y = [0, 0, 0], z = [0, 0, 0];
  var haslattice = false;
  var alchDecouple = false;
  for ( i = 0; i < lines.length; i++ ) {
    line = lines[i].trim();
    if ( line === "" || line.charAt(0) === "#" )
      continue;
    var arr = line.split(/\s+/);
    var key = arr[0].toLowerCase();
    if ( key === "cellbasisvector1" ) {
      x[0] = parseFloat( arr[1] );
      x[1] = parseFloat( arr[2] );
      x[2] = parseFloat( arr[3] );
      haslattice = true;
    } else if ( key === "cellbasisvector2" ) {
      y[0] = parseFloat( arr[1] );
      y[1] = parseFloat( arr[2] );
      y[2] = parseFloat( arr[3] );
      haslattice = true;
    } else if ( key === "cellbasisvector3" ) {
      z[0] = parseFloat( arr[1] );
      z[1] = parseFloat( arr[2] );
      z[2] = parseFloat( arr[3] );
      haslattice = true;
    } else if ( key === "cutoff" ) {
      cutoff = parseFloat( arr[1] );
    } else if ( key === "pmetolerance" ) {
      pmetolerance = parseFloat( arr[1] );
    } else if ( key === "alchdecouple" ) {
      var val = arr[1].toLowerCase();
      alchDecouple = ( val === "on" || val === "yes" || val === "true" );
    }
  }

  if ( haslattice ) {
    // update the lattice type
    if ( x[1] != 0 || x[2] != 0 || y[0] != 0 || y[2] != 0 || z[0] != 0 || z[1] != 0 ) {
      document.getElementById("lattype").selectedIndex = 2;
    } else if ( x[0] != y[1] || y[1] != z[2] ) {
      document.getElementById("lattype").selectedIndex = 1;
    } else {
      document.getElementById("lattype").selectedIndex = 0;
    }
    // update lattice parameters
    document.getElementById("x0").value = x[0];
    document.getElementById("x1").value = x[1];
    document.getElementById("x2").value = x[2];
    document.getElementById("y0").value = y[0];
    document.getElementById("y1").value = y[1];
    document.getElementById("y2").value = y[2];
    document.getElementById("z0").value = z[0];
    document.getElementById("z1").value = z[1];
    document.getElementById("z2").value = z[2];
  }

  // find the default ewaldcoef
  // code adapted from SimParameters.C
  var ewaldcof = 1.0;
  while ( erfc(ewaldcof*cutoff)/cutoff >= pmetolerance ) {
    ewaldcof *= 2.0;
  }
  var ewaldcof_lo = 0, ewaldcof_hi = ewaldcof;
  for ( i = 0; i < 100; ++i ) {
    ewaldcof = 0.5 * ( ewaldcof_lo + ewaldcof_hi );
    if ( erfc(ewaldcof*cutoff)/cutoff >= pmetolerance ) {
      ewaldcof_lo = ewaldcof;
    } else {
      ewaldcof_hi = ewaldcof;
    }
  }
  document.getElementById("kappa").value = ewaldcof;
  document.getElementById("alchDecouple").checked = alchDecouple;
  document.getElementById("hasvcorr").checked = false;
  document.getElementById("hasvreal").checked = false;
}

function update()
{
  var lattype = document.getElementById("lattype").value;
  var ix0 = document.getElementById("x0");
  var ix1 = document.getElementById("x1");
  var ix2 = document.getElementById("x2");
  var iy0 = document.getElementById("y0");
  var iy1 = document.getElementById("y1");
  var iy2 = document.getElementById("y2");
  var iz0 = document.getElementById("z0");
  var iz1 = document.getElementById("z1");
  var iz2 = document.getElementById("z2");
  if ( lattype === "cubic" ) {
    ix1.disabled = true; ix1.value = "0";
    ix2.disabled = true; ix2.value = "0";
    iy0.disabled = true; iy0.value = "0";
    iy2.disabled = true; iy2.value = "0";
    iz0.disabled = true; iz0.value = "0";
    iz1.disabled = true; iz1.value = "0";
    iy1.disabled = true; iy1.value = ix0.value;
    iz2.disabled = true; iz2.value = ix0.value;
  } else if ( lattype === "orthorhombic" ) {
    ix1.disabled = true; ix1.value = "0";
    ix2.disabled = true; ix2.value = "0";
    iy0.disabled = true; iy0.value = "0";
    iy2.disabled = true; iy2.value = "0";
    iz0.disabled = true; iz0.value = "0";
    iz1.disabled = true; iz1.value = "0";
    iy1.disabled = false;
    iz2.disabled = false;
  } else if ( lattype === "triclinic" ) {
    ix1.disabled = false;
    ix2.disabled = false;
    iy0.disabled = false;
    iy1.disabled = false;
    iy2.disabled = false;
    iz0.disabled = false;
    iz1.disabled = false;
    iz2.disabled = false;
  }
  var x0 = parseFloat( ix0.value );
  var x1 = parseFloat( ix1.value );
  var x2 = parseFloat( ix2.value );
  var y0 = parseFloat( iy0.value );
  var y1 = parseFloat( iy1.value );
  var y2 = parseFloat( iy2.value );
  var z0 = parseFloat( iz0.value );
  var z1 = parseFloat( iz1.value );
  var z2 = parseFloat( iz2.value );
  return [[x0, x1, x2], [y0, y1, y2], [z0, z1, z2]];
}

function paint()
{
  var mat = update();
  var vol = getvol(mat);
  if ( vol <= 0 ) {
    return;
  }
  mousescale = parseFloat( document.getElementById("scaleinput").value );
  drawbox(mat, "animationbox", mousescale);
  var kappa = parseFloat( document.getElementById("kappa").value );
  if ( isNaN(kappa) ) kappa = 0;
  var tol = parseFloat( document.getElementById("tol").value );
  if ( isNaN(tol) ) tol = 1e-15;
  var ene = ewald(kappa, tol, mat);
  var K0 = 1, P0 = 1, L0 = 1;
  //var lunit = document.getElementById("lunit").value;
  var lunit = changelunit();
  if ( lunit === "\u212B" ) {
    L0 = 1e-10;
  } else if ( lunit === "nm" ) {
    L0 = 1e-9;
  }
  //console.log(L0, lunit);
  var eunit = document.getElementById("eunit").value;
  var punit = document.getElementById("punit").value;
  if ( eunit === "kcal/mol" ) {
    K0 = echarge * echarge / (4 * Math.PI * eps0) * NA / (L0 * cal * 1e3);
    // 332.06371301869575
    // NAMD uses 332.0636
  } else if ( eunit === "kJ/mol" ) {
    K0 = echarge * echarge / (4 * Math.PI * eps0) * NA / (L0 * 1e3);
  } else if ( eunit === "J" ) {
    K0 = echarge * echarge / (4 * Math.PI * eps0) / L0;
  }
  var L0_4 = L0 * L0 * L0 * L0;
  if ( punit === "atm" ) {
    P0 = echarge * echarge / (4 * Math.PI * eps0) / (L0_4 * atm);
  } else if ( punit === "bar" ) {
    P0 = echarge * echarge / (4 * Math.PI * eps0) / (L0_4 * 1e5);
    // NAMD uses bar as the unit of pressure
    // PRESSUREFACTOR = cal/NA*1e30 = 100 * P0 / K0
    // NAMD's PRESSUREFACTOR is 6.95e6 instead of 6947695.345147256
  } else if ( punit === "Pa" ) {
    P0 = echarge * echarge / (4 * Math.PI * eps0) / L0_4;
  }
  var epsilon = parseFloat( document.getElementById("epsilon").value );
  if ( epsilon <= 0 ) epsilon = 1;
  var RBorn = parseFloat( document.getElementById("RBorn").value );
  var hasvcorr = document.getElementById("hasvcorr").checked;
  var hasvreal = document.getElementById("hasvreal").checked;
  var qtot = parseFloat( document.getElementById("qtot").value );
  var qq = qtot * qtot;
  var efin, evol, ecor, eborn, ecor2;
  var pfin, pvol, pcor, pcor2;
  efin = ene[0] + ene[1] + ene[2] + ene[3];
  evol = ene[3];
  eborn = 2 * Math.PI * RBorn * RBorn / (3 * vol);
  ecor = (efin + eborn) * (1 - 1/epsilon);
  ecor2 = ecor - efin;
  if ( !hasvcorr ) ecor2 += evol;
  if ( !hasvreal ) ecor2 += ene[0];
  pfin = efin / (3 * vol); // efin ~ 1/L ~ 1/V^(1/3)
  pvol = evol / (3 * vol); // evol ~ 1/(kappa^2 V) ~ 1/L
  pcor = (efin + 4 * eborn) / (3 * vol) * (1 - 1/epsilon);
  pcor2 = pcor - pfin;
  if ( !hasvcorr ) pcor2 += pvol;
  if ( !hasvreal ) pcor2 += ene[0] / (3 * vol);
  efin *= K0 * qq;
  evol *= K0 * qq;
  ecor *= K0 * qq;
  ecor2 *= K0 * qq;
  pfin *= P0 * qq;
  pvol *= P0 * qq;
  pcor *= P0 * qq;
  pcor2 *= P0 * qq;
  document.getElementById("efin").innerHTML = numformat(efin);
  document.getElementById("evol").innerHTML = numformat(evol);
  document.getElementById("ecor").innerHTML = numformat(ecor);
  document.getElementById("ecor2").innerHTML = numformat(ecor2);
  document.getElementById("pfin").innerHTML = numformat(pfin);
  document.getElementById("pvol").innerHTML = numformat(pvol);
  document.getElementById("pcor").innerHTML = numformat(pcor);
  document.getElementById("pcor2").innerHTML = numformat(pcor2);

  var hicolor = "#ffffcc", locolor = "white";
  if ( document.getElementById("alchDecouple").checked ) {
    document.getElementById("ecor").style.backgroundColor = hicolor;
    document.getElementById("ecor2").style.backgroundColor = locolor;
    document.getElementById("pcor").style.backgroundColor = hicolor;
    document.getElementById("pcor2").style.backgroundColor = locolor;
  } else {
    document.getElementById("ecor").style.backgroundColor = locolor;
    document.getElementById("ecor2").style.backgroundColor = hicolor;
    document.getElementById("pcor").style.backgroundColor = locolor;
    document.getElementById("pcor2").style.backgroundColor = hicolor;
  }

  var ereal = ene[0] * K0 * qq;
  var erecip = ene[1] * K0 * qq;
  var eself = ene[2] * K0 * qq;
  var eBorn = eborn * K0 * qq;
  document.getElementById("ereal").innerHTML  = numformat(ereal);
  document.getElementById("erecip").innerHTML = numformat(erecip);
  document.getElementById("eself").innerHTML  = numformat(eself);
  document.getElementById("eBorn").innerHTML  = numformat(eBorn);
  var preal = ene[0] / (3 * vol) * P0 * qq;
  var precip = ene[1] / (3 * vol) * P0 * qq;
  var pself = ene[2] / (3 * vol) * P0 * qq;
  var pBorn = eborn / vol * P0 * qq;
  document.getElementById("preal").innerHTML  = numformat(preal);
  document.getElementById("precip").innerHTML = numformat(precip);
  document.getElementById("pself").innerHTML  = numformat(pself);
  document.getElementById("pBorn").innerHTML  = numformat(pBorn);

  document.getElementById("etitle").innerHTML = "Energy<br>(" + eunit + ")";
  document.getElementById("ptitle").innerHTML = "Pressure<br>(" + punit + ")";

  var kappa = ene[4], xm = ene[5], km = ene[6];
  document.getElementById("sinfo").innerHTML = "Notes: "
    + "<i>e</i><sup>2</sup>/(4<i>&pi;&epsilon;</i><sub>0</sub>) = " + numformat(K0) + " " + eunit
    + " = " + numformat(P0) + " " + punit
    + (lunit !== "" ? " " + lunit  + "<sup>4</sup>" : "") + ".<br>"
    + "Ewald coefficient: <i>&kappa;</i> = " + numformat(kappa)
    + (lunit !== "" ? " " + lunit + "<sup>&minus;1</sup>" : "") + ","
    + " as in erfc(&minus;<i>&kappa;<sub> </sub>r</i>) / <i>r</i>.<br>"
    + "Number of neighboring cells in each dimension: "
    + xm[0] + " (<i>x</i>), " + xm[1] + " (<i>y</i>), " + xm[2] + " (<i>z</i>).<br>"
    + "Number of wave vectors in each dimension: "
    + km[0] + " (<i>x</i>), " + km[1] + " (<i>y</i>), " + km[2] + " (<i>z</i>).<br>";
}

function init()
{
  viewmat = [[-0.6211501273837585,  0.7805304954451938, -0.07031831149299674],
             [-0.4395286768639746, -0.2726761673313212,  0.8558400843520434],
             [ 0.6488351573900389,  0.5625120918252239,  0.512438371987373]];
  installmouse("animationbox", "scaleinput");
  paint();
}

window.onload = init;

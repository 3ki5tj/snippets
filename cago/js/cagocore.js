/* core functions of the toy nucleic acids */



"use strict";



/* compute the reference bond lengths, angles, dihedrals and pair distances */
function cago_refgeo(go)
{
  var i, j, n = go.n;
  var dx = [0, 0, 0];

  // calculate reference bond lengths, angles, dihedrals
  go.bref = newarr(n - 1); // bonds
  for ( i = 0; i < n - 1; i++ ) {
    go.bref[i] = vdistx(dx, go.xref[i], go.xref[i + 1]);
  }

  go.aref = newarr(n - 2); // angles
  for (i = 1; i < n - 1; i++) {
    go.aref[i - 1]  = vang(go.xref[i - 1], go.xref[i], go.xref[i + 1]);
  }

  go.dref = newarr(n - 3); // dihedrals
  for (i = 0; i < n - 3; i++) {
    go.dref[i] = vdih(go.xref[i], go.xref[i + 1], go.xref[i + 2], go.xref[i + 3]);
  }

  // reference pair distances
  go.r2ref = newarr2d(n, n);
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      go.r2ref[j][i] = go.r2ref[i][j]
          = vsqr( vdiff(dx, go.xref[i], go.xref[j]) );
    }
  }

  return 0;
}



/* read C-alpha coordinates, and amino acid types from a PDB file */
function cago_loadbasic(go, strpdb)
{
  go.xref = [];
  go.iaa = [];
  go.ires = [];

  var sarr = strpdb.split("\n");
  for ( var i = 0; i < sarr.length; i++ ) {
    var s = sarr[i].trim();

    // we only load a single model
    if ( s.substr(0, 3) === "TER"
      || s.substr(0, 3) === "END" ) {
      break;
    }

    // we are only interested in C-alpha atom
    if ( s.substr(0, 6) !== "ATOM  "
      || s.substr(13, 2) !== "CA" ) {
      continue;
    }

    // discard alternative positions
    if ( s.substr(16, 1) !== ' '
      && s.substr(16, 1) !== 'A' ) {
      continue;
    }

    // load the coordinates
    var x = parseFloat( s.substr(30, 8) );
    var y = parseFloat( s.substr(38, 8) );
    var z = parseFloat( s.substr(46, 8) );
    go.xref.push( [x, y, z] );

    var iaa = res2iaa( s.substr(17, 3).trim() );
    go.iaa.push( iaa );

    var ires = parseInt( s.substr(22, 4).trim() );
    go.ires.push( ires );
  }

  go.n = go.xref.length;
  return go.n;
}



/* make the contact map */
function cago_mkcont(go, strpdb, rc, ctype, nsexcl)
{
  var ires, ir, jr, n = go.n;
  var x, xi, dx = [0, 0, 0], rc2 = rc * rc;

  if ( n === 0 ) {
    console.log("no residue");
  }

  x = newarr(n);
  for ( ir = 0; ir < n; ir++ ) {
    x[ir] = [];
  }

  var sarr = strpdb.split("\n");
  for ( var i = 0; i < sarr.length; i++ ) {
    var s = sarr[i].trim();

    // we only load a single model
    if ( s.substr(0, 3) === "TER"
      || s.substr(0, 3) === "END" ) {
      break;
    }

    // we are only interested in C-alpha atom
    if ( s.substr(0, 6) !== "ATOM  " ) {
      continue;
    }

    // discard alternative positions
    if ( s.substr(16, 1) !== ' '
      && s.substr(16, 1) !== 'A' ) {
      continue;
    }

    if ( ctype === "CA"
      && s.substr(13, 2) !== "CA" ) {
      continue;
    }

    if ( ctype === "Heavy"
      && (s.substr(13, 1) === 'H' || s.substr(12, 1) === 'H') ) {
      continue;
    }

    ires = parseInt( s.substr(22, 4) );

    // find the actual residue index
    for ( ir = 0; ir < n; ir++ ) {
      if ( go.ires[ir] == ires ) {
        break;
      }
    }
    if ( ir >= n ) {
      console.log("unknown residue index ", ires);
      continue;
    }

    // load the coordinates
    xi = [parseFloat( s.substr(30, 8) ),
          parseFloat( s.substr(38, 8) ),
          parseFloat( s.substr(46, 8) )];

    x[ir].push( xi );
  }

  go.iscont = newarr2d(n, n);
  for ( ir = 0; ir < n; ir++ ) {
    for ( jr = ir + nsexcl; jr < n; jr++ ) {
      var ia, ja, isc = false;

      // scan over atoms in the two residues
      for ( ia = 0; ia < x[ir].length && !isc; ia++ ) {
        for ( ja = 0; ja < x[jr].length && !isc; ja++ ) {
          if ( vsqr( vdiff(dx, x[ir][ia], x[jr][ja]) ) < rc2 )
            isc = true;
        }
      }

      go.iscont[ir][jr] = go.iscont[jr][ir] = isc;
    }
  }

  x = null;
  return 0;
}




/* return cago_t from pdb file fnpdb
 * `rcc' is the cutoff radius for defining contacts
 * `ctype' is one of PDB_CONTACT_CA, _HEAVY, _ALL
 * `nsexcl' is the number of successive residues to be excluded as contacts
 * e.g., nsexcl = 4 means `a' and `d' in -a-b-c-d- are excluded */
function CaGo(strpdb, kb, ka, kd1, kd3, nbe, nbc, rcc, ctype, nsexcl)
{
  var i, j;

  // read coordinates, and amino acid types
  if ( cago_loadbasic(this, strpdb) < 0 ) {
    return;
  }

  this.dof = this.n*3 - 6;

  cago_mkcont(this, strpdb, rcc, ctype, nsexcl);

  // count the number of contacts
  for ( this.ncont = 0, i = 0; i < this.n - 1; i++ )
    for ( j = i + 1; j < this.n; j++ )
      this.ncont += this.iscont[i][j];

  // compute the reference bond length, angles, etc.
  cago_refgeo(this);

  // initialize the masses
  this.m = newarr(this.n);
  for ( i = 0; i < this.n; i++ ) {
    this.m[i] = 1.0;
  }

  this.kb = kb;
  this.ka = ka;
  this.kd1 = kd1;
  this.kd3 = kd3;
  this.nbe = nbe;
  this.nbc = nbc;
}



/* bond energy 1/2 k (r - r0)^2 */
function potbond(a, b, r0, k, fa, fb)
{
  var dx = [0, 0, 0], r, dr;

  r = vnorm( vdiff(dx, a, b) );
  dr = r - r0;
  if ( fa ) {
    var amp = k * dr / r;
    vsinc(fa, dx, -amp);
    vsinc(fb, dx,  amp);
  }
  return 0.5 * k * dr * dr;
}



/* harmonic angle 1/2 k (ang - ang0)^2 */
function potang(a, b, c, ang0, k, fa, fb, fc)
{
  var dang, amp, ga = [0, 0, 0], gb = [0, 0, 0], gc = [0, 0, 0];

  if ( fa ) { // compute gradient
    dang = vang(a, b, c, ga, gb, gc) - ang0;
    amp = -k * dang;
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
  } else {
    dang = vang(a, b, c) - ang0;
  }
  return 0.5 * k * dang * dang;
}



/* 1-3 dihedral: k1 * (1 - cos(dang)) + k3 * (1 - cos(3*dang)) */
function potdih13(a, b, c, d, ang0, k1, k3, fa, fb, fc, fd)
{
  var dang, amp, ga = [0, 0, 0], gb = [0, 0, 0], gc = [0, 0, 0], gd = [0, 0, 0];

  if ( fa ) {
    dang = vdih(a, b, c, d, ga, gb, gc, gd) - ang0;
    amp  = -k1 * Math.sin(dang) - 3 * k3 * Math.sin(3*dang);
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
    vsinc(fd, gd, amp);
  } else {
    dang = vdih(a, b, c, d) - ang0;
  }
  return k1 * (1 - Math.cos(dang)) + k3 * (1 - Math.cos(3 * dang));
}



/* 12-10 potential: u = 5(rc/r)^12 - 6(rc/r)^10,
 * the minimum is at r = rc, and u = -1 */
function pot1210(a, b, rc2, eps, fa, fb)
{
  var dx = [0, 0, 0], dr2, invr2, invr4, invr6, invr10, amp;

  dr2 = vsqr( vdiff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr4 = invr2 * invr2;
  invr6 = invr4 * invr2;
  invr10 = invr4 * invr6;
  if ( fa ) {
    amp = 60 * eps * (invr2 - 1) * invr10 / dr2;
    vsinc(fa, dx,  amp);
    vsinc(fb, dx, -amp);
  }
  return eps * (5 * invr2 - 6) * invr10;
}



/* repulsive potential: (rc/r)^12 */
function potr12(a, b, rc2, eps, fa, fb)
{
  var dx = [0, 0, 0], dr2, invr2, invr6, u, amp;

  dr2 = vsqr( vdiff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;
  u = eps * invr6 * invr6;
  if ( fa ) {
    amp = 12 * u / dr2;
    vsinc(fa, dx,  amp);
    vsinc(fb, dx, -amp);
  }
  return u;
}



/* force field from C. Clementi, H. Nymeyer, J. N. Onuchic,
 * J. Mol. Biol, Vol. 298 (2000) 937-953 */
CaGo.prototype.force = function(x, f)
{
  var i, j, id, n = this.n;
  var ene = 0, kb = this.kb, ka = this.ka, kd1 = this.kd1, kd3 = this.kd3;
  var nbe = this.nbe, nbc2 = this.nbc * this.nbc;

  if ( f ) {
    for ( i = 0; i < n; i++ ) {
      vzero(f[i]);
    }
  }

  // bonds
  for ( i = 0; i < n - 1; i++ ) {
    ene += potbond(x[i], x[i + 1], go.bref[i], kb, f[i], f[i + 1]);
  }

  // angles
  for ( i = 0; i < n - 2; i++ ) {
    ene += potang(x[i], x[i + 1], x[i + 2], go.aref[i],
              ka, f[i], f[i + 1], f[i + 2]);
  }

  // dihedrals
  for ( i = 0; i < n - 3; i++ ) {
    ene += potdih13(x[i], x[i + 1], x[i + 2], x[i + 3], go.dref[i],
          kd1, kd3, f[i], f[i + 1], f[i + 2], f[i + 3]);
  }

  // non-bonded
  for ( i = 0; i < n - 4; i++ ) {
    for ( j = i + 4; j < n; j++ ) {
      if ( go.iscont[i][j] ) { // contact pair
        ene += pot1210(x[i], x[j], this.r2ref[i][j], nbe, f[i], f[j]);
      } else {
        ene += potr12(x[i], x[j], nbc2, nbe, f[i], f[j]);
      }
    }
  }

  return ene;
};



/* remove center of mass motion, linear and angular */
CaGo.prototype.rmcom = function(x, v)
{
  rmcom(v, this.m, this.n);
  shiftang(x, v, this.m, this.n);
};



/* compute the kinetic energy */
CaGo.prototype.getekin = function(v)
{
  return md_ekin(go.v, go.m, go.n);
};



/* initialize molecular dynamics
 *  o create an initial structure
 *    if `open', start from a nearly-straight chain,
 *      with a disturbance of `rndamp' in the x, y directions
 *    otherwise start from the reference structure,
 *      with a random disturbance of `rndamp'
 *  o initialize the velocity with the center of mass motion removed
 *  o compute the initial force and energy */
CaGo.prototype.initmd = function(open, rndamp, T0)
{
  var i, n = this.n;
  var s, dx = [0, 0, 0];

  this.f  = newarr2d(n, 3);
  this.v  = newarr2d(n, 3);
  this.x  = newarr2d(n, 3);
  this.x1 = newarr2d(n, 3);

  // initialize position
  if ( open ) { // open chain
    for ( i = 0; i < n - 1; i++ ) {
      dx[0] = 1.0;
      dx[1] = rndamp * randgaus();
      dx[2] = rndamp * randgaus();
      s = sqrt( vsqr( dx ) );
      // x_{i+1} = x_i + dx * bref[i]
      vsadd(this.x[i + 1], this.x[i], dx, this.bref[i] / s);
    }
  } else { // copy from xref, slightly disturb it
    for ( i = 0; i < n; i++ ) {
      dx[0] = rndamp * randgaus();
      dx[1] = rndamp * randgaus();
      dx[2] = rndamp * randgaus();
      vadd(this.x[i], this.xref[i], dx);
    }
  }
  rmcom(this.x, this.m, this.n);
  this.epotref = this.force(this.xref, this.f);

  // initialize velocities
  for (i = 0; i < n; i++) {
    s = Math.sqrt( T0 / this.m[i] );
    this.v[i][0] = s * randgaus();
    this.v[i][1] = s * randgaus();
    this.v[i][2] = s * randgaus();
  }
  this.rmcom(this.x, this.v); // remove center of mass motion
  this.ekin = this.getekin(this.v);
  return 0;
};



/* velocity Verlet */
CaGo.prototype.vv = function(fs, dt)
{
  var i, n = this.n;
  var dth = 0.5 * dt * fs;

  for ( i = 0; i < n; i++ ) { // VV part 1
    vsinc(this.v[i], this.f[i], dth / this.m[i]);
    vsinc(this.x[i], this.v[i], dt);
  }

  this.epot = this.force(this.x, this.f); // calculate force

  for ( i = 0; i < n; i++ ) { // VV part 2
    vsinc(this.v[i], this.f[i], dth / this.m[i]);
  }

  return 0;
};



/* Exact velocity rescaling thermostat */
CaGo.prototype.vrescale = function(v, tp, dt)
{
  return md_vrescale(v, this.m, this.n, this.dof, tp, dt);
};



/* compute the RMSD from the reference structure */
CaGo.prototype.rmsd = function(x, xf)
{
  return vrmsd(x, xf, this.xref, this.m, this.n, 0, null, null);
};



/* count the number of native contacts that are formed in the structure `x'
 * this counting process is independent of the process of defining contacts.
 * here, given a set of defined contacts, we simple observe how many pairs
 *   are close enough to be regarded as contacts
 * a contact is formed if the pair distance is <= gam * native-distance
 * return the number of contacts
 * `*Q' is the ratio of formed contacts / the total number of contacts  */
CaGo.prototype.ncontacts = function(x, gam, mat)
{
  var i, j, nct = 0, n = this.n;
  var dx = [0, 0, 0], gam2;

  // determine gam, multiple of the reference distance
  if ( isNaN(gam) ) {
    gam = 1.2;
  }
  gam2 = gam * gam;

  if ( mat ) {
    for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ ) {
        mat[i][j] = 0;
      }
    }
  }

  for ( i = 0; i < n - 1; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      if ( !this.iscont[i][j] ) { // skip a noncontact pair
        continue;
      }
      if ( vsqr( vdiff(dx, x[i], x[j]) ) < this.r2ref[i][j] * gam2 ) {
        if ( mat ) {
          mat[i][j] = mat[j][i] = 1;
        }
        nct++;
      }
    }
  }

  return nct;
}




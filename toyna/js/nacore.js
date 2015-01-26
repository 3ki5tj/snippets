/* core functions of the toy nucleic acids */



"use strict";



function na_initchain2(na)
{
  var i, ib, is, nr = na.nr;
  var rb = 10.0, rs = 4.0, ang = 2*Math.PI/10, dh = 3.4, th, c, s;

  for ( i = 0; i < nr; i++ ) {
    th = ang * i;
    c = Math.cos(th);
    s = Math.sin(th);
    ib = i*2;
    is = i*2 + 1;
    na.x[ib][2] = rb * c;
    na.x[ib][0] = rb * s;
    na.x[ib][1] = dh * i;
    na.x[is][2] = rs * c;
    na.x[is][0] = rs * s;
    na.x[is][1] = dh * i;
  }
}



function na_initchain3(na)
{
  var i, nr = na.nr;

  var ang = 32.7*Math.PI/180, dh = 2.81;
  var rp = 8.710, thp = -70.502*Math.PI/180, dhp = 3.750;
  var rs = 9.064, ths = -43.651*Math.PI/180, dhs = 2.806;
  var rb = 5.6, thb = -43.651*Math.PI/180, dhb = 0.8;

  for ( i = 0; i < nr; i++ ) {
    var th = ang * i;

    var ip = i*3;
    var cp = Math.cos(th + thp);
    var sp = Math.sin(th + thp);
    na.x[ip][2] = rp * cp;
    na.x[ip][0] = rp * sp;
    na.x[ip][1] = dh * i + dhp;

    var is = i*3 + 1;
    var cs = Math.cos(th + ths);
    var ss = Math.sin(th + ths);
    na.x[is][2] = rs * cs;
    na.x[is][0] = rs * ss;
    na.x[is][1] = dh * i + dhs;

    var ib = i*3 + 2;
    if ( na.seq[i] == 'A' ) {
      rb = 5.458;
      thb = -25.976*Math.PI/180;
    } else if ( na.seq[i] == 'C' ) {
      rb = 5.643;
      thb = -33.975*Math.PI/180;
    } else if ( na.seq[i] == 'G' ) {
      rb = 5.631;
      thb = -22.124*Math.PI/180;
    } else if ( na.seq[i] == 'U' ) {
      rb = 5.633;
      thb = -34.102*Math.PI/180;
    }
    var cb = Math.cos(th + thb);
    var sb = Math.sin(th + thb);
    na.x[ib][2] = rb * cb;
    na.x[ib][0] = rb * sb;
    na.x[ib][1] = dh * i + dhb;
  }
}



/* remove the center of mass motion */
function na_rmcom(x, dim, n)
{
  var i;
  var rc = newarr(dim);

  for ( i = 0; i < n; i++ ) {
    vinc(rc, x[i]);
  }
  vsmul(rc, 1.0 / n);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], rc);
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
function na_shiftang(x, v, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) {
    vinc(xc, x[i]);
  }
  vsmul(xc, 1.0/n);
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vinc(am, ang);
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  var inv = rm3_inv(mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  // ang is the solution of M^(-1) * I
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(dv, ang, xi);
    vinc(v[i], dv);
  }
}




function NA(nr, rc)
{
  var i, d, n;

  if ( isNaN(nr) ) {
    this.seq = nr;
    this.nr = nr = this.seq.length;
    console.log( this.seq, this.nr );
  } else {
    this.nr = nr;
    this.seq = "";
    for ( i = 0; i < nr; i++ ) {
      this.seq += "A";
    }
  }

  this.apr = 3; // atoms per residue
  this.n = n = this.nr * this.apr;
  this.dof = n * D - D * (D + 1) / 2;
  this.rc = rc;
  this.l = 10.0; // TODO
  this.x = newarr2d(n, D);
  this.v = newarr2d(n, D);
  this.f = newarr2d(n, D);

  if ( this.apr == 2 ) {
    na_initchain(this);
  } else {
    na_initchain3(this);
  }

  // initialize random velocities
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ ) {
      this.v[i][d] = randgaus();
    }
  }

  na_rmcom(this.v, D, n);
  na_shiftang(this.x, this.v, n);

  this.epot = 0;
  this.ekin = 0;
}



function na_dist2(dx, a, b, l, invl)
{
  return vsqr( vdiff(dx, a, b) );
}



/* compute force, return energy */
NA.prototype.energy_low = function(x)
{
  var dx = newarr(this.dim), dr2, ir6, ep, rc2 = this.rc2;
  var i, j, npr = 0, n = this.n;

  for (ep = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = vsqr( vdiff(dx, x[i], x[j]) );
      if (dr2 < rc2) {
        dr2 = 1 / dr2;
        ir6 = dr2 * dr2 * dr2;
        ep += 4 * ir6 * (ir6 - 1);
        npr++;
      }
    }
  }
  return ep;
};



NA.prototype.energy = function()
{
  this.epot = this.energy_low(this.x);
  return this.epot;
};



/* compute force, return energy */
NA.prototype.force_low = function(x, f)
{
  var dx = newarr(this.dim), fi = newarr(this.dim);
  var dr2, ir6, fs, ep, rc2 = this.rc2;
  var i, j, npr = 0, n = this.n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }
  for (ep = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = vsqr( vdiff(dx, x[i], x[j]) );
      if (dr2 < rc2) {
        dr2 = 1 / dr2;
        ir6 = dr2 * dr2 * dr2;
        fs = ir6 * (48 * ir6 - 24); // f.r
        fs *= dr2; // f.r / r^2
        vsinc(fi, dx, fs);
        vsinc(f[j], dx, -fs);
        ep += 4 * ir6 * (ir6 - 1);
        npr++;
      }
    }
    vinc(f[i], fi);
  }
  return ep;
};



NA.prototype.force = function()
{
  var ret = this.force_low(this.x, this.f);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
};



/* velocity-verlet */
NA.prototype.vv = function(dt)
{
  var i, n = this.n;
  var dth = dt * 0.5, l = this.l;

  for (i = 0; i < n; i++) { // VV part 1
    vsinc(this.v[i], this.f[i], dth);
    vsinc(this.x[i], this.v[i], dt);
  }
  this.force();
  for (i = 0; i < n; i++) { // VV part 2
    vsinc(this.v[i], this.f[i], dth);
  }
};



/* compute the kinetic energy */
function na_ekin(v, n)
{
  var i;
  var ek = 0;
  for ( i = 0; i < n; i++ ) {
    ek += vsqr( v[i] );
  }
  return ek/2;
}



/* exact velocity rescaling thermostat */
function na_vrescale_low(v, n, dof, tp, dt)
{
  var i;
  var c = (dt < 700) ? Math.exp(-dt) : 0;
  var ek1 = na_ekin(v, n);
  var r = randgaus();
  var r2 = randchisqr(dof - 1);
  var ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * Math.sqrt(c * (1 - c) * ek1 * tp / 2);
  ek2 = Math.max(ek2, 0.0);
  var s = Math.sqrt(ek2/ek1);
  for (i = 0; i < n; i++) {
    vsmul(v[i], s);
  }
  return ek2;
}



NA.prototype.vrescale = function(tp, dt)
{
  return na_vrescale_low(this.v, this.n, this.dof, tp, dt);
};



/* displace a random particle i, return i */
NA.prototype.randmv = function(xi, amp)
{
  var i = Math.floor(rand01() * this.n), d;
  for ( d = 0; d < this.dim; d++ ) {
    xi[d] = this.x[i][d] + (rand01() * 2 - 1) * amp;
  }
  return i;
};



/* compute pair energy */
function na_pair(xi, xj, rc2)
{
  var dx = [0,0,0], dr2, invr2, invr6, u;

  dr2 = vsqr( vdiff(dx, xi, xj) );
  if (dr2 < rc2) {
    invr2 = 1 / dr2;
    invr6 = invr2 * invr2 * invr2;
    u  = 4 * invr6 * (invr6 - 1);
    return [true, u];
  } else {
    return [false, 0.0];
  }
}



/* return the energy change from displacing x[i] to xi */
NA.prototype.depot = function(i, xi)
{
  var j, n = this.n;
  var l = this.l, invl = 1/l, rc2 = this.rc2, u, ret;

  u = 0;
  for ( j = 0; j < n; j++ ) { // pair
    if ( j === i ) {
      continue;
    }
    ret = na_pair(this.x[i], this.x[j], rc2);
    if ( ret[0] ) {
      u -= ret[1];
    }
    ret = na_pair(xi, this.x[j], rc2);
    if ( ret[0] ) {
      u += ret[1];
    }
  }
  return u;
};



/* commit a particle displacement */
NA.prototype.commit = function(i, xi, du)
{
  vcopy(this.x[i], xi);
  this.ep0 += du;
  this.epot += du;
};



/* Metropolis algorithm */
NA.prototype.metro = function(amp, bet)
{
  var acc = 0;
  var xi = newarr(this.dim);

  var i = this.randmv(xi, amp);
  var ret = this.depot(i, xi);
  var du = ret[0];
  if ( du < 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -bet * du ) );
  }
  if ( acc ) {
    this.commit(i, xi, du);
    return 1;
  }
  return 0;
};




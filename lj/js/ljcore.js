/* core functions of the LJ (Lennard-Jones) object */



"use strict";



function lj_initfcc(lj)
{
  if ( lj.dim == 2 ) return lj_initfcc2d(lj);
  else if ( lj.dim == 3 ) return lj_initfcc3d(lj);
}



function lj_gettail(lj, rho, n)
{
  if ( lj.dim == 2 ) return lj_gettail2d(lj, rho, n);
  else if ( lj.dim == 3 ) return lj_gettail3d(lj, rho, n);
}



function lj_setrho(lj, rho)
{
  lj.rho = rho;
  lj.vol = lj.n / rho;
  lj.l = Math.pow(lj.vol, 1.0 / lj.dim);
  lj.rc = Math.min( lj.l / 2, lj.rcdef );
  lj.rc2 = lj.rc * lj.rc;
  var irc = 1 / lj.rc;
  var irc2 = irc * irc;
  var irc6 = irc2 * irc2 * irc2;
  lj.epot_shift = 4 * irc6 * (irc6 - 1);
  var ret = lj_gettail(lj, rho, lj.n); // to be defined in lj2d.h or lj3d.h
  lj.epot_tail = ret[0];
  lj.p_tail = ret[1];
}



/* remove the center of mass motion */
function lj_rmcom(x, dim, n)
{
  var i;
  var rc = newarr(dim);

  for ( i = 0; i < n; i++ )
    vinc(rc, x[i]);
  vsmul(rc, 1./n);
  for ( i = 0; i < n; i++ )
    vdec(x[i], rc);
}



function lj_shiftang(x, v, n)
{
  if ( D === 2 ) lj_shiftang2d(x, v, n);
  else if ( D === 3 ) lj_shiftang3d(x, v, n);
}



function LJ(n, dim, rho, rcdef)
{
  var i, d;

  this.n = n;
  this.dim = dim;
  this.dof = n * dim - dim * (dim + 1) / 2;
  this.rcdef = rcdef;
  this.x = newarr2d(n, dim);
  this.v = newarr2d(n, dim);
  this.f = newarr2d(n, dim);

  lj_setrho(this, rho);
  lj_initfcc(this); /* to be defined later in lj2d.js or lj3d.js */

  // initialize random velocities
  for ( i = 0; i < n; i++ )
    for ( d = 0; d < dim; d++ )
      this.v[i][d] = randgaus();

  lj_rmcom(this.v, dim, n);
  lj_shiftang(this.x, this.v, n); /* to be defined later in lj2d.js or lj3d.js */

  this.epot = 0;
  this.eps = 0;
  this.ep0 = 0;
  this.vir = 0;
  this.ekin = 0;
}



/* OO version of lj_setrho() */
LJ.prototype.setrho = function(rho)
{
  lj_setrho(this, rho);
}



function lj_pbcdist2(dx, a, b, l, invl)
{
  return vsqr( vpbc(vdiff(dx, a, b), l, invl) );
}



/* compute force and virial, return energy */
LJ.prototype.energy_low = function(x)
{
  var dx = newarr(this.dim), dr2, dr6, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1 / l;
  var i, j, npr = 0, n = this.n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if (dr2 > rc2) continue;
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      vir += dr6 * (48 * dr6 - 24); // f.r
      ep += 4 * dr6 * (dr6 - 1);
      npr++;
    }
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
}



LJ.prototype.energy = function()
{
  ret = this.energy_low(this.x);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
}



/* compute force and virial, return energy */
LJ.prototype.force_low = function(x, f)
{
  var dx = newarr(this.dim), fi = newarr(this.dim);
  var dr2, dr6, fs, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1/l;
  var i, j, npr = 0, n = this.n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if (dr2 > rc2) continue;
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48 * dr6 - 24); // f.r
      vir += fs; // f.r
      fs *= dr2; // f.r / r^2
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += 4 * dr6 * (dr6 - 1);
      npr++;
    }
    vinc(f[i], fi);
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
}



LJ.prototype.force = function()
{
  var ret = this.force_low(this.x, this.f);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot
}



/* compute pressure */
LJ.prototype.calcp = function(tp)
{
  return (this.dof * tp + this.vir) / (this.dim * this.vol) + this.p_tail;
}



/* velocity-verlet */
LJ.prototype.vv = function(dt)
{
  var i, n = this.n;
  var dth = dt * .5, l = this.l;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(this.v[i], this.f[i], dth);
    vsinc(this.x[i], this.v[i], dt);
    vwrap(this.x[i], l);
  }
  this.force();
  for (i = 0; i < n; i++) /* VV part 2 */
    vsinc(this.v[i], this.f[i], dth);
}



/* compute the kinetic energy */
function lj_ekin(v, n)
{
  var i;
  var ek = 0;
  for ( i = 0; i < n; i++ ) ek += vsqr( v[i] );
  return ek/2;
}



/* exact velocity rescaling thermostat */
function lj_vrescale_low(v, n, dof, tp, dt)
{
  var i;
  var c = (dt < 700) ? Math.exp(-dt) : 0;
  var ek1 = lj_ekin(v, n);
  var r = randgaus();
  var r2 = randchisqr(dof - 1);
  var ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * Math.sqrt(c * (1 - c) * ek1 * tp / 2);
  if (ek2 < 0) ek2 = 0;
  var s = Math.sqrt(ek2/ek1);
  for (i = 0; i < n; i++) vsmul(v[i], s);
  return ek2;
}



LJ.prototype.vrescale = function(tp, dt)
{
  return lj_vrescale_low(this.v, this.n, this.dof, tp, dt);
}



/* position Langevin barostat, with coordinates only
 * NOTE: the first parameter is the degree of freedom
 * the scaling is r = r*s
 * set cutoff to half of the box */
LJ.prototype.langp0 = function(dt, tp, pext, ensx)
{
  var pint, amp, s, dlnv;
  var i;

  pint = this.calcp(tp);
  amp = Math.sqrt(2 * dt);
  dlnv = ((pint - pext) * this.vol / tp + 1 - ensx) * dt + amp * randgaus();
  s = Math.exp( dlnv / D );
  this.vol *= Math.exp( dlnv );
  this.setrho(this.n / this.vol);
  for ( i = 0; i < this.n; i++ )
    vsmul(this.x[i], s);
  this.force();
}



/* displace a random particle i, return i */
LJ.prototype.randmv = function(xi, amp)
{
  var i = Math.floor(rand01() * this.n), d;
  for ( d = 0; d < D; d++ )
    xi[d] = this.x[i][d] + (rand01() * 2 - 1) * amp;
  return i;
}



/* compute pair energy */
function lj_pair(xi, xj, l, invl, rc2)
{
  var dx = [0,0,0], dr2, invr2, invr6;

  dr2 = lj_pbcdist2(dx, xi, xj, l, invl);
  if (dr2 < rc2) {
    invr2 = 1 / dr2;
    invr6 = invr2 * invr2 * invr2;
    vir = invr6 * (48 * invr6 - 24); /* f.r */
    u  = 4 * invr6 * (invr6 - 1);
    return [true, u, vir];
  } else {
    return [false, 0.0, 0.0];
  }
}



/* return the energy change from displacing x[i] to xi */
LJ.prototype.depot = function(i, xi)
{
  var j, n = this.n;
  var l = this.l, invl = 1/l, rc2 = this.rc2, u, vir, ret;

  u = 0;
  vir = 0.0;
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) continue;
    ret = lj_pair(this.x[i], this.x[j], l, invl, rc2);
    if ( ret[0] ) {
      u -= ret[1];
      vir -= ret[2];
    }
    ret = lj_pair(xi, this.x[j], l, invl, rc2);
    if ( ret[0] ) {
      u += ret[1];
      vir += ret[2];
    }
  }
  return [u, vir];
}



/* commit a particle displacement */
LJ.prototype.commit = function(i, xi, du, dvir)
{
  vcopy(this.x[i], xi);
  this.ep0 += du;
  this.epot += du;
  this.vir += dvir;
}



/* Metropolis algorithm */
LJ.prototype.metro = function(amp, bet)
{
  var acc = 0;
  var xi = newarr(this.dim);

  var i = this.randmv(xi, amp);
  var ret = this.depot(i, xi);
  var du = ret[0];
  var dvir = ret[1];
  if ( du < 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -bet * du ) );
  }
  if ( acc ) {
    this.commit(i, xi, du, dvir);
    return 1;
  }
  return 0;
}





/* core functions of the LJ (Lennard-Jones) object */



"use strict";



/* initialize a fcc lattice */
function lj_initfcc2d(lj)
{
  var i, j, id = 0, n = lj.n;
  var n1 = Math.floor( Math.sqrt(2*n) + 0.999999 ); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;

  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      if ( (i+j) % 2 === 0 ) {
        // add some noise to prevent two atoms happened to
        // be separated precisely by the cutoff distance,
        // which might be half of the box
        lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
        lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
        id++;
      }
    }
  }
}



/* initialize a fcc lattice */
function lj_initfcc3d(lj)
{
  var i, j, k, id, n = lj.n;

  var n1 = Math.floor(Math.pow(2*n, 1.0/3) + 0.999999); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++) {
    for (j = 0; j < n1 && id < n; j++) {
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 === 0) {
          /* add some noise to prevent two atoms happened to
           * be separated by precisely some special cutoff distance,
           * which might be half of the box */
          lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][2] = (k + 0.5) * a + noise * (2*rand01() - 1);
          id++;
        }
      }
    }
  }
}



function lj_initfcc(lj)
{
  if ( lj.dim == 2 ) {
    return lj_initfcc2d(lj);
  } else if ( lj.dim == 3 ) {
    return lj_initfcc3d(lj);
  }
}



/* get the tail correction */
function lj_gettail2d(lj, rho, n)
{
  var irc, irc2, irc6, utail, ptail;

  irc = 1 / lj.rc;
  irc2 = irc * irc;
  irc6 = irc2 * irc2 * irc2;
  utail = Math.PI * rho * n * (0.4*irc6 - 1) * irc2 * irc2;
  ptail = Math.PI * rho * rho * (2.4*irc6 - 3) * irc2 * irc2;
  return [utail, ptail];
}



/* get the tail correction */
function lj_gettail3d(lj, rho, n)
{
  var irc, irc3, irc6, utail, ptail;

  irc = 1 / lj.rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = 8 * Math.PI * rho * n / 9 * (irc6 - 3) * irc3;
  ptail = 32 * Math.PI * rho * rho / 9 * (irc6 - 1.5) * irc3;
  return [utail, ptail];
}



function lj_gettail(lj, rho, n)
{
  if ( lj.dim == 2 ) {
    return lj_gettail2d(lj, rho, n);
  } else if ( lj.dim == 3 ) {
    return lj_gettail3d(lj, rho, n);
  }
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



function LJ(n, dim, rho, rcdef)
{
  var i, d;

  this.n = n;
  this.dim = dim;
  this.dof = n * dim - dim; // - dim * (dim + 1) / 2;
  this.rcdef = rcdef;
  this.x = newarr2d(n, dim);
  this.v = newarr2d(n, dim);
  this.f = newarr2d(n, dim);

  lj_setrho(this, rho);
  lj_initfcc(this);

  // initialize random velocities
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < dim; d++ ) {
      this.v[i][d] = randgaus();
    }
  }

  rmcom(this.v, null, n);
  shiftang(this.x, this.v, null, n);

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
};



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
      if (dr2 < rc2) {
        dr2 = 1 / dr2;
        dr6 = dr2 * dr2 * dr2;
        vir += dr6 * (48 * dr6 - 24); // f.r
        ep += 4 * dr6 * (dr6 - 1);
        npr++;
      }
    }
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
};



LJ.prototype.energy = function()
{
  ret = this.energy_low(this.x);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
};



/* compute force and virial, return energy */
LJ.prototype.force_low = function(x, f)
{
  var dx = newarr(this.dim), fi = newarr(this.dim);
  var dr2, dr6, fs, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1/l;
  var i, j, npr = 0, n = this.n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if (dr2 < rc2) {
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
    }
    vinc(f[i], fi);
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
};



LJ.prototype.force = function()
{
  var ret = this.force_low(this.x, this.f);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
};



/* compute pressure */
LJ.prototype.calcp = function(tp)
{
  return (this.dof * tp + this.vir) / (this.dim * this.vol) + this.p_tail;
};



/* velocity-Verlet with force scaling */
LJ.prototype.vv_fs = function(dt, fs)
{
  var i, n = this.n;
  var dth = dt * 0.5 * fs, l = this.l;

  for (i = 0; i < n; i++) { // VV part 1
    vsinc(this.v[i], this.f[i], dth);
    vsinc(this.x[i], this.v[i], dt);
    vwrap(this.x[i], l);
  }
  this.force();
  for (i = 0; i < n; i++) { // VV part 2
    vsinc(this.v[i], this.f[i], dth);
  }
};




/* velocity-Verlet */
LJ.prototype.vv = function(dt)
{
  this.vv_fs(dt, 1);
};



/* compute the kinetic energy */
function lj_ekin(v, n)
{
  return md_ekin(v, null, n);
}



/* exact velocity rescaling thermostat */
LJ.prototype.vrescale = function(tp, dt)
{
  return md_vrescale(this.v, null, this.n, this.dof, tp, dt);
};



/* adaptive velocity rescaling for an asymptotic microcanonical ensemble */
LJ.prototype.adaptvrescale = function(tp, alpha)
{
  return md_adaptvrescale(this.v, null, this.n, this.dof, tp, alpha);
};



/* Nose-Hoover rescaling thermostat */
LJ.prototype.nhchain = function(tp, dt, zeta, zmass)
{
  return md_nhchain(this.v, null, this.n, this.dof, tp, dt,
      zeta, zmass);
};



/* exact velocity rescaling thermostat */
LJ.prototype.langevin = function(tp, dt)
{
  return md_langevin(this.v, null, this.n, tp, dt);
};



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
  s = Math.exp( dlnv / this.dim );
  this.vol *= Math.exp( dlnv );
  this.setrho(this.n / this.vol);
  for ( i = 0; i < this.n; i++ ) {
    vsmul(this.x[i], s);
  }
  this.force();
};



/* displace a random particle i, return i */
LJ.prototype.randmv = function(xi, amp)
{
  var i = Math.floor(rand01() * this.n), d;
  for ( d = 0; d < this.dim; d++ ) {
    xi[d] = this.x[i][d] + (rand01() * 2 - 1) * amp;
  }
  return i;
};



/* compute pair energy */
function lj_pair(xi, xj, l, invl, rc2)
{
  var dx = [0,0,0], dr2, invr2, invr6, vir, u;

  dr2 = lj_pbcdist2(dx, xi, xj, l, invl);
  if (dr2 < rc2) {
    invr2 = 1 / dr2;
    invr6 = invr2 * invr2 * invr2;
    vir = invr6 * (48 * invr6 - 24); // f.r
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
  for ( j = 0; j < n; j++ ) { // pair
    if ( j === i ) {
      continue;
    }
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
};



/* commit a particle displacement */
LJ.prototype.commit = function(i, xi, du, dvir, l)
{
  vwrap( vcopy(this.x[i], xi), l );
  this.ep0 += du;
  this.epot += du;
  this.vir += dvir;
};



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
    this.commit(i, xi, du, dvir, this.l);
    return 1;
  }
  return 0;
};




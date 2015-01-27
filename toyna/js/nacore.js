/* core functions of the toy nucleic acids */



"use strict";



var TIS_IP = 0;
var TIS_IS = 1;
var TIS_IB = 2;

// bond length parameters
var R0_PS = 4.660;
var K0_PS = 64.0;
var R0_SP = 3.766;
var K0_SP = 23.0;
var R0_SB = [4.700, 4.157, 4.811, 4.163];
var K0_SB = 10.0;

// bond angle parameters
var A0_PSB = [D2R( 94.06), D2R( 86.52), D2R( 98.76), D2R( 86.31)];
var A0_BSP = [D2R(106.69), D2R(108.27), D2R(106.17), D2R(108.28)];
var A0_PSP = D2R(83.53);
var KA_PSB = 5.0;
var KA_BSP = 5.0;
var KA_PSP = 20.0;

// WCA parameters
var WCA_SIG = 3.2;
var WCA_SIG2 = WCA_SIG * WCA_SIG;
var WCA_EPS = 1.0;



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
    na.m[ib] = 100.0; // TODO
    na.x[ib][2] = rb * c;
    na.x[ib][0] = rb * s;
    na.x[ib][1] = dh * i;
    na.m[is] = 100.0; // TODO
    na.x[is][2] = rs * c;
    na.x[is][0] = rs * s;
    na.x[is][1] = dh * i;
  }
}



function na_initchain3(na)
{
  var i, ic, nr = na.nr;

  var ang = 32.7*Math.PI/180, dh = 2.81;
  var rp = 8.710, thp = D2R(-70.502), dhp = 3.750;
  var rs = 9.214, ths = D2R(-41.097), dhs = 2.864;
  var rbarr  = [5.458, 5.643, 5.631, 5.633];
  var thbarr = [D2R(-25.976), D2R(-33.975), D2R(-22.124), D2R(-34.102)];
  var dhbarr = [0.742, 0.934, 0.704, 0.932];
  var rb, thb, dhb, th;
  var bmass = [134.0, 110.0, 150.0, 111.0];

  for ( i = 0; i < nr; i++ ) {
    var th = ang * i;

    na.m[i*3 + TIS_IP] = 95.0;
    na.x[i*3 + TIS_IP][2] = rp * Math.cos(th + thp);
    na.x[i*3 + TIS_IP][0] = rp * Math.sin(th + thp);
    na.x[i*3 + TIS_IP][1] = dh * i + dhp;

    na.m[i*3 + TIS_IS] = 99.0;
    na.x[i*3 + TIS_IS][2] = rs * Math.cos(th + ths);
    na.x[i*3 + TIS_IS][0] = rs * Math.sin(th + ths);
    na.x[i*3 + TIS_IS][1] = dh * i + dhs;

    ic = na.iseq[ i ];
    rb = rbarr[ ic ];
    thb = thbarr[ ic ];
    dhb = dhbarr[ ic ];
    na.m[i*3 + TIS_IB] = bmass[ ic ];
    na.x[i*3 + TIS_IB][2] = rb * Math.cos(th + thb);
    na.x[i*3 + TIS_IB][0] = rb * Math.sin(th + thb);
    na.x[i*3 + TIS_IB][1] = dh * i + dhb;
  }
}



/* remove the center of mass motion */
function na_rmcom(x, m, n)
{
  var i;
  var xc = [0,0,0];
  var mtot = 0;

  for ( i = 0; i < n; i++ ) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c X r
 *   */
function na_shiftang(x, v, m, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  var mtot = 0;

  for (i = 0; i < n; i++) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);

  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vsinc(am, ang, m[i]);

    xx += m[i] * xi[0] * xi[0];
    yy += m[i] * xi[1] * xi[1];
    zz += m[i] * xi[2] * xi[2];
    xy += m[i] * xi[0] * xi[1];
    yz += m[i] * xi[1] * xi[2];
    zx += m[i] * xi[2] * xi[0];
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
  // ang is the solution of M^(-1) * L
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(dv, ang, xi);
    vinc(v[i], dv);
  }
}




function NA(nr, tp, debyel)
{
  var i, d, n;
  var ch2int = {
    "A": 0,
    "C": 1,
    "G": 2,
    "U": 3,
    "T": 3
  };

  if ( isNaN(nr) ) {
    this.seq = nr.toUpperCase();
    this.nr = nr = this.seq.length;
    this.iseq = newarr(nr);

    for ( i = 0; i < nr; i++ ) {
      this.iseq[i] = ch2int[ this.seq[i] ];
    }
    console.log( this.seq, this.nr );
  } else {
    this.nr = nr;
    this.seq = "";
    this.iseq = newarr(nr);
    for ( i = 0; i < nr; i++ ) {
      this.seq += "A";
      this.iseq[i] = 0;
    }
  }

  this.apr = 3; // atoms per residue
  this.n = n = this.nr * this.apr;
  this.dof = n * D - D * (D + 1) / 2;
  this.tp = tp;
  this.debyel = debyel;
  this.m = newarr(n);
  this.x = newarr2d(n, D);
  this.v = newarr2d(n, D);
  this.f = newarr2d(n, D);

  if ( this.apr == 2 ) {
    na_initchain2(this);
  } else {
    na_initchain3(this);
  }

  // initialize random velocities
  for ( i = 0; i < n; i++ ) {
    var amp = Math.sqrt( BOLTZK * tp / this.m[i] );
    for ( d = 0; d < D; d++ ) {
      this.v[i][d] = amp * randgaus();
    }
  }

  na_rmcom(this.v, this.m, n);
  na_shiftang(this.x, this.v, this.m, n);

  this.epot = 0;
  this.ekin = 0;
}



function na_dist2(dx, a, b, l, invl)
{
  return vsqr( vdiff(dx, a, b) );
}



/* compute force, return energy
 * for the TIS model */
NA.prototype.energyTIS_low = function(x, tp, debyel)
{
  var i, j, ic, n = this.n, nr = this.nr;

  var ep = 0;

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = this.iseq[i];
    ep += ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS]);
    ep += ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondlen(R0_SP, K0_SP,
                     x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP]);
    }
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = this.iseq[i];
    ep += ebondang(A0_PSB[ic], KA_PSB,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondang(A0_BSP[ic], KA_BSP,
                     x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP]);
      ep += ebondang(A0_PSP, KA_PSP,
                     x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP]);
    }
  }

  // exclusion volume
  for ( i = 0; i < n - 1; i++ ) {
    for (j = i + 1; j < n; j++) {
      ep += ewca(WCA_SIG2, WCA_EPS, x[i], x[j]);
    }
  }

  // electrostatic interaction
  var eps = getdielecwater(tp);
  var Q = getchargeQ(tp);
  var QQ = Q*Q*KE2/eps;
  for ( i = 0; i < nr - 1; i++ ) {
    for ( j = i + 1; j < nr; j++ ) {
      ep += echargeDH(QQ, debyel,
                      x[i*3 + TIS_IP], x[j*3 + TIS_IP]);
    }
  }
  return ep;
};



NA.prototype.energy = function()
{
  this.epot = this.energyTIS_low(this.x, this.tp, this.debyel);
  return this.epot;
};



/* compute force, return energy */
NA.prototype.forceTIS_low = function(x, f)
{
  var i, j, ic, n = this.n, nr = this.nr;

  var ep = 0;
  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }

  // bond length
  for ( i = 0; i < nr; i++ ) {
    ic = this.iseq[i];
    ep += ebondlen(R0_PS, K0_PS,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS],
                   f[i*3 + TIS_IP], f[i*3 + TIS_IS]);
    ep += ebondlen(R0_SB[ic], K0_SB,
                   x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondlen(R0_SP, K0_SP,
                     x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    }
  }

  // bond angle
  for ( i = 0; i < nr; i++ ) {
    ic = this.iseq[i];
    ep += ebondang(A0_PSB[ic], KA_PSB,
                   x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[i*3 + TIS_IB],
                   f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[i*3 + TIS_IB]);
    if ( i < nr - 1 ) {
      ep += ebondang(A0_BSP[ic], KA_BSP,
                     x[i*3 + TIS_IB], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IB], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
      ep += ebondang(A0_PSP, KA_PSP,
                     x[i*3 + TIS_IP], x[i*3 + TIS_IS], x[(i + 1)*3 + TIS_IP],
                     f[i*3 + TIS_IP], f[i*3 + TIS_IS], f[(i + 1)*3 + TIS_IP]);
    }
  }

  // excluded volume
  for ( i = 0; i < n - 1; i++ ) {
    for (j = i + 1; j < n; j++) {
      ep += ewca(WCA_SIG2, WCA_EPS, x[i], x[j], f[i], f[j]);
    }
  }

  // electrostatic interaction
  var eps = getdielecwater(tp);
  var Q = getchargeQ(tp);
  var QQ = Q*Q*KE2/eps;
  for ( i = 0; i < nr - 1; i++ ) {
    for ( j = i + 1; j < nr; j++ ) {
      ep += echargeDH(QQ, debyel,
                      x[i*3 + TIS_IP], x[j*3 + TIS_IP],
                      f[i*3 + TIS_IP], f[j*3 + TIS_IP]);
    }
  }
  return ep;
};



NA.prototype.force = function()
{
  this.epot = this.forceTIS_low(this.x, this.f, this.tp, this.debyel);
  return this.epot;
};



/* velocity-verlet */
NA.prototype.vv = function(dt)
{
  var i, n = this.n;
  var dth = dt * 0.5;

  for (i = 0; i < n; i++) { // VV part 1
    vsinc(this.v[i], this.f[i], dth / this.m[i]);
    vsinc(this.x[i], this.v[i], dt);
  }
  this.force();
  for (i = 0; i < n; i++) { // VV part 2
    vsinc(this.v[i], this.f[i], dth / this.m[i]);
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
  tp *= BOLTZK;
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
  var u, ret;

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




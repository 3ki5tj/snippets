/* molecular dynamics routines */



"use strict";



/* remove the center of mass motion */
function rmcom(x, m, n)
{
  var i;
  var xc = newarr(D);
  var wt, mtot = 0.0;

  for ( i = 0; i < n; i++ ) {
    wt = m ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



/* annihilate the total angular momentum */
function shiftang2d(x, v, m, n)
{
  var i;
  var am, r2, xc = [0, 0], xi = [0, 0];
  var wt, mtot = 0.0;

  for (i = 0; i < n; i++) {
    wt = m ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

  for (am = r2 = 0.0, i = 0; i < n; i++) {
    wt = m ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    am += wt * vcross2d(xi, v[i]);
    r2 += wt * vsqr(x[i]);
  }

  am = -am / r2;
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am*xi[1];
    v[i][1] +=  am*xi[0];
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
function shiftang3d(x, v, m, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  var wt, mtot = 0.0;

  for (i = 0; i < n; i++) {
    wt = m ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

  for (i = 0; i < n; i++) {
    wt = m ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vsinc(am, ang, wt);

    xx += wt * xi[0] * xi[0];
    yy += wt * xi[1] * xi[1];
    zz += wt * xi[2] * xi[2];
    xy += wt * xi[0] * xi[1];
    yz += wt * xi[1] * xi[2];
    zx += wt * xi[2] * xi[0];
  }
  mat[0][0] = yy + zz;
  mat[1][1] = xx + zz;
  mat[2][2] = xx + yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  var inv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  minv(inv, mat);

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



function shiftang(x, v, m, n)
{
  if ( D === 2 ) {
    shiftang2d(x, v, m, n);
  } else if ( D === 3 ) {
    shiftang3d(x, v, m, n);
  }
}



/* compute the kinetic energy */
function md_ekin(v, m, n)
{
  var ek = 0;

  for ( var i = 0; i < n; i++ ) {
    var wt = m ? m[i] : 1.0;
    ek += wt * vsqr( v[i] );
  }
  return ek * 0.5;
}


/* exact velocity rescaling thermostat */
function md_vrescale(v, m, n, dof, tp, dt)
{
  var c = Math.exp(-dt);
  var ek1 = md_ekin(v, m, n);
  var r = randgaus();
  var r2 = randchisqr(dof - 1);
  var ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
          + 2 * r * Math.sqrt(c * (1 - c) * ek1 * tp / 2);
  if ( ek2 < 0 ) {
    ek2 = 0;
  }
  var s = Math.sqrt(ek2 / ek1);
  for ( var i = 0; i < n; i++ ) {
    vsmul(v[i], s);
  }
  return ek2;
}



/* adaptive velocity rescaling for an asymptotic microcanonical ensemble */
function md_adaptvrescale(v, m, n, dof, tp, alpha, sm1max)
{
  var ek = md_ekin(v, m, n);
  var ekref = dof * tp * 0.5;
  var sm1 = (ekref - ek)/ek * alpha;
  if ( !sm1max ) sm1max = 0.5;
  if ( sm1 > sm1max ) sm1 = sm1max;
  else if ( sm1 < -sm1max ) sm1 = -sm1max;
  var s = Math.sqrt( 1 + sm1 );
  //console.log(alpha, ek, ekref, s);
  for ( var i = 0; i < n; i++ ) {
    vsmul(v[i], s);
  }
  ek *= s * s;
  return ek;
}




/* Nose-Hoover chain thermostat */
function md_nhchain(v, m, n, dof, tp, dt, zeta, zmass)
{
  var i, j, nnhc = zeta.length, s, GQ;

  var mvv = md_ekin(v, m, n) * 2;

  for ( j = nnhc - 1; j >= 0; j-- ) {
    s = ( j == nnhc - 1 ) ? 1 : Math.exp(-zeta[j+1]*dt*0.25);
    GQ = ( j == 0 ) ? (mvv - dof * tp) : (zmass[j-1] * zeta[j-1] * zeta[j-1] - tp);
    zeta[j] = (zeta[j] * s + GQ /zmass[j] * dt*0.5) * s;
  }

  s = Math.exp( -zeta[0] * dt );
  for ( i = 0; i < n; i++ ) {
    vsmul(v[i], s);
  }
  mvv *= s * s;

  for ( j = 0; j < nnhc; j++ ) {
    s = ( j == nnhc - 1 ) ? 1 : Math.exp(-zeta[j+1]*dt*0.25);
    GQ = ( j == 0 ) ? (mvv - dof * tp) : (zmass[j-1] * zeta[j-1] * zeta[j-1] - tp);
    zeta[j] = (zeta[j] * s + GQ /zmass[j] * dt*0.5) * s;
  }

  return mvv * 0.5;
}



/* Langevin-dynamics for the velocities */
function md_langevin(v, m, n, tp, dt)
{
  var s = Math.exp(-dt);
  var ek = 0;
  for ( var i = 0; i < n; i++ ) {
    var mass = ( m ? m[i] : 1.0 );
    var r = Math.sqrt(tp/mass * (1 - s * s));
    for ( var j = 0; j < D; j++ )
      v[i][j] = v[i][j] * s + randgaus() * r;
    ek += 0.5 * mass * vsqr(v[i]);
  }
  return ek;
}



/* randomly swap the velocities of k pairs of particles */
function md_vscramble(v, m, n, k)
{
  var i, j, l, d;
  var vi, vj, smi, smj;

  for ( l = 0; l < k; l++ ) {
    i = Math.floor(rand01() * n);
    j = (i + 1 + Math.floor(rand01() * (n - 1))) % n;

    if ( m ) {
      smi = Math.sqrt( m[i] );
      smj = Math.sqrt( m[j] );
    } else {
      smi = smj = 1;
    }

    for ( d = 0; d < D; d++ ) {
      vi = smi * v[i][d];
      vj = smj * v[j][d];
      v[i][d] = vj / smi;
      v[j][d] = vi / smj;
    }
  }
  return md_ekin(v, m, n);
}





/* two-dimensional Potts model */



"use strict";



/* initialize an lxl Potts model */
function Potts(l, q)
{
  var i, n;

  this.l = l;
  this.n = n = l * l;
  this.q = q;
  this.s = newarr(n);
  for ( i = 0; i < n; i++ ) {
    this.s[i] = 0;
  }
  this.E = -2*n;
  this.proba = newarr(10);
  this.queue = newarr(n);
  this.used = newarr(n);
}



/* set transition probability */
Potts.prototype.setproba = function(bet)
{
  var x = Math.exp(-bet), y = x;
  this.proba[0] = 1;
  this.proba[1] = y; y *= x;
  this.proba[2] = y; y *= x;
  this.proba[3] = y; y *= x;
  this.proba[4] = y;
};



/* pick a random site, count neighbors with different spins
 * return the site */
Potts.prototype.pick = function()
{
  var id, ix, iy, l, n, ixp, ixm, iyp, iym, s, so;

  l = this.l;
  n = this.n;
  id = Math.floor( rand01() * n );
  ix = id % l;
  iy = id - ix;
  ixp = ( ix + 1 ) % l;
  ixm = ( ix + l - 1 ) % l;
  iyp = ( iy + l ) % n;
  iym = ( iy + n - l ) % n;
  so = this.s[id];
  this.sn = (so + 1 + Math.floor( (this.q - 1) * rand01() )) % this.q;
  this.h = 0;
  s = this.s[iy  + ixp]; this.h += (s == so) - (s == this.sn);
  s = this.s[iy  + ixm]; this.h += (s == so) - (s == this.sn);
  s = this.s[iyp + ix ]; this.h += (s == so) - (s == this.sn);
  s = this.s[iym + ix ]; this.h += (s == so) - (s == this.sn);
  return id;
}



/* flip site id, with h being the energy before the flip */
Potts.prototype.flip = function(id, sn, h)
{
  this.s[id] = sn;
  this.E += h;
  return this.E;
}


/* compute total energy and magnetization */
Potts.prototype.energy = function()
{
  var l, n, i, j, e;

  e = 0;
  l = this.l;
  n = l * l;
  for ( i = 0; i < n; i += l ) {
    for ( j = 0; j < l; j++ ) {
      var id = i + j;
      var idr = i + (j + 1) % l;
      var idu = (i + l) % n + j;
      var s = this.s[id];
      var su = this.s[idu];
      var sr = this.s[idr];
      e += (s == su) + (s == sr);
    }
  }
  return this.E = -e;
}



/* add spin j to the queue if s[j] is different from s
 * return the spin */
Potts.prototype.addtoqueue = function(j, so, sn, r)
{
  var sj = this.s[j];

  if ( sj == so && !this.used[j] && rand01() < r ) {
    this.queue[ this.cnt ] = j;
    this.cnt++;
    this.used[j] = 1;
  }
  return (sj == so) - (sj == sn);
}



/* Wolff algorithm,
 * padd should be precomputed as 1 - exp(-beta) */
Potts.prototype.wolff = function(padd)
{
  var l = this.l, n = this.n, i, ix, iy, h = 0;

  // randomly selected a seed
  var id = Math.floor ( rand01() * n );
  var so = this.s[id];
  var sn = (so + 1 + Math.floor(rand01() * (this.q - 1))) % this.q;
  this.cnt = 0;
  this.queue[ this.cnt++ ] = id;
  for ( i = 0; i < n; i++ ) {
    this.used[i] = 0;
  }
  this.used[id] = 1;

  // go through spins in the queue
  for ( i = 0; i < this.cnt; i++ ) {
    id = this.queue[i];
    this.s[id] = sn;
    // add neighbors of i with the same spins
    ix = id % l;
    iy = id - ix;
    h += this.addtoqueue(iy + (ix + 1) % l,     so, sn, padd);
    h += this.addtoqueue(iy + (ix + l - 1) % l, so, sn, padd);
    h += this.addtoqueue((iy + l) % n + ix,     so, sn, padd);
    h += this.addtoqueue((iy + n - l) % n + ix, so, sn, padd);
  }

  this.E += h;
  return 0;
}




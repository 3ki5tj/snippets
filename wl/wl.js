


"use strict";



var WL_FLATNESSABS = 0x0010; // use (hmax-hmin)/(hmax+hmin) for flatness



function WL(n0, n1, lnf0, flatness, frac, c, flags)
{
  this.n0 = n0;
  this.n1 = n1;
  this.n = n1 - n0;
  this.h = newarr(this.n);
  this.v = newarr(this.n);
  this.lnf = lnf0;
  this.tot = 0;
  this.isinvt = 0;
  this.flatness = flatness;
  this.frac = frac;
  this.c = c;
  this.flags = flags;
}



/* compute the total of the histogram */
function wl_gethtot(h, n)
{
  var htot = 0, i;

  for ( i = 0; i < n; i++ ) {
    htot += h[i];
  }
  return htot;
}



/* clear the histogram */
function wl_clearh(h, n)
{
  for ( var i = 0; i < n; i++ ) {
    h[i] = 0.0;
  }
}



/* trim the bottom of the potential */
function wl_trimv(v, n)
{
  var i, vmin = v[0];

  for ( i = 1; i < n; i++ ) {
    if ( v[i] < vmin ) {
      vmin = v[i];
    }
  }
  for ( i = 0; i < n; i++ ) {
    v[i] -= vmin;
  }
}



WL.prototype.trimv = function()
{
  wl_trimv(this.v, this.n);
}



/* add an entry, update the histogram and potential */
WL.prototype.add = function(i)
{
  i -= this.n0;
  if ( i < 0 || i >= this.n ) {
    console.log("wl: out of range", i, this.n);
    return -1;
  }
  this.h[i] += 1.0;
  this.v[i] += this.lnf;
  this.tot += 1.0;
  return 0;
};



/* compute the histogram flatness from the absolute value */
function wl_getflatnessabs(h, n)
{
  var hmin = hmax = h[0];
  for ( var i = 1; i < n; i++ ) {
    if ( h[i] > hmax ) {
      hmax = h[i];
    } else if ( h[i] < hmin ) {
      hmin = h[i];
    }
  }
  return hmax > hmin ? (hmax - hmin) / (hmax + hmin) : 1.0;
}



/* compute the histogram flatness from the standard deviation */
function wl_getflatnessstd(h, n)
{
  var y, sh = 0, shh = 0;

  for ( var i = 0; i < n; i++ ) {
    y = h[i];
    sh += y;
    shh += y * y;
  }
  if ( sh <= 0 ) return 1.0;
  sh /= n;
  shh = shh / n - sh * sh;
  if ( shh < 0 ) shh = 0;
  return Math.sqrt(shh) / sh;
}



WL.prototype.getflatness = function()
{
  if ( this.flags & WL_FLATNESSABS ) {
    return wl_getflatnessabs(this.h, this.n);
  } else {
    return wl_getflatnessstd(this.h, this.n);
  }
};



/* lnf = 1/t */
WL.prototype.lnfinvt = function()
{
  return this.c * this.n / this.tot;
};



/* update lnf, return 1 if the Wang-Landau stage is switched */
WL.prototype.updatelnf = function()
{
  var flatness, nlnf, lnfinvt;

  if ( this.isinvt ) {
    this.lnf = this.lnfinvt();
    return 0;
  }

  flatness = this.getflatness();
  if ( flatness < this.flatness ) {
    nlnf = this.lnf * this.frac;
    lnfinvt = this.lnfinvt();
    console.log(this.c, this.n, this.tot, lnfinvt, nlnf);
    if ( nlnf < lnfinvt ) {
      console.log("changing lnf from " + this.lnf.toExponential(3) + " to " + lnfinvt.toExponential(3) + " (1/t), flatness " + roundto(flatness*100.0, 2) + "%");
      this.isinvt = 1;
      this.lnf = lnfinvt;
    } else {
      console.log("changing lnf from " + this.lnf.toExponential(3) + " to " + nlnf.toExponential(3) + " (1/t " + lnfinvt.toExponential(3) + "), flatness " + roundto(flatness*100.0, 2) + "%");
      this.lnf = nlnf;
      wl_clearh(this.h, this.n);
    }
    return 1;
  }
  return 0;
}





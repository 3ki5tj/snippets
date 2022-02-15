// radial distribution function

"use strict";

function Rdf(rmax, dr, dim, nstDeposit, sel)
{
  this.rmax = rmax;
  this.dr = dr;
  this.sel = sel;
  this.dim = dim;
  this.nstDeposit = nstDeposit;
  this.reset();
}
  
Rdf.prototype.reset = function() {
  this.n = Math.floor(this.rmax/this.dr);
  //console.log(this.n, this.rmax, this.dr);
  this.hist = new Array(this.n + 1);
  for (let i = 0; i < this.n; i++) {
    this.hist[i] = 0;
  }
  this.plot = null;
  this.gr = new Array(this.n);
  this.count = 0;
  this.vsphr = new Array(this.n);
  for (let i = 0; i < this.n; i++) {
    let r = this.dr * i;
    if (this.dim === 3) {
      this.vsphr[i] = 4 * Math.PI/3 * (Math.pow(r + this.dr, 3) - Math.pow(r, 3));
    } else {
      this.vsphr[i] = Math.PI * (Math.pow(r + this.dr, 2) - Math.pow(r, 2));
    }
  }
  this.step = 0;
};

Rdf.prototype.add = function(lj) {
  if (++this.step % this.nstDeposit !== 0) {
    return;
  }

  var dx = newarr(lj.dim), dr2, r, rc = lj.l*0.5;
  var l = lj.l, invl = 1 / l, x = lj.x;
  var n = lj.n;
  var dr = this.dr, hist = this.hist, ir;

  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      r = Math.sqrt(dr2);
      if (r > rc) {
        continue;
      }
      ir = Math.floor(r/dr);
      hist[ir] += 1;
    }
  }

  this.count += 1;
};

Rdf.prototype.get = function(lj) {
  let n = this.n;
  let fac = this.count * 0.5 * lj.n * (lj.n - 1) / lj.vol;
  for (let i = 0; i < n; i++) {
    this.gr[i] = this.hist[i] / this.vsphr[i] / fac;
  }
  return this.gr;
};  
  
// update the RDF plot
Rdf.prototype.updatePlot = function(lj)
{
  let data = "r,g(r),1\n";
  let gr = this.get(lj);
  for (let i = 0; i < this.n; i++ ) {
    let r = (i + 0.5)*this.dr;
    data += "" + roundto(r, 6) + "," + gr[i] + ",1\n";
  }

  if (this.plot === null) {
    let options = LJUI.getDefaultDygraphOption();
    Object.assign(options, {
      xlabel: '<small>r</small>',
      ylabel: '<small>g(r)</small>',
      includeZero: true,
      drawPoints: true,
      xRangePad: 2,
    });    
    this.plot = new Dygraph(document.querySelector(this.sel), data, options);
  } else {
    this.plot.updateOptions({ file: data });
  }
}

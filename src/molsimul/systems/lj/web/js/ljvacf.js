"use strict";


function VacfVecSeq(n)
{
  this.n = n;
  this.reset();
}


VacfVecSeq.prototype.reset = function()
{
  this.seq = new Array();
  this.ave = new Array(this.n + 1);
  for (let i = 0; i <= this.n; i++) {
    this.ave[i] = new Averager();
  }
}

VacfVecSeq.prototype.add = function(v)
{
  let n = this.n;
  let seq = this.seq;
  let len = seq.length;
  let ave = this.ave;

  // update the averages
  this.ave[0].add( vdot(v, v) );
  for (let i = 0; i < len; i++) {
    ave[len-i].add( vdot(seq[i], v) );
  }

  // update the sequence
  if (len >= n) {
    seq.shift();
  }
  // add the new vector
  seq.push(v);
}


VacfVecSeq.prototype.getAcf = function()
{
  let n = this.n;
  let x = new Array(n + 1);
  for (let i = 0; i <= n; i++) {
    x[i] = this.ave[i].getMean();
  }
  return x;
}


// sel: DOM selector
function Vacf(np, dt, tmax, sel)
{
  this.np = np;
  this.dt = dt;
  this.tmax = tmax;
  this.itmax = Math.floor(tmax/dt);
  this.sel = sel;
  this.reset();
}


Vacf.prototype.reset = function() {
  this.vecSeq = new Array(this.np);
  for (let j = 0; j < this.np; j++) {
    this.vecSeq[j] = new VacfVecSeq(this.itmax);
  }
  this.plot = null;
};


Vacf.prototype.add = function(lj) {
  let np = this.np;
  let v = lj.v;
  let dim = lj.dim;

  for (let j = 0; j < np; j++ ) {
    let vj;
    if (dim === 2) {
      vj = [v[j][0], v[j][1]];
    } else {
      vj = [v[j][0], v[j][1], v[j][2]];
    }
    this.vecSeq[j].add(vj);
  }
};


Vacf.prototype.autoTruncate = function(vacf, threshold, multiple)
{
  let vacf0 = vacf[0];
  let n = this.itmax + 1;

  for (let i = 1; i < n; i++) {
    if (vacf[i] < vacf0*threshold) {
      let icutoff = Math.min(i * multiple, n);
      return vacf.slice(0, icutoff);
    }
  }
  return vacf;
}


Vacf.prototype.get = function()
{
  let np = this.np;
  let itmax = this.itmax;
  let seqs = new Array(np);

  for (let j = 0; j < np; j++) {
    seqs[j] = this.vecSeq[j].getAcf();
  }

  let vacf = new Array(itmax + 1);
  for (let i = 0; i <= itmax; i++) {
    // for any time interval, do the particle average
    let vv = 0;
    for (let j = 0; j < np; j++) {
      vv += seqs[j][i];
    }
    vv /= np;

    vacf[i] = vv;
  }

  // truncate the tail of the ACF
  return this.autoTruncate(vacf, 0.2, 8);
};

Vacf.prototype.updatePlot = function()
{
  let data = "Time,Velocity,0\n";
  let vacf = this.get();

  let vacf0 = vacf[0];
  for (let i = 0; i < vacf.length; i++ ) {
    data += "" + roundto(i*this.dt, 6) + ","
               + roundto(vacf[i]/vacf0, 4) + ",0\n";
  }

  //console.log(vlen, icutoff);

  if ( this.plot === null ) {
    let options = LJUI.getDefaultDygraphOption();
    Object.assign(options, {
      xlabel: '<small>Time</small>',
      ylabel: '<small>Velocity auto-correlation function</small>',
      includeZero: true,
      xRangePad: 2,
    });
    this.plot = new Dygraph(document.querySelector(this.sel), data, options);
  } else {
    this.plot.updateOptions({ file: data });
  }
}

"use strict";

function Averager()
{
  this.count = 0.0;
  this.sum = 0.0;
  this.sqr = 0.0;
}

Averager.prototype.add = function(x) {
  this.count += 1;
  this.sum += x;
  this.sqr += x * x;
};

Averager.prototype.getMean = function(x) {
  if (this.count > 0) {
    return this.sum / this.count;
  } else {
    return 0.0;
  }
};


Averager.prototype.getVar = function(x) {
  if (this.count > 0) {
    var xx = this.sqr / this.count;
    var mean = this.getMean();
    return xx  - mean * mean;
  } else {
    return 0.0;
  }
};

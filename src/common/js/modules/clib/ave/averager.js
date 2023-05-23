"use strict";

function Averager()
{
  this.count = 0.0;
  this.sum = 0.0;
  this.sum2 = 0.0;
}

Averager.prototype.add = function(x)
{
  this.count += 1;
  this.sum += x;
  this.sum2 += x * x;
};


Averager.prototype.getMean = function()
{
  if (this.count > 0) {
    return this.sum / this.count;
  } else {
    return 0.0;
  }
};


Averager.prototype.getVar = function()
{
  if (this.count > 0) {
    var mean = this.getMean();
    var x2 = this.sum2 / this.count;
    return x2  - mean * mean;
  } else {
    return 0.0;
  }
};


Averager.prototype.getStd = function()
{
  return Math.sqrt(this.getVar());
};


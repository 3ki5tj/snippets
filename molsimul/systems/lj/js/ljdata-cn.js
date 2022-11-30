"use strict"

function LJData(dim, n, dof, rho, tp, isMC, isMicrocanonical)
{
  this.dim = dim;
  this.n = n;
  this.dof = dof;
  this.rho = rho;
  this.tp = tp;
  this.isMC = isMC;
  this.isMicrocanonical = isMicrocanonical;

  // reference values from the equation of state
  if ( this.dim === 3 ) {
    let ref = lj_eos3dPVEhBH(rho, tp);
    this.Uref = ref[0];
    this.Pref = ref[1];
  } else {
    this.Uref = null;
    this.Pref = null;
  }

  this.reset();
}

LJData.prototype.reset = function()
{
  this.count = 0;
  this.aveU = new Averager();
  this.aveP = new Averager();
  if (this.isMC) {
    this.aveAcc = new Averager();
  } else {
    this.aveK = new Averager();
    if (this.isMicrocanonical) {
      this.aveBeta = new Averager();
      this.aveDbde = new Averager();  
    }  
  }
};

LJData.prototype.add = function(U, P, ekin, acc)
{
  this.aveU.add(U);
  this.aveP.add(P);
  if (this.isMC) {
    this.aveAcc.add(acc);
  } else {
    this.aveK.add(2*ekin/this.dof);
    if (this.isMicrocanonical) {
      let bet = (this.dof * 0.5 - 1) / ekin;
      let dbde = -bet / ekin;
      this.aveBeta.add(bet);
      this.aveDbde.add(dbde);
    }  
  }
  this.count++;
};

LJData.prototype.getTemperature = function()
{
  if ( this.isMicrocanonical ) {
    let aveBeta = this.aveBeta.getMean();
    return 1.0/Math.max(aveBeta, 1e-5);
  } else {
    return this.aveK.getMean();
  }
};

// estimate beta'(E) for the microcanonical ensemble
LJData.prototype.getDbde = function()
{
  if ( this.count < 3 ) {
    return -1.0 / (this.tp * this.tp * (0.5 * this.dof - 1));
  }
  
  let betaVar = this.aveBeta.getVar();
  let dbde1 = this.aveDbde.getMean();
  let dbde = dbde1 + betaVar;
  return Math.min(dbde, dbde1*0.05);
}

LJData.prototype.getReportUP = function()
{
  let s;

  let aveU = this.aveU.getMean();
  s = '势能/粒子数：' + roundto(aveU, 3);
  if ( this.Uref ) {
    s += ' （<span class="tag">参考值:</span> ' + roundto(this.Uref, 3) + "）";
  }
  s += "<br>";

  let aveP = this.aveP.getMean();
  s += '压强: ' + roundto(aveP, 3);
  if ( this.Pref ) {
    s += ' （<span class="tag">参考值:</span> ' + roundto(this.Pref, 3) + "）";
  }

  s += "<br>收集到 " + this.count + " 个样本";

  return s;
};

LJData.prototype.getReport = function()
{
  let sinfo;

  if (this.isMC) {
    let aveAcc = this.aveAcc.getMean();
    sinfo = "蒙特卡罗接受率：" + roundto(100.0 * aveAcc, 2) + "%<br>"
          + this.getReportUP();
  } else {
    let tpMeasured = this.getTemperature();

    sinfo = '实测温度：'
            + roundto(tpMeasured, 3)
            + ' （<span class="tag">设定值：</span> ' + roundto(this.tp, 3) + '）<br>';

    if ( this.isMicrocanonical ) {
        sinfo += '<span class="math"><i>N d&beta;/dE</i>：</span>'
                 + roundto(this.n * this.getDbde(this.tp, this.dof), 4) + '<br>';
    }
    sinfo += this.getReportUP();
  }
  return sinfo;
};
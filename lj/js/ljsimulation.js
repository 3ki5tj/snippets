"use strict";



function LJSimulation()
{
  this.readSystemParams();
  this.readDataParams();
}


LJSimulation.prototype.readSystemParams = function()
{
  let n = getInt("#n", 216);
  let dim = getInt("#dimension", 2);
  if ( dim === 2 || dim === 3 ) {
    // set the global variable for dimension
    D = dim;
  }
  let rho = getFloat("#density", 0.8);

  let rcdef = getFloat("#rcutoff", 1000.0);

  let simulMethod = qSel("#simul-method").value;
  this.simulMethod = simulMethod;
  this.isMC = (simulMethod === "mc");

  let lj = new LJ(n, D, rho, rcdef);
  this.lj = lj;

  lj.force();

  qSel("#rcutoff-actual").innerHTML = roundto(lj.rc, 4);
};


LJSimulation.prototype.readDataParams = function()
{
  let lj = this.lj;

  let tp = getFloat("#temperature", 1.5);
  this.tp = tp;

  let simulMethod = this.simulMethod;

  let mddt = getFloat("#md-dt", 0.002);
  if (mddt < 0.0001) {
    mddt = 0.0001;
    setTimeout(function() {
      qSel("#md-dt").value = mddt;
    }, 5000);
  };
  this.mddt = mddt;

  this.mcamp = getFloat("#mc-amp", 0.2);

  let thtype = qSel("#thermostat-type").value;
  let vrDamp = getFloat("#vr-damp", 10.0);
  let nhclen = getInt("#nhc-len", 10);
  let nhcmass1 = getFloat("#nhc-mass1", 10.0);
  let nhcmass2 = getFloat("#nhc-mass2", 1.0);
  let langDamp = getFloat("#langevin-damp", 1.0);

  this.thstat = new LJThermostat(thtype, this.tp, this.mddt*0.5,
                                 vrDamp,
                                 nhclen, nhcmass1, nhcmass2,
                                 langDamp);

  this.dof = ( thtype === "langevin" ) ? lj.dim * lj.n : lj.dof;

  this.isMicrocanonical = (thtype === "adaptive-v-rescaling");

  this.statData = new LJData(lj.dim, lj.n, this.dof, lj.rho, tp, this.isMC, this.isMicrocanonical);

  // initialize the RDF plot
  let rdfMcNstDeposit = getInt("#rdf-mc-nst-deposit", 2000);
  let rdfMdNstDeposit = getInt("#rdf-md-nst-deposit", 20);
  let rdfNstDeposit = (simulMethod === "md") ? rdfMdNstDeposit : rdfMcNstDeposit;
  this.rdf = new Rdf(lj.l*0.5, 0.02, lj.dim, rdfNstDeposit, "#rdf-plot");

  // initialize the VACF plot
  if (simulMethod === "md") {

    let vacfN = getInt("#vacf-np", 1);
    if (vacfN >= lj.n) {
      vacfN = lj.n;
    } else if (vacfN < 1) {
      vacfN = 1;
    }

    let vacfTmax = getFloat("#vacf-tmax", 10.0);
    this.vacf = new Vacf(vacfN, this.mddt, vacfTmax, "#vacf-plot");

  } else {

    this.vacf = null;
    qSel("#vacf-plot").innerHTML = "";

  }  
}


LJSimulation.prototype.doMC = function(nsteps)
{
  let lj = this.lj;

  for ( let istep = 0; istep < nsteps; istep++ ) {
    let acc = lj.metro(this.mcamp, 1.0 / this.tp);
    this.statData.add(lj.epot/lj.n, lj.calcp(this.tp), null, acc);
    this.rdf.add(lj);
  }
};

LJSimulation.prototype.doMD = function(nsteps)
{
  let lj = this.lj;
  let thstat = this.thstat;

  for (let istep = 0; istep < nsteps; istep++) {
    thstat.apply(lj, this.statData);
    lj.vv(this.mddt);
    let ekin = thstat.apply(lj, this.statData);

    this.statData.add(lj.epot/lj.n, lj.calcp(this.tp), ekin);
    this.rdf.add(lj);
    this.vacf.add(lj);
  }
};

LJSimulation.prototype.pulse = function(nsteps)
{
  if ( this.simulMethod === "md" ) {
    this.doMD(nsteps);
  } else if ( this.simulMethod === "mc" ) {
    this.doMC(nsteps);
  }

  qSel("#msg-info").innerHTML = this.statData.getReport();
};

LJSimulation.prototype.resetData = function()
{
  if (this.statData) {
    this.statData.reset();
  }

  if (this.rdf) {
    this.rdf.reset();
  }

  if (this.vacf) {
    this.vacf.reset();
  }
};
"use strict";

function LJThermostat(type, tp, dt,
                      vrDamp,
                      nhclen, nhcmass1, nhcmass2,
                      langDamp)
{
  this.type = type;
  this.tp = tp;
  this.dt = dt;
  
  this.vrDamp = vrDamp;

  this.zeta = new Array(nhclen);
  this.zmass = new Array(nhclen);
  for (let i = 0; i < nhclen; i++ ) {
    this.zeta[i] = 0.0;
    this.zmass[i] = ( i === 0 ) ? nhcmass1 : nhcmass2;
  }

  this.langDamp = langDamp;
}

LJThermostat.prototype.apply = function(lj, statData)
{
  let type = this.type;
  let tp = this.tp;
  let dt = this.dt;

  if ( type === "v-rescaling" ) {
    return lj.vrescale(tp, dt * this.vrDamp);
  } else if ( type === "nh" ) {
    return lj.nhchain(tp, dt, this.zeta, this.zmass);
  } else if ( type === "langevin" ) {
    return lj.langevin(tp, dt * this.langDamp);
  } else if ( type == "adaptive-v-rescaling" ) {
    let thalpha = 0.5 / (statData.count + 1); // 1/t factor, 0.5 for half step
    let dbde = statData.getDbde();
    return lj.adaptvrescale(tp, thalpha, dbde);
  }
}

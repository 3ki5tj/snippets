"use strict";

function LJSimulationController()
{
  this.simul = null;

  // paint settings
  this.timer = null;
  this.timerInterval = 100; // in milliseconds
  this.paintCount = 0;

  let self = this;

  // initialize the controller
  LJUI.init();

  this.mouse3dController = new Mouse3DController("#ljbox", "#ljscale",
    function() {
      self.paint();
    },
    function() {
      self.pauseOrResume();
    });

  document.querySelector("#start").addEventListener("click", function(){
    self.start();
  });
  document.querySelector("#pause").addEventListener("click", function(){
    self.pause();
  });
  document.querySelector("#resume").addEventListener("click", function(){
    self.resume();
  });
  document.querySelector("#stop").addEventListener("click", function(){
    self.stop();
  });
  document.querySelector("#reset-data").addEventListener("click", function(){
    self.resetData();
  });

  // show Pause instead of Resume initially
  self.setPauseOrResume(true);

  // respond to changes of critical parameters: restart simulation
  let onCriticalParamChange = function() {
    if ( self.timer !== null ) {
      self.stop();
      self.start();
    }
  };
  document.querySelectorAll(".c-param").forEach(function(item) {
    item.addEventListener("change", function(){
      onCriticalParamChange();
    });
    item.addEventListener("input", function(){
      onCriticalParamChange();
    });
  });

  // respond to changes of data-critical parameters:
  // reset data, but keep the simulation running
  let onDataParamChange = function() {
    self.resetData();
    let simul = self.simul;
    if (simul) {
      let isRunning = (self.timer !== null);
      if (isRunning) {
        self.pause();
      }
      simul.readDataParams();
      if (isRunning) {
        self.resume();
      }
    }
  };
  document.querySelectorAll(".d-param").forEach(function(item) {
    item.addEventListener("change", function(){
      onDataParamChange();
    });
    item.addEventListener("input", function(){
      onDataParamChange();
    });
  });

}

LJSimulationController.prototype.paint = function()
{
  let simul = this.simul;

  if (!simul) {
    return;
  }

  let lj = simul.lj;
  if (this.mouse3dController) {
    let canvasScale = this.mouse3dController.getScale();
    let viewMatrix = this.mouse3dController.getViewMatrix();
    let canvas = document.querySelector("#ljbox");
  
    if ( lj.dim === 2 ) {
      ljdraw2d(lj, canvas, canvasScale);
    } else if ( lj.dim === 3 ) {
      ljdraw3d(lj, canvas, canvasScale, viewMatrix);
    }
  }

  if (++this.paintCount % 10 === 0) {
    simul.rdf.updatePlot(lj);
    if (simul.vacf !== null) {
      simul.vacf.updatePlot();
    }
  }
};


LJSimulationController.prototype.setPauseOrResume = function(showPause)
{
  if (showPause) {
    qSel("#resume").style.display = "none";
    qSel("#pause").style.display = "";
  } else {
    qSel("#resume").style.display = "";
    qSel("#pause").style.display = "none";
  }
};


LJSimulationController.prototype.startTimer = function()
{
  let self = this;

  this.timer = setInterval(function(){

    let simul = self.simul;
    if (simul) {
      let nstepsPerPulse;
      if (simul.isMC) {
        let mcNstepsPerSecond = getInt("#mc-nsteps-ps", 10000);
        nstepsPerPulse = mcNstepsPerSecond * self.timerInterval / 1000;
      } else {
        let mdNstepsPerSecond = getInt("#md-nsteps-ps", 1000);
        nstepsPerPulse = mdNstepsPerSecond * self.timerInterval / 1000;
      }
      simul.pulse(nstepsPerPulse);
      self.paint();
    }

  }, this.timerInterval);

  this.setPauseOrResume(true);
};

LJSimulationController.prototype.stopTimer = function()
{
  if (this.timer) {
    clearInterval(this.timer);
    this.timer = null;  
  }
  this.setPauseOrResume(false);
}


LJSimulationController.prototype.resetData = function()
{
  if (!this.simul) {
    return;
  }

  this.simul.resetData();
  this.paintCount = 0;
};


LJSimulationController.prototype.stop = function()
{
  this.stopTimer();
  this.resetData();
  this.mouse3dController.resetViewMatrix();
};


LJSimulationController.prototype.pause = function()
{
  if (!this.simul) {
    return;
  }

  this.stopTimer();
};



LJSimulationController.prototype.resume = function()
{
  if (!this.simul) {
    return;
  }

  if (this.timer === null) {
    this.startTimer();
  }
};


LJSimulationController.prototype.start = function()
{
  if (this.timer) {
    this.stop();
  }
  this.simul = new LJSimulation();
  this.startTimer();
  if (this.mouse3dController) {
    this.mouse3dController.reset();
  }
};


LJSimulationController.prototype.pauseOrResume = function()
{
  if (!this.simul) {
    this.start();
  } else {
    if (this.timer) {
      this.pause();
    } else {
      this.resume();
    }
  }
};

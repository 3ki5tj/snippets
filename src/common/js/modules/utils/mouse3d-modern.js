/* mouse control for rotating and scaling a 3D object
 * Note: onclick is not captured */

"use strict";


// in seconds, shorter than this, trigger a click
let _clickDuration = 0.1;


function MouseRotation3D(target, paintFunc, clickFunc) {
  this.dragStarted = false;
  this.x = -1;
  this.y = -1;
  this.viewMatrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

  this.paintFunc = paintFunc;
  this.clickFunc = clickFunc;

  let self = this;
  target.addEventListener("mousedown", function(e){
    self.mousedown(e);
  });
  target.addEventListener("mousemove", function(e){
    self.mousemove(e);
  });
  target.addEventListener("mouseup", function(e){
    self.mouseup(e);
  });
}

MouseRotation3D.prototype.mousedown = function(e) {
  e = e || window.event;
  this.x = parseInt(e.clientX);
  this.y = parseInt(e.clientY);
  this.dragStarted = true;
  this.dragStartTime = new Date().getTime();
  e.preventDefault();
};

MouseRotation3D.prototype.mousemove = function(e) {
  if (!this.dragStarted) {
    return;
  }
  e = e || window.event;
  if (this.x >= 0 && this.y >= 0) {
    let target = e.target ? e.target : e.srcElement;
    this.viewMatrix = mxrot3d(this.viewMatrix, 180.0*(e.clientY - this.y)/target.height);
    this.viewMatrix = myrot3d(this.viewMatrix, 180.0*(e.clientX - this.x)/target.width);
    if (this.paintFunc) {
      this.paintFunc();
    }
  }

  this.x = e.clientX;
  this.y = e.clientY;
  //console.log("touchmove", this.x, this.y, touchObj);
  e.preventDefault();
}

MouseRotation3D.prototype.mouseup = function(e) {
  if (!this.dragStarted) {
    return;
  }
  this.dragStarted = false;
  let duration = (new Date().getTime() - this.dragStartTime)*0.001;
  //console.log(duration, this.clickFunc);
  if (duration < _clickDuration && this.clickFunc) {
    this.clickFunc();
  }
  this.x = -1;
  this.y = -1;
  e.preventDefault();
}



function TouchRotation3D(target, paintFunc, clickFunc) {
  this.dragStarted = false;
  this.x = -1;
  this.y = -1;
  this.viewMatrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

  this.scale = 1.0;
  this.scaler = null;
  this.paintFunc = paintFunc;
  this.clickFunc = clickFunc;

  let self = this;
  target.addEventListener("touchstart", function(e){
    self.touchstart(e);
  });
  target.addEventListener("touchmove", function(e){
    self.touchmove(e);
  });
  target.addEventListener("touchend", function(e){
    self.touchend(e);
  });
}


TouchRotation3D.prototype.touchstart = function(e) {
  // http://www.javascriptkit.com/javatutors/touchevents.shtml
  if (e.changedTouches.length === 2) {
    return;
  }
  let touchObj = e.changedTouches[0]; // reference first touch point (ie: first finger)
  this.x = parseInt(touchObj.clientX); // get x position of touch point relative to left edge of browser
  this.y = parseInt(touchObj.clientY);
  this.dragStarted = true;
  this.dragStartTime = new Date().getTime();
  //console.log("touchstart", this.x, this.y, touchObj);
  e.preventDefault();
};

TouchRotation3D.prototype.touchmove = function(e) {
  if (e.changedTouches.length === 2) {
    return;
  }
  if (!this.dragStarted) {
    return;
  }
  let touchObj = e.changedTouches[0];
  if (this.x >= 0 && this.y >= 0) {
    let target = e.target ? e.target : e.srcElement;
    this.viewMatrix = mxrot3d(this.viewMatrix, 180.0*(touchObj.clientY - this.y)/target.height);
    this.viewMatrix = myrot3d(this.viewMatrix, 180.0*(touchObj.clientX - this.x)/target.width);
    if (this.paintFunc) {
      this.paintFunc();
    }
  } 
  this.x = touchObj.clientX;
  this.y = touchObj.clientY;
  //console.log("touchmove", this.x, this.y, touchObj);
  e.preventDefault();
}

TouchRotation3D.prototype.touchend = function(e) {
  if (e.changedTouches.length === 2) {
    return;
  }
  if (!this.dragStarted) {
    return;
  }
  let duration = (new Date().getTime() - this.dragStartTime)*0.001;
  //console.log(duration);
  if (duration < _clickDuration && this.clickFunc) {
    this.clickFunc();
  }
  this.dragStarted = false;
  this.x = -1;
  this.y = -1;
  e.preventDefault();
}


CanvasScalar.prototype.applyScale = function(scale)
{
  let scalar = this.scalar;
  this.scale = scale;
  if (scalar) {
    scalar.value = this.scale;
    if (scale > scalar.max) {
      scalar.value = this.scale = scalar.max;
    } else if (scale < scalar.min) {
      scalar.value = this.scale = scalar.min;
    }
  }

  // repaint if necessary
  if (this.paintFunc) {
    this.paintFunc();
  }
}


CanvasScalar.prototype.installPinch = function()
{
  let self = this;
  let target = this.target;

  target.addEventListener("touchstart", function(e) {
    if (e.touches.length === 2) {
      self.pinchStarted = true;
      self.pinchDist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY);
    }
  });

  target.addEventListener("touchmove", function(e) {
    if (!self.pinchStarted) {
      return;
    }
    if (e.touches.length === 2) {
      let newPinchDist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY);
      let scale = newPinchDist / self.pinchDist;
      scale = Math.min(1.1, Math.max(0.9, scale))
      //console.log(self.pinchDist, newPinchDist, scale);
      self.applyScale(scale * self.scale);
      self.pinchDist = newPinchDist;
    }
  });

  target.addEventListener("touchend", function(e) {
    if (!self.pinchStarted) {
      return;
    }
    // for some reason we can't control the touches.length on touch end so precisely
    //if (e.touches.length === 2) {
      self.pinchStarted = false;
      self.pinchDist = 100;
    //}
  });

}

// install the mouse wheel event
CanvasScalar.prototype.installWheel = function()
{
  let target = this.target;
  let self = this;
  let wheelHandler = function(e) {
    let delta = 0; // positive for scrolling up
    e = e || window.event;
    if ( e.wheelDelta ) { // IE/Opera
      delta = e.wheelDelta / 120;
    } else if ( e.detail ) { // Firefox
      delta = -e.detail / 3;
    }

    let scale = self.scale;
    if ( delta > 0 ) {
      scale *= 1.05;
    } else if ( delta < 0 ) {
      scale *= 0.95;
    }
    self.applyScale(scale);

    //console.log("wheel", delta);
    e.preventDefault();
    e.returnValue = false;
  };

  if ( target.addEventListener ) {
    // for IE9+, Chrome, Safari, Opera
    target.addEventListener('mousewheel', wheelHandler, false);
    // for Firefox
    target.addEventListener('DOMMouseScroll', wheelHandler, false);
  } else { // for IE 6/7/8
    target.attachEvent("onmousewheel", wheelHandler);
  }
}



function CanvasScalar(target, scalar)
{
  this.target = target;
  this.scalar = scalar;
  this.scale = scalar.value;
  this.installWheel();

  this.installPinch();

  let self = this;
  scalar.addEventListener("change", function(e) {
    self.scale = self.scalar.value;
    if (this.paintFunc) {
      paintFunc();
    }
  });
}


function Mouse3DController(boxSel, scalarSel, paintFunc, clickFunc)
{
  let target = document.querySelector(boxSel);
  let scalar = document.querySelector(scalarSel);

  this.scaleController = new CanvasScalar(target, scalar, paintFunc);

  if ("ontouchstart" in document.documentElement) {
    this.rotationController = new TouchRotation3D(target, paintFunc, clickFunc);
  } else {
    this.rotationController = new MouseRotation3D(target, paintFunc, clickFunc);
  }
}

Mouse3DController.prototype.getViewMatrix = function() {
  return this.rotationController.viewMatrix;
}

Mouse3DController.prototype.resetViewMatrix = function() {
  munit(this.rotationController.viewMatrix);
}

Mouse3DController.prototype.getScale = function() {
  return this.scaleController.scale;
}

Mouse3DController.prototype.reset = function() {
  this.resetViewMatrix();
}

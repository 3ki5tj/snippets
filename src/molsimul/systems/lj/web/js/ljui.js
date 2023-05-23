"use strict";

let LJUI = (function() {
  function qSel(sel)
  {
    return document.querySelector(sel);
  }

  function qSelAll(sel)
  {
    return document.querySelectorAll(sel);
  }

  // show/hide MC/MD parameters
  function onChangeSimulMethod()
  {
    // hide all thermostat parameters first
    qSelAll(".simul-method-param").forEach(function(item){
      item.style.display = "none";
    });

    // enable the specific simulation method parameters
    let simulMethod = qSel("#simul-method").value;
    qSelAll("." + simulMethod + "-param").forEach(function(item){
      item.style.display = "";
    });
  }

  // show/hide parameters for a specific thermostat type
  function onChangeThermostatType()
  {
    // hide all thermostat parameters first
    qSelAll(".thstat-param").forEach(function(item){
      item.style.display = "none";
    });

    // enable the specific thermostat parameters
    let thtype = qSel("#thermostat-type").value;
    qSelAll("." + thtype + "-param").forEach(function(item){
      item.style.display = "";
    });
  }


  function adjustCanvasSize()
  {
    let el = document.querySelector(".canvas");
    let w = el.offsetWidth;
    el.style.height = w + "px";
    // set the size of the context
    el.width = w;
    el.height = w;
  }

  function getDefaultDygraphOption()
  {
    let options = {
      axisLabelFontSize: 10,    
    };

    // determine the width of dygraph plots
    let winWidth = window.innerWidth || screen.width;
    const dygraphMargin = 10*2; // *2 for left and right
    if (winWidth >= 1000) { // plots are float to the right
      options.width = 480;
    } else if (winWidth >= 600) { // two columns
      options.width = Math.floor(winWidth/2 - dygraphMargin);
    } else { // a single colume
      options.width = winWidth - dygraphMargin;
    }
  
    return options;
  }

  function initFolders()
  {
    qSelAll(".folder").forEach(function(item){
      let target = qSel(item.dataset.target);
      let initHide = item.dataset.initHide;
      //console.log(target, initHide, item.dataset);

      // hide the target initially
      if (initHide == "true") {
        target.style.display = "none";
      }

      item.addEventListener("click", function(){
        if (target.style.display === "none") {
          target.style.display = "";
        } else {
          target.style.display = "none";
        }
      });
    });
  }

  function addQrCode()
  {
    // generate a QR code of the page
    let pageUrl = "https://codepen.io/3ki5tj/full/ZEJdeXV";

	// dynamically change the page url
	let m = location.href.match(/cdpn\.io\/(.*)\/full.*\/(.*)/);
	if (m) {
		pageUrl = "https://codepen.io/" + m[1] + "/full/" + m[2];
		console.log(pageUrl);
	}

    setTimeout(function() {
      let url = 'https://api.qrserver.com/v1/create-qr-code/?data='
              + encodeURIComponent(pageUrl)
              + '&size=200x200';
      qSel("#qr-page").src = url;
    }, 1000);
  }

  function init()
  {
    qSel("#simul-method").addEventListener("change", function(){
      onChangeSimulMethod();
    });
    onChangeSimulMethod();

    qSel("#thermostat-type").addEventListener("change", function(){
      onChangeThermostatType();
    });
    onChangeThermostatType();

    window.addEventListener("resize", function(){
      adjustCanvasSize();
    });
    adjustCanvasSize();

    initFolders();

    addQrCode();
  
    // disable Dygraph's touch motions on smartphones
    // https://stackoverflow.com/a/11365986
    Dygraph.defaultInteractionModel.touchend = function() {};
    Dygraph.defaultInteractionModel.touchmove = function() {};
    Dygraph.defaultInteractionModel.touchstart = function() {};
  }
  
  return {
    init: init,
    getDefaultDygraphOption: getDefaultDygraphOption,
  };
})();
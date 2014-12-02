
var datarr = [];



function rotate(x, y, th)
{
  var c = Math.cos(th), s = Math.sin(th);
  //return [ c*x + s*y, -s*x + c*y ];
  return [x, y];
}



function mkpoints()
{
  var ngaus = get_int("ngaus", 1);

  datarr = [];
  for ( var i = 0; i < ngaus; i++ ) {
    var i1 = i + 1;
    var xc = get_float("xc_" + i1, 0.0);
    var yc = get_float("yc_" + i1, 0.0);
    var a = get_float("a_" + i1, 1.0);
    var b = get_float("b_" + i1, 1.0);
    var theta = get_float("theta_" + i1, 0.0);
    var th = theta * Math.PI / 180.0;
    var npt = get_int("npt_" + i1, 100);
    for ( var j = 0; j < npt; j++ ) {
      var xy = rotate(a*gaussrand(), b*gaussrand(), th);
      datarr.push( [xc + xy[0], yc + xy[1]] );
    }
  }
}



/* approximately draw an ellipse */
function drawEllipse(ctx, xc, yc, w, h, theta) {
  var th = theta * Math.PI / 180;
  var a = rotate(0, -h/2, th);
  var b = rotate(w*2/3, -h/2, th);
  var c = rotate(w*2/3, +h/2, th);
  ctx.beginPath();
  ctx.moveTo(xc + a[0], yc + a[1]);
  ctx.bezierCurveTo(xc + b[0], yc + b[1], xc + c[0], yc + c[1], xc - a[0], yc - a[1]);
  ctx.bezierCurveTo(xc - b[0], yc - b[1], xc - c[0], yc - c[1], xc + a[0], yc + a[1]);
  ctx.closePath();
  ctx.stroke();
  ctx.fill();
}



function gengaus()
{
  var c = grab("demobox");
  var ctx = c.getContext("2d");
  var w = c.width;
  var h = c.height;
  var r = Math.min(w/2, h/2) - 1;

  // draw the background
  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect(0, 0, w, h);

  var s = 100.0; // real to screen
  ctx.fillStyle = "#c0c0c0";
  ngaus = get_int("ngaus");
  for ( var i = 0; i < ngaus; i++ ) {
    var i1 = i + 1;
    var xc = get_float("xc_" + i1, 0.0);
    var yc = get_float("yc_" + i1, 0.0);
    var a = get_float("a_" + i1, 1.0);
    var b = get_float("b_" + i1, 1.0);
    var theta = get_float("theta_" + i1, 0.0);
    drawEllipse(ctx, w/2 + xc*s, h/2 + yc*s, a*s, b*s, theta);
  }

  // draw random dots
  ctx.fillStyle = "black";
  mkpoints();

  for ( var i = 0; i < datarr.length; i++ ) {
    var x = datarr[i][0], y = datarr[i][1];
    // a dot === a one-by-one rectangle
    ctx.fillRect( w/2 + x * s, h/2 + y * s, 1, 1);
  }
}



function change_params()
{
  var ngaus = grab("ngaus").value;
  if ( !is_int(ngaus) ) return;
  ngaus = parseInt(ngaus);
  var tab = grab("gausTable");
  var tbody = tab.lastChild;
  var rowid, row, td;

  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    rowid = "gaus_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = ngaus; i < nrows; i++ ) {
    rowid = "gaus_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create non-existing elements
  for ( var i = nrows; i < ngaus; i++ ) {
    var i1 = "" + (i+1);
    rowid = "gaus_row_" + i1;
    row = document.createElement("tr");

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // x_c
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0.0" id="xc_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // y_c
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0.0" id="yc_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // a
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="1.0" id="a_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // b
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="1.0" id="b_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // theta
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="theta_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // npt
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="300" id="npt_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);
  }

}



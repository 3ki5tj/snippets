"use strict";

var seq_g = [];
var atomls_g = [];
var length_g = 30;
var pdbout_g = "";
var emout_g = "";
var equilout_g = "";

// adding the function trim() to strings
if (typeof(String.prototype.trim) === "undefined") {
  String.prototype.trim = function() {
    return String(this).replace(/^\s+|\s+$/g, '');
  };
}

// grab sequence
function readseq(s)
{
  // to upper case
  s = s.toUpperCase();
  // remove leading and trailing spaces
  s = s.replace(/^\s+|\s+$/g, '');
  // convert spaces to tabs
  s = s.replace(/\s+/g, "\t");
  var arr = s.split("\t");
  var i, out = [];
  var aa_1to3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR"
  };
  var aa_names = ["ALA", "CYS", "CYX", "CYS2", "ASP", "ASPP", "ASN", "ASH",
    "GLU", "GLH", "GLUP", "PHE", "GLY",
    "HIS", "HID", "HIE", "HIP", "HS2", "HSD", "HSE", "HSP",
    "ILE", "LYS", "LYN", "LEU", "MET", "ASN",
    "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"];

  // detect if the format is 1-letter or 3-letter
  for ( i = 0; i < arr.length; i++ ) {
    var tok = arr[i];
    if ( tok.length === 4 ) {
      if ( aa_names.indexOf(tok) >= 0 )
        continue;
      var ter = tok.slice(0, 1);
      if ( ter === "C" || ter === "N" ) {
        tok = tok.slice(-3);
      }
    }
    if ( tok.length !== 3 || aa_names.indexOf(tok) < 0 )
      break;
  }

  if ( i === arr.length ) { // 3-letter sequence
    for ( var i = 0; i < arr.length; i++ ) {
      var ch = arr[i].slice(0, 1); // first letter
      if ( arr[i].length === 4 && (ch === "N" || ch === "C") )
        arr[i] = arr[i].slice(-3);
    }
  } else { // 1-letter sequence
    s = s.replace(/\s+/g, "");
    arr = [];
    for ( i = 0; i < s.length; i++ ) {
      arr.push( aa_1to3[ s[i] ] );
    }
    console.log( arr );
  }

  return arr;
}

function mkatom(atomid, atomname, resname, resid, x, ele)
{
  atomid = ("   " + atomid).slice(-4);
  atomname = (atomname + "  ").slice(0, 3);
  resname = (resname + "    ").slice(0, 4);
  resid = ("   " + resid).slice(-4);
  var x0 = ("     " + x[0].toFixed(3)).slice(-8);
  var x1 = ("     " + x[1].toFixed(3)).slice(-8);
  var x2 = ("     " + x[2].toFixed(3)).slice(-8);
  var ele = atomname.slice(0, 1);
  return "ATOM   " + atomid + "  " + atomname + " "
       + resname + " " + resid + '    ' + x0 + x1 + x2
       + "  1.00  1.00          " + ele + "\n";
}

function pushatom(atomls, res)
{
  // we have to deep copy the coordinates
  atomls.push( [res[0], res[1], res[2].slice()] );
}

function mkter(atomid, resname, resid)
{
  atomid = ("   " + atomid).slice(-4);
  resname = (resname + "    ").slice(0, 4);
  resid = ("   " + resid).slice(-4);
  return "TER    " + atomid + "      " + resname + " " + resid + "\n";
}

function gchoose(gd, g1, g2, g3, rotamer)
{
  var arr = [gd, g1, g2, g3];
  return arr[rotamer];
}

// shift the center of structure
function shiftcenter(atomls)
{
  var n = atomls.length, xmin, xmax, x0, x1, r0 = 100000, r1 = -10000;

  for ( var i = 0; i < n; i++ ) {
    var x = atomls[i][2];
    var resid = atomls[i][1];
    if ( resid < r0 ) r0 = resid;
    if ( resid > r1 ) r1 = resid;
    if ( ["CA", "CH3", "CAY", "CAT"].indexOf(atomls[i][0]) >= 0 ) {
      if ( x0 === undefined ) x0 = x.slice(0);
      x1 = x.slice(0);
    }
    if ( i == 0 ) {
      xmin = x.slice(0);
      xmax = x.slice(0);
      continue;
    }
    for ( var d = 0; d < 3; d++ ) {
      if ( x[d] < xmin[d] ) xmin[d] = x[d];
      if ( x[d] > xmax[d] ) xmax[d] = x[d];
    }
  }
  var xc = [0, 0, 0], sz = [0, 0, 0];
  for ( d = 0; d < 3; d++ ) {
    xc[d] = (xmin[d] + xmax[d]) / 2;
    sz[d] = xmax[d] - xmin[d];
  }
  // shift the center
  for ( var i = 0; i < n; i++ )
    vdec(atomls[i][2], xc);

  document.getElementById("info").innerHTML = "Info: " +
    (r1 - r0 + 1) + " residues; " +
    n + " atoms. " +
    "Dimension: "
    + "<i>x</i>: " + sz[0].toFixed(2) + "&#8491;, "
    + "<i>y</i>: " + sz[1].toFixed(2) + "&#8491;, "
    + "<i>z</i>: " + sz[2].toFixed(2) + "&#8491;. "
    + "End-to-end distance: "
    + vdist(x0, x1).toFixed(2) + "&#8491;.";

  return sz;
}

function mkpdb(seq)
{
  var nter = document.getElementById("n-caps").value;
  var cter = document.getElementById("c-caps").value;
  var format = document.getElementById("format").value;
  var resId0 = parseInt(document.getElementById("res-id-start").value);

  var nres = seq.length;
  if ( nter === "ACE" ) {
    seq.unshift("ACE");
  } else {
    seq.unshift("");
  }

  var D2R = Math.PI/180;
  var R2D = 180/Math.PI;
  // standard bond lengths in angstroms
  var B_CAC     = 1.53;
  var B_CN_PEP  = 1.33;
  var B_NCA     = 1.46;
  var B_CC      = 1.54;
  var B_CN      = 1.48; // peptide bond
  var B_CO      = 1.22; // carbonyl
  var B_CS      = 1.81
  var BH_OHOH   = 2.80
  var BH_NHOH   = 2.90
  var BH_OHOC   = 2.80
  var B_CC_RING = 1.40

  // cosine and sine values
  var c12 = Math.cos(12 * D2R), s12 = Math.sin(12 * D2R);
  var c30 = Math.cos(30 * D2R), s30 = Math.sin(30 * D2R);
  var c36 = Math.cos(36 * D2R), s36 = Math.sin(36 * D2R);
  var c72 = Math.cos(72 * D2R), s72 = Math.sin(72 * D2R);

  var rotate = parseFloat( document.getElementById("rotate").value );
  var swing = parseFloat( document.getElementById("swing").value );
  var rise = parseFloat( document.getElementById("rise").value );

  var rotang, swgang, risang;
  var useDefaultAngles = false;
  if ( useDefaultAngles ) {
    rotang = 110.0 / Math.sqrt(nres) * D2R;
    swgang = 50.0 * D2R;
    risang = (38.0 / Math.sqrt(nres) + 81.0 / nres) * D2R;
  } else {
    rotang = rotate * D2R;
    swgang = swing * D2R;
    risang = rise * D2R;
  }

  var cr = Math.cos(risang), sr = Math.sin(risang);
  var q = Math.cos(swgang) * cr * cr + sr * sr - 1.0/3;
  q /= 1 + Math.cos(swgang);
  if ( q < 1e-6 ) q = 1e-6;
  var vang = Math.asin( Math.sqrt(q) );

  var thp = risang + vang; // for even-index residues
  var c1p = Math.cos(thp), s1p = Math.sin(thp);
  var c2p = Math.cos(thp - Math.PI/3), s2p = Math.sin(thp - Math.PI/3);
  var thm = risang - vang; // for odd-index residues
  var c1m = Math.cos(thm), s1m = Math.sin(thm);
  var c2m = Math.cos(thm + Math.PI/3), s2m = Math.sin(thm + Math.PI/3);

  var phi = 0.5 * Math.acos(-1.0/3);
  var c3 = Math.cos(phi), s3 = Math.sin(phi);

  console.log("angles", rotang, swgang, risang, "res", nres,
      "theta+:", thp, "theta-", thm, "vang", vang, "q", q);

  var resid = 0;
  var resnm = "";
  var n = seq.length;
  var xyang = 0;
  var os = [0, 0, 0];
  var atomls = [];
  var dir_nca = [0, 0, 0];
  var sgn;
  var c1, s1, c2, s2, cp, sp;
  var xca = [0, 0, 0];
  var xc = [0, 0, 0];
  var xo = [0, 0, 0];
  var xn = [0, 0, 0];
  var xnp = [0, 0, 0];
  var xcb = [0, 0, 0];
  var xg = [0, 0, 0], xz = [0, 0, 0], xh = [0, 0, 0];
  var xog1 = [0, 0, 0], xcg1 = [0, 0, 0];
  var xog2 = [0, 0, 0], xcg2 = [0, 0, 0];
  var xog3 = [0, 0, 0], xcg3 = [0, 0, 0];
  var xd1 = [0, 0, 0], xd2 = [0, 0, 0], xd = [0, 0, 0];
  var xe1 = [0, 0, 0], xe2 = [0, 0, 0], xe = [0, 0, 0];
  var dir_cac = [0, 0, 0], dir_nca = [0, 0, 0], dx = [0, 0, 0];
  var u = [0, 0, 0], v = [0, 0, 0], w = [0, 0, 0], p = [0, 0, 0], q = [0, 0, 0];

  var imax = n;
  // we need the coordinates of the next CA if the C-terminal is NMET
  if ( cter === "NMET" ) imax += 1;

  for ( var i = 0; i < imax; i++ ) {
    if ( i % 2 == 0 ) {
      sgn = 1;
      c1 = c1p, s1 = s1p, c2 = c2p, s2 = s2p;
    } else {
      sgn = -1;
      c1 = c1m, s1 = s1m, c2 = c2m, s2 = s2m;
    }

    // build the coordinates of the backbone atoms, C, CA and O
    xyang += rotang;
    cp = Math.cos(xyang);
    sp = Math.sin(xyang);

    xca = os.slice(0); // clone: xca = os

    // position of C from CA
    dir_cac = [c1 * cp, c1 * sp, s1];
    vsadd(xc, xca, dir_cac, B_CAC);

    // position of O from C
    dx = [-sr * cp, -sr * sp, cr];
    vsadd(xo, xc, dx, B_CO * sgn);

    if ( i >= n ) break;

    if ( i > 0 ) {
      // add CB
      // p is pointing to the side of CA opposite to N and C
      // q is perpendicular q
      // w is the position of CB that makes residue left-handed
      // w is the direction CA->CB
      vnormalize( vdiff(p, dir_nca, dir_cac) );
      vnormalize( vcross3d(q, dir_cac, dir_nca) );
      w = [p[0]*c3 + q[0]*s3, p[1]*c3 + q[1]*s3, p[2]*c3 + q[2]*s3];
      vnormalize(w);
      vsadd(xcb, xca, w, B_CC);

      // compute the three CG positions
      // G1 position, opposite to C in the C=0 group
      // viewed along the direction of CB-CA
      vnormalize( vdiff(u, xca, xc) );
      vsadd(xog1, xcb, u, B_CO);
      vsadd(xcg1, xcb, u, B_CC);
      // G2 position, opposite to N in the N-H group
      vnormalize( vdiff(v, xca, xn) );
      vsadd(xog2, xcb, v, B_CO);
      vsadd(xcg2, xcb, v, B_CC);
      // G3 position, p, q and u are perpendicular to w
      vperpen(p, u, w);
      vperpen(q, v, w);
      vadd(u, p, q);
      // v is the direction of CB -> CG3
      v = [-u[0]+w[0]/3, -u[1]+w[1]/3, -u[2]+w[2]/3];
      vnormalize(v);
      vsadd(xog3, xcb, v, B_CO);
      vsadd(xcg3, xcb, v, B_CC);
    } else {
    }

    xnp = xn.slice(0);
    // position of N from C
    dx = [c2 * cp, c2 * sp, s2];
    vsadd(xn, xc, dx, B_CN_PEP);
    dir_nca = dir_cac.slice(0);

    // set the origin for the next CA
    vsadd(os, xn, dir_nca, B_NCA);

    xyang += swgang * sgn;
    cp = Math.cos(xyang);
    sp = Math.sin(xyang);

    resnm = seq[i];

    if ( i > 0 ) {
      pushatom(atomls, ["N", resid, xnp] );
      pushatom(atomls, ["CA", resid, xca] );
    } else if ( resnm !== "" ) { // ACE
      var atnm = (format === "CHARMM") ? "CAY" : "CH3";
      pushatom(atomls, [atnm, resid, xca] );
    }

    if ( i > 0 || resnm !== "" ) {
      //console.log( i, resnm, xc, xo);
      var cname = "C", oname = "O";
      if (format === "CHARMM" && i === 0) {
        cname = "CY";
        oname = "OY";
      }
      pushatom(atomls, [cname, resid, xc.slice(0)] );
      if ( i === n - 1 ) {
        if ( cter === "" ) {
          var oc1 = "OC1", oc2 = "OC2";
          if ( format === "CHARMM" ) {
            oc1 = "OT1";
            oc2 = "OT2";
          } else if ( format === "AMBER" ) {
            oc1 = "O";
            oc2 = "OXT";
          }
          pushatom(atomls, [oc1, resid, xo] );
          vdiff(u, xc, xo);
          vnormalize(u);
          vdiff(v, xc, xca);
          vnormalize(v);
          dx = [B_CO*(u[0]+v[0]), B_CO*(u[1]+v[1]), B_CO*(u[2]+v[2])];
          vadd(p, xc, dx);
          pushatom(atomls, [oc2, resid, p] );
        } else {
          pushatom(atomls, [oname, resid, xo] );
        }
      } else {
        pushatom(atomls, [oname, resid, xo] );
      }
    }

    if ( resnm === "" || resnm === "GLY" ) {
      resid += 1;
      continue;
    } else if ( resnm === "ACE" ) {
      if ( format !== "CHARMM" ) {
        resid += 1;
      }
      continue;
    }

    pushatom(atomls, ["CB", resid, xcb] );

    var rotamer = 0;
    if ( resnm === "" || resnm === "ALA" ) {
      ;
    var xg = [0, 0, 0];
    } else if ( resnm === "SER" ) {
      xg = gchoose(xog2, xog1, xog2, xog3, rotamer);
      pushatom(atomls, ["OG", resid, xg] );

    } else if ( resnm.slice(0, 2) === "CY" ) {
      xg = gchoose(xog2, xog1, xog2, xog3, rotamer);
      pushatom(atomls, ["SG", resid, xg] );

    } else if ( resnm === "VAL" ) {
      pushatom(atomls, ["CG1", resid, xcg1] );
      pushatom(atomls, ["CG2", resid, xcg2] );

    } else if ( resnm === "LEU" ) {
      xg = gchoose(xcg1, xcg1, xcg2, xcg3, rotamer);
      var xdc = [0, 0, 0], xd23 = [0, 0, 0];
      vdiff(q, xcb, xca);
      vadd(xd1, xg, q);
      vdiff(u, xg, xcb);
      vsadd(v, q, u, -1./3);
      vsadd(xdc, xg, u, 1./3);
      vsadd(xd23, xdc, v, -1./2);
      vnormalize( vcross3d(w, q, u) );
      vsadd(xd2, xd23, w, Math.sqrt(2.0/3) * B_CC);
      pushatom(atomls, ["CG", resid, xg] );
      pushatom(atomls, ["CD1", resid, xd1] );
      pushatom(atomls, ["CD2", resid, xd2] );

    } else if ( resnm === "ILE" ) {
      var xg1 = xcg1, xg2 = xcg2;
      if ( rotamer === 3 ) {
        xg1 = xcg3;
        xg2 = xcg1;
      }
      vdiff(dx, xcb, xca);
      vadd(xd, xg1, dx);
      pushatom(atomls, ["CG1", resid, xg1] );
      pushatom(atomls, ["CD",  resid, xd]  );
      pushatom(atomls, ["CG2", resid, xg2] );

    } else if ( resnm === "TRP" ) {
      xg = gchoose(xcg1, xcg1, xcg2, xcg3, rotamer);
      vnormalize( vdiff(v, xg, xcb) );
      vdiff(dx, xca, xcb);
      vnormalize( vcross3d(u, v, dx) );
      if ( xg === xcg2 ) {
        u = vneg(u);
      }
      var x8 = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],
                [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]];
      vsadd(x8[0], xg,    vlincomb2(dx, u, v,  c36,  s36), B_CC_RING);
      vsadd(x8[1], xg,    vlincomb2(dx, u, v, -c36,  s36), B_CC_RING);
      vsadd(x8[2], x8[0], vlincomb2(dx, u, v, -c72,  s72), B_CC_RING);
      vsadd(x8[3], x8[1], vlincomb2(dx, u, v,  c72,  s72), B_CC_RING);
      vsadd(x8[4], x8[1], vlincomb2(dx, u, v, -c12, -s12), B_CC_RING);
      vsadd(x8[5], x8[3], vlincomb2(dx, u, v, -c72,  s72), B_CC_RING);
      vsadd(x8[7], x8[5], vlincomb2(dx, u, v, -c12, -s12), B_CC_RING);
      vsadd(x8[6], x8[7], vlincomb2(dx, u, v, -c72, -s72), B_CC_RING);
      var at8 = ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"];
      pushatom(atomls, ["CG", resid, xg] );
      for ( var k = 0; k < 8; k++ ) {
        pushatom(atomls, [at8[k], resid, x8[k]] );
      }

    } else if ( resnm.slice(0, 2) === "HI" ) {
      xg = gchoose(xcg2, xcg1, xcg2, xcg3, rotamer);
      vnormalize( vdiff(v, xg, xcb) );
      vnormalize( vcross3d(w, vdiff(dx, xcb, xca), v) );
      vadd(xd1, xg,  vlincomb2(dx, v, w, s36 * B_CC_RING,  c36 * B_CC_RING));
      vadd(xe1, xd1, vlincomb2(dx, v, w, s72 * B_CC_RING, -c72 * B_CC_RING));
      vadd(xd2, xg,  vlincomb2(dx, v, w, s36 * B_CC_RING, -c36 * B_CC_RING));
      vadd(xe2, xd2, vlincomb2(dx, v, w, s72 * B_CC_RING,  c72 * B_CC_RING));
      pushatom(atomls, ["CG",  resid, xg ] );
      pushatom(atomls, ["ND1", resid, xd1] );
      pushatom(atomls, ["CE1", resid, xe1] );
      pushatom(atomls, ["CD2", resid, xd2] );
      pushatom(atomls, ["NE2", resid, xe2] );

    } else if ( resnm === "PHE" || resnm === "TYR" ) {
      xg = gchoose(xcg2, xcg1, xcg2, xcg3, rotamer);
      vnormalize( vdiff(v, xg, xcb) );
      vdiff(dx, xcb, xca);
      vnormalize( vcross3d(w, dx, v) );
      vlincomb2(dx, v, w, s30 * B_CC_RING, c30 * B_CC_RING);
      vadd(xd1, xg, dx);
      vsadd(xe1, xd1, v, B_CC_RING);
      var xz = [0, 0, 0], xh = [0, 0, 0];
      vsadd(xz, xg, v, 2 * B_CC_RING);
      vsadd(xh, xz, v, B_CO);
      vlincomb2(dx, v, w, s30 * B_CC_RING, -c30 * B_CC_RING);
      vadd(xd2, xg, dx);
      vsadd(xe2, xd2, v, B_CC_RING);

      pushatom(atomls, ["CG",  resid, xg ] );
      pushatom(atomls, ["CD1", resid, xd1] );
      pushatom(atomls, ["CE1", resid, xe1] );
      pushatom(atomls, ["CZ",  resid, xz ] );
      if ( resnm === "TYR" ) {
        pushatom(atomls, ["OH", resid, xh] );
      }
      pushatom(atomls, ["CD2", resid, xd2] );
      pushatom(atomls, ["CE2", resid, xe2] );

    } else if ( resnm === "THR" ) {
      pushatom(atomls, ["OG1", resid, xog1] );
      pushatom(atomls, ["CG2", resid, xcg2] );

    } else if ( resnm === "GLU" || resnm === "GLN" ) {
      xg = gchoose(xcg2, xcg1, xcg2, xcg3, rotamer);
      vdiff(v, xcb, xca);
      var xcd = [0, 0, 0];
      vadd(xcd, xg, v);
      vnormalize(v);
      vnormalize( vcross3d(u, v, vdiff(dx, xcb, xg) ) );
      vlincomb2(dx, u, v,  c30, s30);
      vsadd(xe1, xcd, dx, B_CO);
      vlincomb2(dx, u, v, -c30, s30);
      vsadd(xe2, xcd, dx, B_CO);
      pushatom(atomls, ["CG",  resid, xg ] );
      pushatom(atomls, ["CD",  resid, xcd] );
      pushatom(atomls, ["OE1", resid, xe1] );
      var atme2 = (resnm === "GLU") ? "OE2" : "NE2";
      pushatom(atomls, [atme2, resid, xe2] );

    } else if ( resnm === "ASP" || resnm === "ASN" ) {
      xg = gchoose(xcg1, xcg1, xcg2, xcg3, rotamer);
      vnormalize( vdiff(v, xg, xcb) );
      vperpen(u, vdiff(dx, xca, xcb), v);
      vlincomb2(dx, u, v,  c30, s30);
      vsadd(xd1, xg, dx, B_CO);
      vlincomb2(dx, u, v, -c30, s30);
      vsadd(xd2, xg, dx, B_CO);
      pushatom(atomls, ["CG",  resid, xg ] );
      pushatom(atomls, ["OD1", resid, xd1] );
      var atmd2 = (resnm === "ASP") ? "OD2" : "ND2";
      pushatom(atomls, [atmd2, resid, xd2] );

    } else if ( resnm === "LYS" ) {
      xg = gchoose(xcg2, xcg1, xcg2, xcg3, rotamer);
      vadd(xd, xg, vdiff(dx, xcb, xca));
      vadd(xe, xd, vdiff(dx, xg,  xcb));
      vnormalize( vdiff(dx, xcb, xca) );
      vsadd(xz, xe, dx, B_CN);
      pushatom(atomls, ["CG", resid, xg] );
      pushatom(atomls, ["CD", resid, xd] );
      pushatom(atomls, ["CE", resid, xe] );
      pushatom(atomls, ["NZ", resid, xz] );

    } else if ( resnm === "ARG" ) {
      var xh1 = [0, 0, 0], xh2 = [0, 0, 0];
      xg = gchoose(xcg2, xcg1, xcg2, xcg3, rotamer);
      vdiff(q, xcb, xca);
      vadd(xd, xg, q);
      vnormalize( vdiff(v, xg, xcb) );
      vsadd(xe, xd, v, B_CN);
      vnormalize( vcross3d(w, v, q) );
      vlincomb2(dx, w, v, c30,  s30);
      vsadd(xz, xe, dx, B_CN);
      vlincomb2(dx, w, v, c30, -s30);
      vsadd(xh1, xz, dx, B_CN);
      vsadd(xh2, xz, v, B_CN);
      pushatom(atomls, ["CG",  resid, xg ] );
      pushatom(atomls, ["CD",  resid, xd ] );
      pushatom(atomls, ["NE",  resid, xe ] );
      pushatom(atomls, ["CZ",  resid, xz ] );
      pushatom(atomls, ["NH1", resid, xh1] );
      pushatom(atomls, ["NH2", resid, xh2] );

    } else if ( resnm === "MET" ) {
      xg = gchoose(xcg1, xcg1, xcg2, xcg3, rotamer);
      vnormalize( vdiff(dx, xcb, xca) );
      vsadd(xd, xg, dx, B_CS);
      vnormalize( vdiff(dx, xg,  xcb) );
      vsadd(xe, xd, dx, B_CS);
      pushatom(atomls, ["CG", resid, xg] );
      pushatom(atomls, ["SD", resid, xd] );
      pushatom(atomls, ["CE", resid, xe] );

    } else if ( resnm === "PRO" ) {
      vnormalize( vdiff(u, xcb, xnp) );
      vdiff(p, xcb, xca);
      vdiff(q, xnp, xca);
      vnormalize( vadd(v, p, q) );
      vlincomb2(dx, u, v, -c72, s72);
      vsadd(xd, xcb, dx, B_CC);
      vlincomb2(dx, u, v,  c72, s72);
      vsadd(xg, xnp, dx, B_CC);
      pushatom(atomls, ["CD", resid, xd] );
      pushatom(atomls, ["CG", resid, xg] );

    } else {
      console.log("unknown residue", resnm);
    }

    resid += 1;
  }

  if ( format === "CHARMM" ) { // CHARMM
    if ( cter === "NH2" || cter === "NMET" ) {
      resid -= 1;
      if ( cter === "NH2" ) {
        pushatom(atomls, ["NT", resid, xn] );
      } else {
        pushatom(atomls, ["NT", resid, xn] );
        pushatom(atomls, ["CAT", resid, xca] );
      }
      seq.push( resnm );
      resid += 1;
    }
  } else if ( format === "CHARMM-GMX" ) { // CHARMM-GMX
    if ( cter === "NH2" || cter === "NMET" ) {
      if ( cter === "NH2" ) {
        // don't do anything, force field does not have CT2
        //pushatom(atomls, ["N", resid, xn] );
        //seq.push( "CT2" );
      } else { // NMET
        pushatom(atomls, ["N", resid, xn] );
        pushatom(atomls, ["CH3", resid, xca] );
        seq.push( "CT3" );
        resid += 1;
      }
    }
  } else { // AMBER-GMX
    if ( cter === "NH2" || cter === "NMET" ) {
      if ( cter === "NH2" ) {
        pushatom(atomls, ["N", resid, xn] );
        seq.push( "NH2" );
      } else { // NMET
        pushatom(atomls, ["N", resid, xn] );
        pushatom(atomls, ["CH3", resid, xca] );
        seq.push( "NME" );
      }
      resid += 1;
    }
  }

  // shift the coordinates to make them nonnegative
  var sz = shiftcenter(atomls);

  // write PDB
  var boxsize = parseFloat( document.getElementById("boxsize").value );
  var src = "";
  resid = 0;

  var resIdOffset = 0;
  if ( nter !== "" ) {
    resIdOffset = 1;
  }
  resIdOffset += resId0 - 1;

  var nn = atomls.length;
  var xyz = [], x = [];
  for ( var k = 0; k < nn; k++ ) {
    var trio = atomls[k];
    resid = trio[1];
    resnm = seq[resid];
    if ( resid === 0 && nter !== "" && format === "CHARMM" ) {
      resnm = seq[1];
    }
    resid += resIdOffset;
    var x0 = trio[2];
    for ( var d = 0; d < 3; d++ )
      x[d] = x0[d] + boxsize * 0.5;
    src += mkatom(k + 1, trio[0], resnm, resid, x);
    xyz.push( x0 );
  }
  src += mkter(k + 1, resnm, resid + resIdOffset);

  // write script
  var script = "", runem = "", runequil = "", emname = "", equilname = "";
  var outname = document.getElementById("outname").value;
  var out0name = outname + "_raw.pdb";
  document.getElementById("out0name").innerHTML = out0name;
  var nstemin = document.getElementById("nstemin").value;
  var nstequil = document.getElementById("nstequil").value;
  var out1 = outname, out2 = "", cellsrc = "";
  var watermodel = (boxsize > 0) ? "tip3p" : "none";

  if ( format === "CHARMM" ) {
    script += "package require psfgen\n" +
              "topology top_all27_prot_lipid.inp\n" +
              "pdbalias residue HIS HSE\n" +
              "pdbalias atom ILE CD1 CD\n";
    script += "segment U {pdb " + out0name;
    if ( nter === "ACE" ) {
      script += "\nfirst ACE";
    }
    if ( cter === "NH2" ) {
      script += "\nlast CT2";
    } else if ( cter === "NMET" ) {
      script += "\nlast CT3";
    }
    script += "}\n";
    script += "coordpdb " + out0name + " U\n" +
              "guesscoord\n" +
              "writepdb " + outname + ".pdb\n" +
              "writepsf " + outname + ".psf\n";

  } else if ( format === "AMBER" ) {
    script += "# To use this script on the command line:\n"
           +  "# tleap -f " + document.getElementById("scriptname").innerHTML.trim() + "\n";
    script += "source leaprc.protein.ff14SB\n";
    script += 'mol = loadpdb "' + out0name + '"\n';
    if ( boxsize ) {
      var bs = boxsize.toFixed(1);
      script += "source leaprc.water.tip3p\n";
      script += "set default FlexibleWater on\n";
      script += "solvate mol TIP3PBOX {" + bs + " " + bs + " " + bs + "}\n";
      out1 += "_wb";
    }
    script += "savepdb mol " + out1 + ".pdb\n";
    script += "saveamberparm mol " + out1 + ".prmtop " + out1 + ".inpcrd\n";
    script += "quit\n";

  } else if ( format === "CHARMM-GMX" ) {
    script += "gmx pdb2gmx -f " + out0name + " -o "
            + outname + ".gro -ff charmm27 -water " + watermodel + " -ignh -ter\n";
    if ( cter === "" || cter === "NMET" ) {
      script += "# select None for both terminals\n";
    } else if ( cter === "NH2" ) {
      script += "# select None for N-terminal, CT2 for C-terminal\n";
    }

  } else if ( format === "AMBER-GMX" ) {
    var amberver = document.getElementById("amberver").value;
    script += "gmx pdb2gmx -f " + out0name + " -o "
            + outname + ".gro -ff " + amberver + " -water " + watermodel + " -ignh\n";
  }

  if ( format === "CHARMM"   // NAMD running script
    || format === "AMBER" ) {
    emname = "em.conf";
    equilname = "equil.conf";

    // commands for solvation
    if ( boxsize > 0 ) {
      out1 = outname + "_wb";
      var boxs = ("   " + boxsize.toFixed(1)).slice(-5);
      var boxh = (boxsize*0.5).toFixed(1);
      if ( format === "CHARMM" ) {
        script += "\n# adding water into the box\n"
        script += "package require solvate\n";
        script += "solvate " + outname + ".psf " + outname + ".pdb" +
         " -minmax {{0 0 0} {" + boxs + " " + boxs + " " + boxs +
         "}} -o " + out1 + "\n";
      }

      cellsrc +=
       "cellBasisVector1   " + boxs + "   0.0   0.0\n" +
       "cellBasisVector2     0.0 " + boxs + "   0.0\n" +
       "cellBasisVector3     0.0   0.0 " + boxs + "\n" +
       "cellOrigin          " + boxh + "  " + boxh + "  " + boxh + "\n";
      cellsrc +=
       "PME                 yes\n" +
       "PMEGridSpacing      1.0\n";
    }
    out2 = outname + "_init";
    runem = "";
    if ( format === "CHARMM" ) {
      runem += "paraTypeCharmm      on\n" +
               "parameters          par_all27_prot_lipid.inp\n";
    } else {
      runem += "amber               on\n" +
               "parmfile            " + out1 + ".prmtop\n";
    }
    runem += "" +
     "set inname          " + out1 + "\n" +
     "set outname         " + out2 + "\n" +
     "set temp            300\n" +
     "structure           " + out1 + ".psf\n" +
     "coordinates         $inname.pdb\n" +
     "temperature         $temp\n" +
     "exclude             scaled1-4\n" +
     "1-4scaling          1.0\n" +
     "switching           on\n" +
     "timestep            2.0\n" +
     "rigidBonds          all\n" +
     "nonbondedFreq       1\n" +
     "fullElectFrequency  1\n" +
     "stepspercycle       4\n";
    if ( boxsize > 0 ) {
      runem += cellsrc;
      runem += "" +
       "switchdist          10.0\n" +
       "cutoff              12.0\n" +
       "pairlistdist        14.0\n";
    } else {
      runem += "# implicit solvent parameters\n" +
       "GBIS                on\n" +
       "ionConcentration    0.2\n" +
       "SASA                on\n" +
       "surfaceTension      0.005\n" +
       "alphaCutoff         14.0\n";
      runem += "" +
       "switchdist          15.0\n" +
       "cutoff              16.0\n" +
       "pairlistdist        18.0\n";
    }
    runem += "" +
     "wrapAll             on\n" +
     "langevin            on\n" +
     "langevinDamping     1.0\n" +
     "langevinHydrogen    off\n" +
     "langevinTemp        $temp\n" +
     "outputname          $outname\n" +
     "binaryoutput        no\n" +
     "binaryrestart       no\n" +
     "DCDfreq             1000000\n" +
     "restartfreq         1000000\n" +
     "xstfreq             1000000\n" +
     "outputEnergies      1000000\n" +
     "outputPressure      1000000\n" +
     "outputTiming        1000000\n" +
     "minimize            " + nstemin + "\n" +
     "run                 " + nstequil + "\n";
    runem = "# To use this script on the command line:\n"
          + "# namd2 " + emname + "\n\n" + runem;

    if ( format === "CHARMM" ) {
      script = "# To use this script on the command line:\n"
             + "#   vmd -dispdev text -eofexit < input.tcl\n\n"
             + script;
    }

  } else { // GROMACS .mdp files
    var out1 = outname;
    if ( boxsize > 0 ) {
      var boxgro = outname + "_box.gro";
      out1 = outname + "_wb";
      script += "\n# adding water to the box\n";
      script += "gmx editconf -f " + outname + ".gro -o "
              + boxgro + " -box " + (boxsize*0.1).toFixed(2) + "\n";
      script += "gmx solvate -cp " + boxgro + " -cs spc216.gro -o " + out1 + " -p topol.top\n";
    }
    var gmxrun = "" +
     "dt = 0.002\n" +
     "nstxtcout = 1000000\n" +
     "nstxout = 0\n" +
     "nstvout = 0\n" +
     "nstfout = 0\n" +
     "nstcalcenergy = 10\n" +
     "nstcomm = 10\n" +
     "nstlog = 1000\n" +
     "nstenergy = 1000\n" +
     "xtc-grps = System\n" +
     "tc-grps = System\n" +
     "energygrps = System\n" +
     "tau-t = 0.1\n" +
     "ref-t = 300\n";
    if ( boxsize > 0 ) { // explicit solvent parameters
      gmxrun += "" +
       "ns-type = grid\n" +
       "nstlist = 10\n" +
       "cutoff-scheme = verlet\n" +
       "rvdw = 1.2\n" +
       "rvdw-switch = 0.9\n" +
       "rlist = 1.2\n" +
       "rcoulomb = 1.2\n" +
       "vdwtype = shift\n" +
       "coulombtype = PME\n" +
       "fourierspacing = 0.144\n" +
       "pme-order = 4\n" +
       "ewald-rtol = 1e-5\n";
    } else { // implicit solvent parameters
      gmxrun += "" +
       "comm-mode = angular\n" +
       "ns-type = simple\n" +
       "nstlist = 0\n" +
       "pbc = no\n" +
       "cutoff-scheme = group\n" +
       "rlist = 0.0\n" +
       "rvdw = 0.0\n" +
       "rcoulomb = 0.0\n" +
       "implicit-solvent = GBSA\n" +
       "gb-algorithm = OBC\n" +
       "rgbradii = 0.0\n" +
       "nstgbradii = 1\n" +
       "vdwtype = Cut-off\n" +
       "coulombtype = Cut-off\n" +
       "sa-surface-tension = 2.25936\n";
    }
    // for energy minimization
    var gmxem = "" +
     "define -DFLEX_SPC\n" +
     "integrator = steep\n" +
     "constraints = none\n" +
     "emstep = 0.01\n" +
     "emtol = 2000\n" +
     "gen-vel = no\n" +
     "nsteps = " + nstemin + "\n" +
     "gen-seed = " + Math.floor(Math.random() * 1000000000 + 1)+ "\n" +
      gmxrun;
    // for the equilibration run
    var gmxequil = "" +
     "integrator = md\n" +
     "constraints = hbonds\n" +
     "nsteps = " + nstequil + "\n" +
     "gen-seed = " + Math.floor(Math.random() * 1000000000 + 1) + "\n" +
      gmxrun;
    runem = "### Energy minimization (em.mdp) #########\n" + gmxem;
    runequil = "### Equilibration run (equil.mdp) ###########\n" + gmxequil;
    emname = "em.mdp";
    equilname = "equil.mdp";

    // energy minimization script
    script += "\n# Energy minimization (first prepare em.mdp from the right box)\n"
            + "gmx grompp -f em.mdp -c " + out1 + " -o em.tpr\n"
            + "gmx mdrun -v -deffnm em -c " + outname + "_em -nt 1\n";

    // equlibration script
    script += "\n# Equilibration run (first prepare equil.mdp from the right box)\n"
            + "gmx grompp -f equil.mdp -c " + outname + "_em -o equil.tpr\n"
            + "gmx mdrun -v -deffnm equil -c " + outname + "_init -nt 1\n";
  }
  document.getElementById("emscript").innerHTML = emname;
  document.getElementById("equilscript").innerHTML = equilname;

  return [src, atomls, sz, script, runem, runequil];
}

function mkspx(refresh)
{
  var format = document.getElementById("format").value;

  //document.getElementById("amberver").disabled = ( format.slice(0, 5) !== "AMBER" );
  document.getElementById("amberver_wrapper").style.visibility
    = ( format === "AMBER-GMX" ) ? "visible" : "hidden";

  if ( format.slice(-3) === "GMX" ) {
    document.getElementById("scriptname").innerHTML = "input.sh";
    document.getElementById("andpart").style.visibility = "visible";
    document.getElementById("emoutput").rows = "7";
    document.getElementById("equiloutput").style.visibility = "visible";
  } else if ( format === "AMBER" ) {
    document.getElementById("scriptname").innerHTML = "input.leap";
    document.getElementById("andpart").style.visibility = "hidden";
    document.getElementById("emoutput").rows = "15";
    document.getElementById("equiloutput").style.visibility = "hidden";
  } else {
    document.getElementById("scriptname").innerHTML = "input.tcl";
    document.getElementById("andpart").style.visibility = "hidden";
    document.getElementById("emoutput").rows = "15";
    document.getElementById("equiloutput").style.visibility = "hidden";
  }

  if ( refresh || atomls_g.length === 0 ) {
    seq_g = readseq( document.getElementById("aainput").value );
    var ret = mkpdb(seq_g);
    pdbout_g = ret[0];
    atomls_g = ret[1];
    emout_g = ret[4];
    equilout_g = ret[5];
    document.getElementById("pdboutput").value = pdbout_g;
    document.getElementById("scriptoutput").value = ret[3];
    document.getElementById("emoutput").value = emout_g;
    document.getElementById("equiloutput").value = equilout_g;
    var sz = ret[2];
    length_g = 0.5 * Math.max(sz[0], sz[1], sz[2]);
  }

  mousescale = parseFloat( document.getElementById("scaleinput").value );

  var ballscale = parseFloat( document.getElementById("ballScaleInput").value );
  console.log(m2str(viewmat), "scale", mousescale, "ball scale", ballscale, "length", length_g);

  var boxsize = parseFloat( document.getElementById("boxsize").value );

  pdbdraw(seq_g, atomls_g, length_g, boxsize,
      "animation-box", mousescale, ballscale, false, false);
}

function paint()
{
  mkspx(false);
}

function mapchange(a, b)
{
  document.getElementById(b).value = document.getElementById(a).value;
  mkspx(true);
}

function init()
{
  //viewmat = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ];
  //viewmat = [ [0, 1, 0], [0, 0, 1], [1, 0, 0] ];
  viewmat = [[-0.6211501273837585,0.7805304954451938,-0.07031831149299674],
             [-0.4395286768639746,-0.2726761673313212,0.8558400843520434],
             [ 0.6488351573900389,0.5625120918252239,0.512438371987373]];
  installmouse("animation-box", "scaleinput");
  mkspx(true);
}

init();

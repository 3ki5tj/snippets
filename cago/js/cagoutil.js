


"use strict";



var PI = 3.141592653589793;



var aanames = [
  "GLY", "ALA", "VAL", "LEU", "ILE",
  "PRO", "SER", "THR", "CYS", "MET",
  "ASN", "GLN", "ASP", "GLU",
  "LYS", "ARG", "HIS",
  "PHE", "TYR", "TRP"];

var aaletters = [
  "G", "A", "V", "L", "I",
  "P", "S", "T", "C", "M",
  "N", "Q", "D", "E",
  "K", "R", "H",
  "F", "Y", "W"];

var aacolors = [
  "#404040", "#606060", "#8080a0", "#80a0a0", "#80a080",
  "#a08080", "#a080a0", "#a0a080", "#a0a000", "#808000",
  "#a040a0", "#e080e0", "#a00000", "#e00000",
  "#0000a0", "#0000e0", "#8080e0",
  "#00a0a0", "#00a000", "#804000"];

var aaradii = [
  0.90, 1.00, 1.10, 1.20, 1.20,
  1.10, 1.20, 1.20, 1.20, 1.40,
  1.20, 1.30, 1.25, 1.35,
  1.35, 1.40, 1.40,
  1.40, 1.40, 1.50];

/* residue name to integer of amino acid */
function res2iaa(res)
{
  var i;

  for ( i = 0; i < aanames.length; i++ ) {
    if ( res === aanames[i] ) {
      return i;
    }
  }
  return -1;
}




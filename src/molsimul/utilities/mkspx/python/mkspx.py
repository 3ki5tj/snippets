#!/usr/bin/env python3

''' build an extended chain (beta-strand) conformation
    adapted from the C program of the same name

  Copyright (C) 2010-2023  Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  make a spiral-like configuration
  use -h for more information about usage

  TODO:
  *  more advanced rotamer specification
  *  collision detection  '''

import sys, os, subprocess, getopt, shutil, re, random, glob
from math import *
import rv3



def usage():
  ''' print usage and die '''

  print("""  Usage:\n
    %s [OPTIONS] input""" % sys.argv[0])

  print("""
  Generate a PDB from a given sequence or PDB
  Copyright (C) 2010-2023 Cheng Zhang

  The `input' can be a PDB file, a PDB code, or a sequence file (see below)

  OPTIONS
  -------

    -o, --output=   the output PDB file, in which the coordinates assume
                    an extended helical spiral
    -R, --rotate=   rotation angle in degrees in the horizontal x-y plane
    -W, --swing=    the swinging angle in degrees between successive
                    peptide planes. If `rotate' is zero, the projection of
                    backbone trace on x-y plane is a zig-zag line around
                    the straight line, deflected only at alpha-carbons (CA)
                    The angles of segments at the joints is `swing'
    -S, --rise=     angle in degrees of rising in the vertical z-axis
    -T, --ter=      terminals caps, can be "N" (ACE), "C" (NH2) or "NC"
    --O0            set the N-terminal at the origin
    --ver=          specify GROMACS version, e.g., "4.5"
    -v, --verbose=  be verbose
    -h, --help      help


  Format of the sequence file

  The program can read a .seq file, which specifies the 3- or 4-letter amino acid sequence,
  such as NMET ALA LEU ....  The separators can be space, comma, etc.
  The entire SEQRES section of a PDB file would also work
  There can be an optional first line in the sequence file, which looks like
    # rotate swing rise letter-seq
  where,
    - all angles are in degrees
    - rotate is the angle of helical rotation in the x-y plane
    - swing is the swinging angle between successive CO-NH peptide planes
    - rise is the angle of the spiral rising along the z-axis
    - letter-seq is the one-letter amino acid sequence (no spaces)
      which is used to verify the sequence

  After the sequence, optional lines can be added to specify rotamer positions, e.g.
    # 93 G3
  specifies that C-gamma of the 93rd residue will placed at gamma 3 position
  G1 opposes CO; G2 opposes N; G3 is the remaining direction
  """)

  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvo:R:W:S:",
        [ "help", "verbose=", "ver=", "version=",
          "rotate=", "swing=", "rise=", "ter=",
          "O0", "ver=",
          "output=", ])
  except getopt.GetoptError as err:
    print(str(err))
    usage()

  gmxver = None
  start0 = False
  verbose = 0
  sver = fninp = fnout = None
  rotang = swgang = risang = None
  ter = ""
  for o, a in opts:
    if o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-o", "--output",):
      fnout = a
    elif o in ("-R", "--rotate",):
      rotang = float(a) * pi/180
    elif o in ("-W", "--swing",):
      swgang = float(a) * pi/180
    elif o in ("-S", "--rise",):
      risang = float(a) * pi/180
    elif o in ("-T", "--ter"):
      ter = a
    elif o in ("--O0",):
      start0 = True
    elif o in ("--ver", "--version",):
      sver = a
    elif o in ("-h", "--help"):
      usage()

  if sver:  # GROMACS version
    if "." in sver:
      sver = sver.replace(".", "0")
    if len(sver) < 5: sver += "0" * (5 - len(sver))
    gmxver = int(sver)

  if len(args) > 0:
    fninp = args[0]
  else:
    print("need an input")
    usage()

  print("input", fninp, " output", fnout, " GROMACS", gmxver, " verbose", verbose, end=None)
  print(" angles(rad) ", [rotang, swgang, risang], " ter ", ter)
  return fninp, fnout, [rotang, swgang, risang], ter, start0, gmxver, verbose



def main():
  fninp, fnout, angs, ter, start0, gmxver, verbose = doargs()

  if not os.path.exists(fninp) and re.match("[a-zA-Z0-9]{4}$", fninp):
    fninp += ".pdb"
    if os.system("wget http://www.rcsb.org/pdb/files/" + fninp) != 0:
      print("cannot download", fninp)
      raise Exception

  if fninp.lower().endswith(".pdb"):
    seq = pdb2seq(fninp)
  else:
    seq = readseq(fninp, angs)

  src = mkpdb(seq, angs, ter, start0, gmxver)

  if fnout == None: fnout = "out.pdb"
  print("writing", fnout)
  open(fnout, "w").write(src)



# amino acid letters
aaletter = {
  "ALA" : "A",
  "ARG" : "R",
  "ASN" : "N",
  "ASP" : "D",
  "CYS" : "C",
  "GLU" : "E",
  "GLN" : "Q",
  "GLY" : "G",
  "HIS" : "H",
  "ILE" : "I",
  "LEU" : "L",
  "LYS" : "K",
  "MET" : "M",
  "PHE" : "F",
  "PRO" : "P",
  "SER" : "S",
  "THR" : "T",
  "TRP" : "W",
  "TYR" : "Y",
  "VAL" : "V",
  }



def readseq(fn, angs):
  ''' read sequence file, return a sequence of [amino-acid, rotamer] '''

  s = open(fn).readlines()

  # tag information line
  i0 = 0
  if s[0].startswith("#"):
    ls = s[0][1:].strip().split()
    if len(ls) >= 3:
      # read parameters for the angles
      for j in range(3):
        if angs[j] == None:
          angs[j] = float(ls[j]) * pi/180
    i0 += 1

  # read the residue names
  i = i0
  seq = []
  while i < len(s):
    if s[i].startswith("#"): break
    ls = s[i].strip().split()
    # filter noises
    for a in ls:
      aa = a.strip()
      if aa == "SEQRES" or not aa.isalpha():
        continue

      # remove the terminal label
      if len(aa) >= 4 and aa[0] in "NC" and not aa.startswith("CY"):
        aa = aa[1:]

      # for a list of tuples (residue, rotamer)
      seq += [ [aa, 0], ]  # `0' means default
    i += 1

  # read rotamer positions, each line looks like: `# 93 G3'
  # where
  # `93' is the residue index starting from 1
  # `G3' is the rotamer position for the gamma carbon
  # G1 opposes CO; G2 opposes N; G3 is the remaining direction
  while i < len(s) and s[i].startswith("#"):
    ln = s[i].strip()
    m = re.match("^#\s*([0-9]+)\s*G\s*([1-3])\s*", ln)
    if not m:
      print("Warning: line %s of %s is corrupted: `%s'" % (i, fn, ln.strip()))
    sid = int( m.group(1) ) - 1
    rot = int( m.group(2) )
    if sid < 0 or sid >= len(seq):
      print("index %d out of range %s, %s" % (sid, len(seq), ln.strip()))
    seq[sid] = [ seq[sid][0], rot ]
    i += 1

  return seq



def pdb2seq(fnpdb):
  ''' gather sequence information from pdb '''

  # 1. determine the sequence from the PDB
  seq = []
  s = []
  # since the PDB file can be large, we only read one line at a time
  for ln in open(fnpdb).readlines():
    if ln.startswith("ENDMDL") or ln.startswith("TER"):
      break # stop after the first model
    if not ln.startswith("ATOM"):
      continue
    # save lines of atom coordinates
    s += [ ln, ]
    if not re.search("^ATOM\s*[0-9]+\s*CA ", ln):
      continue
    res = ln[17:21].strip()
    resid = ln[22:26].strip() # keep it as a string
    seq += [ [res, resid], ] # use resid as a temporary tag

  # 2. determine the rotamers of the CG atom
  for i in range(len(seq)):
    res, resid = seq[i]

    # collect lines with the same `resid'
    si = [ln for ln in s if ln.startswith("ATOM") and ln[22:26].strip() == resid]

    # collect coordinates
    xn = xc = xca = xcb = xcg = None
    for ln in si:
      atnm = ln[12:16].strip()
      # get the coordinates
      r = [ ln[30:38], ln[38:46], ln[46:54] ]
      r = [ float(x.strip()) for x in r ]
      if atnm == "CA": xca = r
      elif atnm == "N": xn = r
      elif atnm == "C": xc = r
      elif atnm == "CB": xcb = r
      elif atnm == "CG": xcg = r

    seq[i] = (res, 0)
    if xn and xc and xca and xcb and xcg:
      # compute the CG rotamer
      dir_can  = rv3.normalize(rv3.diff(xn, xca))
      dir_cac  = rv3.normalize(rv3.diff(xc, xca))
      dir_cgcb = rv3.normalize(rv3.diff(xcb, xcg))
      dot1 = rv3.dot(dir_cac, dir_cgcb)
      dot2 = rv3.dot(dir_can, dir_cgcb)
      # G1 position, opposite to C in the C=O group
      # G2 position, opposite to N in the N-H group
      # G3 position, p, q and u are perpendicular to w
      if dot1 < 0 and dot2 < 0: rotamer = 3
      elif dot1 > dot2: rotamer = 1
      else: rotamer = 2
      #print(i, resid, xn, xc, xca, xcb, xcg)
      seq[i] = (res, rotamer)
  return seq



def mkatom(atomid, atomname, resname, resid, x, ele = ""):
  if ele:
    if ele[0] != atomname[0]:
      print("element name %s disagrees with atom name %s" % (ele, atomname))
      raise Exception
  else:
    ele = atomname[0]
  return "ATOM   %4d  %-3s %-4s %4d    %8.3f%8.3f%8.3f  1.00  1.00          %s\n" % (
    atomid, atomname, resname, resid, x[0], x[1], x[2], ele)



def getminmax(r):
  ''' compute the system dimension '''

  xmin = [1000,] * 3
  xmax = rv3.neg(xmin)
  for d in range(3):
    ls = [ r[k][d] for k in range(len(r)) ]
    xmin[d] = min(ls)
    xmax[d] = max(ls)
  return xmin, xmax



def mkpdb(seq, angs, ter, start0, gmxver):
  ''' write a .pdb file with the given sequence information '''

  nres = len(seq)
  if not ter: ter = ""  # in case `ter == None'
  if "N" in ter: # add the N-terminal cap
    seq = [ ["ACE", None], ] + seq
  else: # add an empty N-terminal cap
    seq = [ [None, None], ] + seq

  D2R = pi/180.0
  R2D = 1.0/D2R

  # standard bond lengths in angstroms
  B_CAC     = 1.53
  B_CN_PEP  = 1.33
  B_NCA     = 1.46
  B_CC      = 1.54
  B_CN      = 1.48 # peptide bond
  B_CO      = 1.22 # carbonyl
  B_CS      = 1.81
  BH_OHOH   = 2.80
  BH_NHOH   = 2.90
  BH_OHOC   = 2.80
  B_CC_RING = 1.40

  # cosine and sine values
  c12, s12 = cos(12 * D2R), sin(12 * D2R)
  c30, s30 = cos(30 * D2R), sin(30 * D2R)
  c36, s36 = cos(36 * D2R), sin(36 * D2R)
  c72, s72 = cos(72 * D2R), sin(72 * D2R)

  r''' compute angles of the spiral
  * A peptide plane is defined as the plane of CA-C-N-CA
  * Below, all peptide planes are parallel to the vertical z-axis.
    Thus, the projection of each CA-C-N-CA on the horizontal x-y plane
    is a straight line segment.
  * The angle of rotation in the horizontal x-y plane is `rotang'
    More precisely, it is the angle of CA--CA--CA of three CA
  * Note, successive peptide planes switch directions (side view)
           C   CA  ...  CA   N
          / \ /          \ / \
        CA   N            C   CA
    Viewed from the top, the four atoms CA-C-N-CA are still aligned
    on the same straight line
  * If both `rotang' and `swgang' = 0, then the projection of the whole
    chain on the horizontal x-y plane is a straight line:
      CA--CA--CA--CA--CA--
    But if `rotang = 0' and `swgang != 0', then the projection on to the
    horizontal x-y plane like this:
      CA    CA    CA
        \  /  \  /  \
         CA    CA
    with `swgang' being the outer angle of CA--CA--CA.
    Here, the C and N atoms are omitted.
  * The use of both `rotang' and `swgang' is to make sure that
    the angle of N-CA-C is roughly 109 28'
  * `risang' is vertical tilting angle '''
  rotang, swgang, risang = angs
  # default values
  if rotang == None and swgang == None and risang == None:
    # the resulting size:
    # x-y width = sqrt(nres) * 3.5 + 14.5
    rotang = 110.0 / sqrt(nres) * D2R
    swgang = 50.0 * D2R
    # the formula makes height = 0.707 * width
    risang = (38.0 / sqrt(nres) + 81.0 / nres) * D2R
    #risang = rotang * 0.37 # z height  = sqrt(nres) * 2.7 + 3.3
  else:
    # partially missing
    if rotang == None: rotang = 10 * D2R
    if swgang == None: swgang = 50 * D2R
    if risang == None: risang = 9 * D2R

  ''' compute the average angle between N-CA with the horizontal plane
      The calculation uses the approximation of `rotang = 0'
      it assumes that the angle N-CA-C is 109.28
      the bond lengths are irrelevant. Let a = swgang, b = risang,
      u and v = angles of N-CA and CA-C with the horizontal plane
      so u = ang + b, v = ang - b
      For convenience, assuming the distances |N-CA| = |CA-C| = 1
        (if not true, relocate N and C at the respective line segments),
      we want |N-C| = sqrt(8/3) to make cos(N-CA-C) = -1/3. But
      |N-C|^2 = [cos^2 u + cos^2 v + 2 cos u cos v cos a] + (sin u - sin v)^2
              = 2 + cos(u+v)(1 + cos a) + cos(u-v)(cos a - 1)
      By u + v = 2 ang, and u - v = 2 b, `ang' can be solved as
        sin^2 ang = [sin^2 b cos a + sin^2 b - 1/3]/(1 + cos a) '''
  cr, sr = cos(risang), sin(risang)
  q = cos(swgang) * cr**2 + sr**2 - 1./3
  q /= 1 + cos(swgang)
  if q < 0:
    print("turning angle %s too large" % (swgang * R2D))
    raise Exception
  vang = asin( sqrt(q) )

  thp = vang + risang # for even-index residues
  c1p, s1p = cos(thp), sin(thp) # for CA-A or N-CA
  c2p, s2p = cos(thp - pi/3), sin(thp - pi/3) # for C-N

  thm = risang - vang # for odd-index residues
  c1m, s1m = cos(thm), sin(thm) # for CA-A or N-CA
  c2m, s2m = cos(thm + pi/3), sin(thm + pi/3) # for C-N

  phi = .5 * acos(-1./3)
  c3, s3 = cos(phi), sin(phi)

  print("angles", rv3.vround([rotang, swgang, risang], 5), "res", nres, end=None)
  print("theta+:", round(thp, 3), "theta-:", round(thm, 3))

  resid = 0
  n = len(seq)
  xyang = 0 # angle in the x-y plane
  os = [0, 0, 0] # current position
  atomls = [] # entries
  dir_nca = [0, 0, 0]
  xn = [0, 0, 0]

  for i in range(0, n):
    sgn = -(i % 2) * 2 + 1
    if sgn > 0: # even-index residues
      c1, s1, c2, s2 = c1p, s1p, c2p, s2p
    else: # odd-index residues
      c1, s1, c2, s2 = c1m, s1m, c2m, s2m

    # build the coordinates of the backbone atoms C, CA and O
    xyang += rotang
    cp, sp = cos(xyang), sin(xyang)

    # the peptide plane is vertical: parallel to the z-axis
    xca = os[:]

    # position of C
    dir_cac = [c1 * cp, c1 * sp, s1] # dir_cac is from CA to C
    xc = rv3.sadd(xca, dir_cac, B_CAC)

    # position of O
    xo = rv3.sadd(xc, [-sr * cp, -sr * sp, cr], B_CO * sgn)

    print(i, xyang, xc, xo)
    if i > 0:
      # position of CB
      # p is pointing to the side of CA opposite to N and C
      # q is perpendicular q
      # w is the position of CB that makes residue left-handed
      # w is the direction CA->CB
      p = rv3.normalize(rv3.diff(dir_nca, dir_cac))
      q = rv3.normalize(rv3.cross(dir_cac, dir_nca))
      w = rv3.normalize(rv3.lincomb2(p, q, c3, s3)) # CA->CB
      xcb = rv3.sadd(xca, w, B_CC)

      # compute the three CG positions
      # G1 position, opposite to C in the C=O group
      # viewed along the direction of CB-CA
      u = rv3.normalize( rv3.diff(xca, xc) )
      xog1 = rv3.sadd(xcb, u, B_CO)
      xcg1 = rv3.sadd(xcb, u, B_CC)
      # G2 position, opposite to N in the N-H group
      v = rv3.normalize( rv3.diff(xca, xn) )
      xog2 = rv3.sadd(xcb, v, B_CO)
      xcg2 = rv3.sadd(xcb, v, B_CC)
      # G3 position, p, q and u are perpendicular to w
      p = rv3.perpen(u, w)
      q = rv3.perpen(v, w)
      u = rv3.add(p, q)
      # v is the direction of CB->CG3
      v = rv3.normalize( rv3.lincomb2(u, w, -1., 1./3) )
      xog3 = rv3.sadd(xcb, v, B_CO)
      xcg3 = rv3.sadd(xcb, v, B_CC)
    else: # CB and CG atoms for the 0th residue, ACE or None
      xcb = None
      xog1 = xog2 = xog3 = xcg1 = xcg2 = xcg3 = None

    xnp = xn[:]
    xn = rv3.sadd(xc, [c2 * cp, c2 * sp, s2], B_CN_PEP)
    dir_nca = dir_cac[:]

    # set the origin for the next CA
    os = rv3.sadd(xn, dir_nca, B_NCA)

    xyang += swgang * sgn  # flip around the circle
    cp, sp = cos(xyang), sin(xyang)

    resnm = seq[i][0]

    if i > 0:
      atomls.append( ["N", resid, xnp] )
      atomls.append( ["CA", resid, xca] )
    elif resnm != None: # ACE as the 0th residue
      atomls.append( ["CH3", resid, xca] )

    if i > 0 or resnm != None: # we may skip the N-terminal ACE
      atomls.append( ["C", resid, xc] )
      if i == n - 1 and not "C" in ter: # no C-terminal cap
        atomls.append( ["OC1", resid, xo] )
        u = rv3.normalize( rv3.diff(xc, xo) )
        v = rv3.normalize( rv3.diff(xc, xca) )
        p = rv3.sadd(xc, rv3.add(u, v), B_CO)
        atomls.append( ["OC2", resid, p] )
      else:
        atomls.append( ["O", resid, xo] )

    if not resnm or resnm in ("GLY", "ACE",): # no CB
      resid += 1
      continue

    # put CB
    atomls.append( ["CB", resid, xcb] )

    # if seq[i][1] is not specified (0), then return the default `gd'
    # otherwise, return one of g1, g2 or g3
    GCHOOSE = lambda gd, g1, g2, g3: (gd, g1, g2, g3)[ seq[i][1] ]
    # put side chains
    if not resnm:
      pass
    elif resnm == "ALA":
      pass
    elif resnm == "SER":
      atomls.append( ["OG", resid, GCHOOSE(xog2, xog1, xog2, xog3)] )

    elif resnm.startswith("CY"):
      atomls.append( ["SG", resid, GCHOOSE(xog2, xog1, xog2, xog3)] )

    elif resnm == "VAL":
      atomls.append( ["CG1", resid, xcg1] )
      atomls.append( ["CG2", resid, xcg2] )

    elif resnm == "LEU":
      xg = GCHOOSE(xcg1, xcg1, xcg2, xcg3)
      q = rv3.diff(xcb, xca)
      xd1 = rv3.add(xg, q)
      u = rv3.diff(xg, xcb)
      v = rv3.sadd(q, u, -1./3)
      xdc = rv3.sadd(xg, u, 1./3) # center of xd1, xd2, xd3
      # reversely extend away from v, center of xd2 and xd3
      xd23 = rv3.sadd(xdc, v, -1./2)
      w = rv3.normalize(rv3.cross(q, u))
      xd2 = rv3.sadd(xd23, w, sqrt(2./3) * B_CC)

      atomls.append( ["CG",  resid, xg] )
      atomls.append( ["CD1", resid, xd1] )
      atomls.append( ["CD2", resid, xd2] )

    elif resnm == "ILE":
      xg1, xg2 = xcg1, xcg2
      if seq[i][1] == 3: xg1, xg2 = xcg3, xcg1
      xd = rv3.add(xg1, rv3.diff(xcb, xca))
      atomls.append( ["CG1", resid, xg1] )
      atomls.append( ["CD",  resid, xd] )
      atomls.append( ["CG2", resid, xg2] )

    elif resnm == "TRP":
      xg = GCHOOSE(xcg1, xcg1, xcg2, xcg3)
      v = rv3.normalize( rv3.diff(xg, xcb) )
      u = rv3.normalize( rv3.cross(v, rv3.diff(xca, xcb) ) )
      if xg == xcg2: u = rv3.neg(u)
      x8 = [[0,0,0],] * 8
      x8[0] = rv3.sadd(xg,    rv3.lincomb2(u, v,  c36,  s36), B_CC_RING) # CD1
      x8[1] = rv3.sadd(xg,    rv3.lincomb2(u, v, -c36,  s36), B_CC_RING) # CD2
      x8[2] = rv3.sadd(x8[0], rv3.lincomb2(u, v, -c72,  s72), B_CC_RING) # NE1
      x8[3] = rv3.sadd(x8[1], rv3.lincomb2(u, v,  c72,  s72), B_CC_RING) # CE2
      x8[4] = rv3.sadd(x8[1], rv3.lincomb2(u, v, -c12, -s12), B_CC_RING) # CE3
      x8[5] = rv3.sadd(x8[3], rv3.lincomb2(u, v, -c72,  s72), B_CC_RING) # CZ2
      x8[7] = rv3.sadd(x8[5], rv3.lincomb2(u, v, -c12, -s12), B_CC_RING) # CH2
      x8[6] = rv3.sadd(x8[7], rv3.lincomb2(u, v, -c72, -s72), B_CC_RING) # CZ3
      at8 = ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
      atomls.append( ["CG", resid, xg] )
      for k in range(8):
        atomls.append( [at8[k], resid, x8[k]] )

    elif resnm.startswith("HI"):
      xg = GCHOOSE(xcg2, xcg1, xcg2, xcg3)
      v = rv3.normalize( rv3.diff(xg, xcb) )
      w = rv3.normalize( rv3.cross( rv3.diff(xcb, xca), v) )
      xd1 = rv3.add(xg,  rv3.lincomb2(v, w, s36 * B_CC_RING,  c36 * B_CC_RING))
      xe1 = rv3.add(xd1, rv3.lincomb2(v, w, s72 * B_CC_RING, -c72 * B_CC_RING))
      xd2 = rv3.add(xg,  rv3.lincomb2(v, w, s36 * B_CC_RING, -c36 * B_CC_RING))
      xe2 = rv3.add(xd2, rv3.lincomb2(v, w, s72 * B_CC_RING,  c72 * B_CC_RING))
      atomls.append( ["CG",  resid, xg] )
      atomls.append( ["ND1", resid, xd1] )
      atomls.append( ["CE1", resid, xe1] )
      atomls.append( ["CD2", resid, xd2] )
      atomls.append( ["NE2", resid, xe2] )

    elif resnm in ("PHE", "TYR"):
      xg = GCHOOSE(xcg2, xcg1, xcg2, xcg3)
      v = rv3.normalize( rv3.diff(xg, xcb) )
      w = rv3.normalize( rv3.cross(rv3.diff(xcb, xca), v) )
      xd1 = rv3.add(xg, rv3.lincomb2(v, w, s30 * B_CC_RING, c30 * B_CC_RING))
      xe1 = rv3.sadd(xd1, v, B_CC_RING)
      xz = rv3.sadd(xg, v, 2 * B_CC_RING)
      xh = rv3.sadd(xz, v, B_CO)
      xd2 = rv3.add(xg, rv3.lincomb2(v, w, s30 * B_CC_RING, -c30 * B_CC_RING))
      xe2 = rv3.sadd(xd2, v, B_CC_RING)

      atomls.append( ["CG",  resid, xg] )
      atomls.append( ["CD1", resid, xd1] )
      atomls.append( ["CE1", resid, xe1] )
      atomls.append( ["CZ",  resid, xz] )
      if resnm == "TYR":
        atomls.append( ["OH", resid, xh] )
      atomls.append( ["CD2", resid, xd2] )
      atomls.append( ["CE2", resid, xe2] )

    elif resnm == "THR":
      # make sure CB is right-handed
      atomls.append( ["OG1", resid, xog1] )
      atomls.append( ["CG2", resid, xcg2] )

    elif resnm in ("GLU", "GLN",):
      xg = GCHOOSE(xcg2, xcg1, xcg2, xcg3)
      v = rv3.diff(xcb, xca)
      xcd = rv3.add(xg, v)
      v = rv3.normalize(v)
      u = rv3.normalize( rv3.cross(v, rv3.diff(xcb, xg)) )
      xe1 = rv3.sadd(xcd, rv3.lincomb2(u, v, c30, s30), B_CO)
      xe2 = rv3.sadd(xcd, rv3.lincomb2(u, v, -c30, s30), B_CO)
      nme2 = "NE2"
      if resnm == "GLU": nme2 = "OE2"
      atomls.append( ["CG", resid, xg] )
      atomls.append( ["CD", resid, xcd] )
      atomls.append( ["OE1", resid, xe1] )
      atomls.append( [nme2, resid, xe2] )

    elif resnm in ("ASP", "ASN",):
      xg = GCHOOSE(xcg1, xcg1, xcg2, xcg3)
      v = rv3.normalize(rv3.diff(xg, xcb))
      u = rv3.perpen(rv3.diff(xca, xcb), v)
      xd1 = rv3.sadd(xg, rv3.lincomb2(u, v,  c30, s30), B_CO)
      xd2 = rv3.sadd(xg, rv3.lincomb2(u, v, -c30, s30), B_CO)
      nmd2 = "ND2"
      if resnm == "ASP": nmd2 = "OD2"
      atomls.append( ["CG",  resid, xg] )
      atomls.append( ["OD1", resid, xd1] )
      atomls.append( [nmd2,  resid, xd2] )

    elif resnm == "LYS":
      xg = GCHOOSE(xcg2, xcg1, xcg2, xcg3)
      xd = rv3.add(xg, rv3.diff(xcb, xca))
      xe = rv3.add(xd, rv3.diff(xg, xcb))
      xz = rv3.sadd(xe, rv3.normalize(rv3.diff(xcb, xca)), B_CN)
      atomls.append( ["CG", resid, xg] )
      atomls.append( ["CD", resid, xd] )
      atomls.append( ["CE", resid, xe] )
      atomls.append( ["NZ", resid, xz] )

    elif resnm == "ARG":
      xg = GCHOOSE(xcg2, xcg1, xcg2, xcg3)
      q = rv3.diff(xcb, xca)
      xd = rv3.add(xg, q)
      v = rv3.normalize(rv3.diff(xg, xcb))
      xe = rv3.sadd(xd, v, B_CN)
      w = rv3.normalize( rv3.cross(v, q) )
      xz = rv3.add(xe, rv3.lincomb2(w, v, c30 * B_CN, s30 * B_CN))
      xh1 = rv3.add(xz, rv3.lincomb2(w, v, c30 * B_CN, -s30 * B_CN))
      xh2 = rv3.sadd(xz, v, B_CN)
      atomls.append( ["CG", resid, xg] )
      atomls.append( ["CD", resid, xd] )
      atomls.append( ["NE", resid, xe] )
      atomls.append( ["CZ", resid, xz] )
      atomls.append( ["NH1", resid, xh1] )
      atomls.append( ["NH2", resid, xh2] )

    elif resnm == "MET":
      xg = GCHOOSE(xcg1, xcg1, xcg2, xcg3)
      xd = rv3.sadd(xg, rv3.normalize(rv3.diff(xcb, xca)), B_CS)
      xe = rv3.sadd(xd, rv3.normalize(rv3.diff(xg,  xcb)), B_CS)
      atomls.append( ["CG", resid, xg] )
      atomls.append( ["SD", resid, xd] )
      atomls.append( ["CE", resid, xe] )

    elif resnm == "PRO":
      u = rv3.normalize(rv3.diff(xcb, xnp))
      p = rv3.diff(xcb, xca)
      q = rv3.diff(xnp, xca)
      v = rv3.normalize(rv3.add(p, q))
      xd = rv3.add(xcb, rv3.lincomb2(u, v, -c72 * B_CC, s72 * B_CC))
      xg = rv3.add(xnp, rv3.lincomb2(u, v,  c72 * B_CC, s72 * B_CC))
      atomls.append( ["CD", resid, xd] )
      atomls.append( ["CG", resid, xg] )

    else:
      print("residue %d: %s is not supported" % (i, resnm))

    resid += 1

  if "C" in ter:
    atomls.append( ["N", resid, xn] )
    seq += [ ["NH2", 0], ]
    resid += 1

  r = [ a[2] for a in atomls ]
  xmin, xmax = getminmax(r)
  print(len(seq)-1, "residues, System dimension:", rv3.vround(rv3.diff(xmax, xmin), 2), " rising", (risang/rotang))

  # shift the coordinates to make them nonnegative
  #if not start0:
  #  # add 2-angstrom margin
  #  xshift = rv3.add(xmin, [-2, -2, -2])

  #  for k in range(len(atomls)): # shift coordinates
  #    atnm, resid, x = atomls[k]
  #    atomls[k] = [atnm, resid, rv3.diff(x, xshift)]

  # write a PDB file
  src = ""
  resid = 0
  offset = 0
  if "N" in ter: offset = 1
  for k in range(len(atomls)):
    atnm, resid, x = atomls[k]
    resnm = seq[resid][0]
    src += mkatom(k + 1, atnm, resnm, resid + offset, x)
  src += "TER    %4d      %-4s %4d%52s\n" % (
      k + 2, resnm, resid + offset, " ")

  return src



if __name__ == "__main__":
  main()


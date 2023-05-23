#!/usr/bin/env python



''' make the three-interaction-site of a given PDB '''



import os, sys
from math import *
from vct import *



def R2D(x):
  return x * 180 / pi



def D2R(x):
  return x * pi / 180



class RNARes:
  def __init__(self):
    self.rP = None
    self.rS = None
    self.rB = None
    self.rwS = [] # corrinates for the sugar
    self.rwB = [] # corrdinates/weight for the base
    self.resname = "  X"
    self.vS = None



def a2w(a):
  if a == "C":
    w = 12
  elif a == "N":
    w = 14
  elif a == "H":
    w = 1
  elif a == "O":
    w = 16
  else:
    print "unknown atom ", a
    raise Error
  return w



def parsePDB(fnpdb):
  ''' parse PDB into residues
      see http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
      for the PDB format
  '''
  s = open(fnpdb).readlines()

  resarr = [int(ln[22:26]) for ln in s if ln.startswith("ATOM ")]
  resmin = min(resarr)
  resmax = max(resarr)

  nres = resmax - resmin + 1
  if "ideal" in fnpdb:
    res = [None] * (nres * 2)
  else:
    res = [None] * nres

  for i in range(len(res)):
    res[i] = RNARes()

  if "ideal" in fnpdb:
    i0 = 0
  else:
    i0 = -(resmin - 1)

  for ln in s:
    if ln.startswith("TER") and "ideal" in fnpdb:
      i0 += nres
    if not ln.startswith("ATOM  "):
      continue
    i = int(ln[22:26]) - 1 + i0
    res[i].resid = int(ln[22:26])  # the nominal residue id

    atomname = ln[12:16].strip()
    # a hack for HO'
    if ln[16] == "'": atomname += "'"

    x = [ float(ln[30:38]), float(ln[38:46]), float(ln[46:54]) ]
    if atomname == "P":
      res[i].rP = x
    elif atomname in ("C1'", "H1'",
                      "C2'", "H2'", "O2'", "HO2", "HO2'",
                      "C3'", "H3'",
                      "C4'", "O4'", "H4'",
                      "C5'", "1H5'", "2H5'") and "H" not in atomname:

      if atomname == "C1'":
        res[i].rC1p = x

      # note that we exclude hydrogen atoms in computing the com
      w = a2w( ln[13] )
      res[i].rwS += [(x, w), ]
    elif not atomname.endswith("'") and not atomname.endswith("P"):
      if atomname == "O2":
        res[i].rO2 = x
      elif atomname == "C2":
        res[i].rC2 = x
      elif atomname == "N1":
        res[i].rN1 = x
      elif atomname == "N3":
        res[i].rN3 = x

      w = a2w( ln[13] )
      res[i].rwB += [(x, w), ]
    res[i].resname = ln[17:20]

  # compute the center of mass of the sugar and base
  for i in range(len(res)):
    rS = [0, 0, 0]
    rwS = 0.0
    for x, w in res[i].rwS:
      rS[0] += x[0] * w
      rS[1] += x[1] * w
      rS[2] += x[2] * w
      rwS += w
    if rwS <= 0.0:
      print "residue %d has no sugar" % (i+1)
      continue
    else:
      rS[0] /= rwS
      rS[1] /= rwS
      rS[2] /= rwS
      res[i].rS = rS

    rB = [0, 0, 0]
    rwB = 0.0;
    for x, w in res[i].rwB:
      rB[0] += x[0] * w
      rB[1] += x[1] * w
      rB[2] += x[2] * w
      rwB += w
    if rwB <= 0.0:
      print "residue %d has no base" % (i+1)
      continue
    else:
      rB[0] /= rwB
      rB[1] /= rwB
      rB[2] /= rwB
      res[i].rB = rB
    print "%4d %s %6.1f %6.1f" % (i + 1, res[i].resname, rwS, rwB)
    #raw_input()

    resnm = res[i].resname.strip()
    print resnm, i
    if resnm == "C":
      res[i].rtip = res[i].rO2
      res[i].rtip2 = res[i].rN3
    elif resnm == "G":
      res[i].rtip = res[i].rC2
      res[i].rtip2 = res[i].rN1
    elif resnm == "A":
      res[i].rtip = res[i].rC2
      res[i].rtip2 = res[i].rN1
    elif resnm == "U" or resnm == "T":
      res[i].rtip = res[i].rO2
      res[i].rtip2 = res[i].rN3
    res[i].vS = vdiff(res[i].rtip, res[i].rC1p)
  return res



def writetisPDB(res, fn):
  ''' write residues to `fn` '''
  nres = len(res)
  atomid = 1
  pdbfmt = "ATOM  %5d %4s %-3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
  s = ""
  for i in range(nres):
    if not res[i].rP or not res[i].rS or not res[i].rB:
      continue
    resnm = res[i].resname
    resid = res[i].resid
    s += pdbfmt % (atomid, "P  ", resnm, resid, res[i].rP[0], res[i].rP[1], res[i].rP[2])
    atomid += 1
    s += pdbfmt % (atomid, "S  ", resnm, resid, res[i].rS[0], res[i].rS[1], res[i].rS[2])
    atomid += 1
    s += pdbfmt % (atomid, "B  ", resnm, resid, res[i].rB[0], res[i].rB[1], res[i].rB[2])
    atomid += 1
  print "writing %s" % fn
  open(fn, "w").write(s)



def helixparams(r1, r2, r3):
  '''
      r1^2 = 4 R^2 sin^2(phi/2)   +   a^2 = 2 R^2 (1 - sin(phi)  ) +   a^2
      r2^2 = 4 R^2 sin^2(phi*3/2) + 4 a^2 = 2 R^2 (1 - sin(3*phi)) + 4 a^2
      r3^2 = 4 R^2 sin^2(phi*5/2) + 9 a^2 = 2 R^2 (1 - sin(5*phi)) + 9 a^2
      so
        (r3^2 - 9 r1^2) / (r2^2 - 4 r1^2)
      = 6 - 4 sin^2(phi/2)
  '''
  rr = max( 6 - (r3*r3 - 9.0*r1*r1) / (r2*r2 - 4.0*r1*r1), 0 )
  phi = 2.0 * asin( 0.5 * sqrt(rr) )
  R = a = 0
  R = sqrt( 4*r1*r1 - r2*r2) / rr
  a = sqrt(r1*r1 - R*R*rr)
  return phi, R, a



def measure_ideal(res):
  nres = len(res)

  # coordinates of the phosphate
  print "Phosphate"
  for i in range(nres):
    x, y, z = res[i].rP
    r = sqrt(x*x + y*y)
    ang = atan2(y, x)
    print "%2d %8.3f %8.3f %8.3f %8.3f %8.3f(%8.3f)" % (i+1,
        x, y, z, sqrt(x*x+y*y), ang, ang*180/pi)

  # coordinates of the sugar
  print "Sugar"
  for i in range(nres):
    x, y, z = res[i].rS
    r = sqrt(x*x + y*y)
    ang = atan2(y, x)
    print "%2d %8.3f %8.3f %8.3f %8.3f %8.3f(%8.3f)" % (i+1,
        x, y, z, sqrt(x*x+y*y), ang, ang*180/pi)

  # coordinates of the base
  print "Base"
  for i in range(nres):
    x, y, z = res[i].rB
    r = sqrt(x*x + y*y)
    ang = fmod(atan2(y, x) - i*32.7*pi/180 - pi*200, 2*pi)
    print "%3s BASE %2d %8.3f %8.3f %8.3f %8.3f %8.3f(%8.3f)" % (
        res[i].resname, i+1,
        x, y, z-2.81*i, sqrt(x*x+y*y), ang, ang*180/pi)



def measure(res):
  nres = len(res)

  # bond lengths
  print "bond len. PS       SP       SB       PP       P2P      P3P"
  for i in range(nres):
    if not res[i].rP or not res[i].rS or not res[i].rB:
      continue
    dSB = vdist(res[i].rS, res[i].rB)
    dPS = vdist(res[i].rP, res[i].rS)
    if i < nres - 1 and res[i+1].rP:
      dSP = vdist(res[i].rS, res[i+1].rP)
      dPP = vdist(res[i].rP, res[i+1].rP)
    else:
      dSP = 0
      dPP = 0
    if i < nres - 2 and res[i+2].rP:
      dP2P = vdist(res[i].rP, res[i+2].rP)
    else:
      dP2P = 0
    if i < nres - 3 and res[i+3].rP:
      dP3P = vdist(res[i].rP, res[i+3].rP)
    else:
      dP3P = 0

    #if i < nres/2 - 3: # compute the helix parameters of P
    #  apt, rad, rise = helixparams(dPP, dP2P, dP3P)
    #  print "ang per turn %8.5f (%8.2f), radius %8.3f, rise %8.3f" % (apt, apt*180/pi, rad, rise)

    print  "%2d %s %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f" % (
        i + 1, res[i].resname, dPS, dSP, dSB, dPP, dP2P, dP3P)

  # bond angles
  print "bond ang. PSB                BSP                PSP                PPP"
  for i in range(nres):
    if not res[i].rP or not res[i].rS or not res[i].rB:
      continue
    aPSB = vang(res[i].rP, res[i].rS, res[i].rB)
    if i < nres - 1 and res[i+1].rP:
      aBSP = vang(res[i].rB, res[i].rS, res[i+1].rP)
      aPSP = vang(res[i].rP, res[i].rS, res[i+1].rP)
    else:
      aBSP = 0
      aPSP = 0

    if i < nres - 2 and res[i+1].rP and res[i+2].rP:
      aPPP = vang(res[i].rP, res[i+1].rP, res[i+2].rP)
    else:
      aPPP = 0

    print  "%2d %s %8.5f(%8.2f) %8.5f(%8.2f) %8.5f(%8.2f) %8.5f(%8.2f)" % (
        i + 1, res[i].resname,
        aPSB, aPSB*180/pi,
        aBSP, aBSP*180/pi,
        aPSP, aPSP*180/pi,
        aPPP, aPPP*180/pi)

  # stacking interaction
  print "stack     r0       phi1               phi2"
  for i in range(nres - 1):
    if (not res[i].rP or not res[i].rS or not res[i].rB
     or not res[i+1].rP or not res[i+1].rS or not res[i+1].rB):
      continue
    if i == nres/2 - 1 and nres == 22:
      continue
    dis = vdist(res[i].rB, res[i+1].rB)
    dih1 = vdih(res[i].rP, res[i].rS, res[i+1].rP, res[i+1].rS)
    if ( i < nres - 2 and i != nres/2 - 2
      and res[i+1].rP and res[i+1].rS and res[i+2].rP ):
      dih2 = vdih(res[i].rS, res[i+1].rP, res[i+1].rS, res[i+2].rP)
    else:
      dih2 = 0
    print "%2d %s-%s %8.5f %8.5f(%8.2f) %8.5f(%8.2f)" % (
        i+1, res[i].resname, res[i+1].resname.strip(),
        dis, dih1, dih1*180/pi, dih2, dih2*180/pi)



def measure_hbonds(res, fn):
  nres = len(res)

  if "ideal" in fn:
    ip = 1
  else:
    ip = 2
  dmax = 5
  d2max = 4
  dotmax = -0.9
  dotimax = -0.8

  print "Hydrogen bonds"
  print "   i-j    i-j    R0     TH10               TH20               PHI0               PHI10              PHI20"
  for i in range(nres):
    if not res[i].rB or not res[i].rS:
      continue
    for j in range(i + ip, nres):
      if not res[j].rB or not res[j].rS:
        continue
      dvi = res[i].vS # vdiff(res[i].rB, res[i].rS)
      dvj = res[j].vS # vdiff(res[j].rB, res[j].rS)
      dot = vdot(dvi, dvj) / sqrt( vsqr(dvi) * vsqr(dvj) )
      vij = vdiff(res[i].rtip, res[j].rtip)
      doti = vdot(vij, dvi) / sqrt( vsqr(vij) * vsqr(dvi) )
      dij = vdist(res[i].rtip, res[j].rtip)
      dij2 = vdist(res[i].rtip2, res[j].rtip2)
      if dij < dmax and dij2 < d2max and dot < dotmax and doti < dotimax:
        rbb = vdist(res[i].rB, res[j].rB)
        th1 = vang(res[i].rS, res[i].rB, res[j].rB)
        th2 = vang(res[i].rB, res[j].rB, res[j].rS)
        phi = vdih(res[i].rS, res[i].rB, res[j].rB, res[j].rS)
        phi1 = phi2 = 0
        if res[i].rP:
          phi1 = vdih(res[i].rP, res[i].rS, res[i].rB, res[j].rB)
        if res[j].rP:
          phi2 = vdih(res[i].rB, res[j].rB, res[j].rS, res[j].rP)
        print "%4d-%-4d %s-%-s %8.3f %8.5f(%8.3f) %8.5f(%8.3f) %8.5f(%8.3f) %8.5f(%8.3f) %8.5f(%8.3f)" % (
            res[i].resid, res[j].resid,
            res[i].resname.strip(), res[j].resname.strip(),
            rbb, th1, R2D(th1), th2, R2D(th2),
            phi, R2D(phi), phi1, R2D(phi1), phi1, R2D(phi2))
        #print "%4d-%-4d %s-%-s %8.3f %8.3f %8.3f %8.3f" % (
        #    res[i].resid, res[j].resid,
        #    res[i].resname.strip(), res[j].resname.strip(),
        #    dij, dij2, dot, doti)



def mktis(fnpdb):
  ''' make the three-interaction-site model '''
  res = parsePDB(fnpdb)

  if "ideal" in fnpdb:
    measure_ideal(res)
  measure(res)
  measure_hbonds(res, fnpdb)

  fx = os.path.splitext(fnpdb)
  fnout = fx[0] + ".tis" + fx[1]
  writetisPDB(res, fnout)



if __name__ == "__main__":
  fnpdb = "refs/ideal_arns.pdb"
  if len(sys.argv) > 1:
    fnpdb = sys.argv[1]
  mktis(fnpdb)


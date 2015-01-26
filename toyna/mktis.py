#!/usr/bin/env python



''' make the three-interaction-site of a given PDB '''



import os, sys
from math import *
from vct import *



class RNARes:
  def __init__(self):
    self.rwS = [None]*6 # corrinates for the sugar
    self.rwB = []  # corrdinates/weight for the base
    pass



def parsePDB(fnpdb):
  ''' parse PDB into residues
      see http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
      for the PDB format
  '''
  s = open(fnpdb).readlines()
  nres = max( int(ln[22:26]) for ln in s if ln.startswith("ATOM "))
  res = [None] * (nres * 2)
  for i in range(len(res)):
    res[i] = RNARes()

  i0 = 0
  for ln in s:
    if ln.startswith("TER"):
      i0 += nres
    if not ln.startswith("ATOM  "):
      continue
    i = int(ln[22:26]) - 1 + i0

    atomname = ln[12:16].strip()
    # a hack for HO'
    if ln[16] == "'": atomname += "'"

    x = [ float(ln[30:38]), float(ln[38:46]), float(ln[46:54]) ]
    if atomname == "P":
      res[i].rP = x
    elif atomname == "C1'":
      res[i].rwS[0] = (x, 12)
    elif atomname == "C2'":
      res[i].rwS[1] = (x, 12)
    elif atomname == "C3'":
      res[i].rwS[2] = (x, 12)
    elif atomname == "C4'":
      res[i].rwS[3] = (x, 12)
    elif atomname == "C5'":
      res[i].rwS[4] = (x, 12)
    elif atomname == "O4'":
      res[i].rwS[5] = (x, 16)
    elif not atomname.endswith("'") and not atomname.endswith("P"):
      a = ln[13]
      if a == "C":
        w = 12
      elif a == "N":
        w = 14
      elif a == "H":
        w = 1
      elif a == "O":
        w = 16
      else:
        print "unknown atom ", a, ln
        raise Error
      res[i].rwB += [(x, w),]
    res[i].resname = ln[17:20]

  # compute the center of mass of the sugar and base
  for i in range(len(res)):
    rS = [0, 0, 0]
    rw = 0.0
    for x, w in res[i].rwS:
      rS[0] += x[0] * w
      rS[1] += x[1] * w
      rS[2] += x[2] * w
      rw += w
    rS[0] /= rw
    rS[1] /= rw
    rS[2] /= rw
    res[i].rS = rS

    rB = [0, 0, 0]
    rw = 0.0;
    for x, w in res[i].rwB:
      rB[0] += x[0] * w
      rB[1] += x[1] * w
      rB[2] += x[2] * w
      rw += w
    rB[0] /= rw
    rB[1] /= rw
    rB[2] /= rw
    res[i].rB = rB
    #print i, rS, res[i].rwS
    #print i, rB, res[i].rwB
    #raw_input()
  return res



def writetisPDB(res, fn):
  ''' write residues to `fn` '''
  nres = len(res)
  atomid = 1
  pdbfmt = "ATOM  %5d %4s %-3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
  s = ""
  for i in range(nres):
    s += pdbfmt % (atomid, "P  ", res[i].resname, i + 1, res[i].rP[0], res[i].rP[1], res[i].rP[2])
    atomid += 1
    s += pdbfmt % (atomid, "C  ", res[i].resname, i + 1, res[i].rS[0], res[i].rS[1], res[i].rS[2])
    atomid += 1
    s += pdbfmt % (atomid, "N  ", res[i].resname, i + 1, res[i].rB[0], res[i].rB[1], res[i].rB[2])
    atomid += 1
  print "writing %s" % fn
  open(fn, "w").write(s)



def measure(res):
  nres = len(res)

  # bond lengths
  print "bond len. SB       PS       SP"
  for i in range(nres):
    dSB = vdist(res[i].rS, res[i].rB)
    dPS = vdist(res[i].rP, res[i].rS)
    if i < nres - 1:
      dSP = vdist(res[i].rS, res[i+1].rP)
    else:
      dSP = 0
    print  "%2d %s %8.5f %8.5f %8.5f" % (
        i + 1, res[i].resname, dPS, dSP, dSB)

  # bond angles
  print "bond ang. PSB                BSP                PSP"
  for i in range(nres):
    aPSB = vang(res[i].rP, res[i].rS, res[i].rB)
    if i < nres - 1:
      aBSP = vang(res[i].rB, res[i].rS, res[i+1].rP)
      aPSP = vang(res[i].rP, res[i].rS, res[i+1].rP)
    else:
      aBSP = 0
      aPSP = 0
    print  "%2d %s %8.5f(%8.2f) %8.5f(%8.2f) %8.5f(%8.2f)" % (
        i + 1, res[i].resname,
        aPSB, aPSB*180/pi,
        aBSP, aBSP*180/pi,
        aPSP, aPSP*180/pi)

  # stacking interaction
  print "stack     r0       phi1     phi2"
  for i in range(nres - 1):
    if i == nres/2 - 1: continue
    dis = vdist(res[i].rB, res[i+1].rB)
    dih1 = vdih(res[i].rP, res[i].rS, res[i+1].rP, res[i+1].rS)
    if i < nres - 2 and i != nres/2 - 2:
      dih2 = vdih(res[i].rS, res[i+1].rP, res[i+1].rS, res[i+2].rP)
    else:
      dih2 = 0
    print "%2d %s-%s %8.5f %8.5f %8.5f" % (
        i+1, res[i].resname, res[i+1].resname.strip(),
        dis, dih1, dih2)



def mktis(fnpdb):
  ''' make the three-interaction-site model '''
  res = parsePDB(fnpdb)

  measure(res)

  fx = os.path.splitext(fnpdb)
  fnout = fx[0] + ".tis" + fx[1]
  writetisPDB(res, fnout)



if __name__ == "__main__":
  if len(sys.argv) > 1:
    fnpdb = sys.argv[1]
  mktis(fnpdb)

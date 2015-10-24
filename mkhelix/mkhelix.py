#!/usr/bin/env python



'''
make a PDB file for an ideal alpha helix
'''


import os, sys
from math import *


# number of residues
nres = 10
fnout = "helix.pdb"
cter = 1


def doargs():
    pass



def mkxyz(i, r, z0, dz, th0, dth):
    ''' return the coordinates from helix parameters '''
    th = (th0 + i * dth) * pi / 180
    z = z0 + i * dz
    return [ r * cos(th), r * sin(th), z]


def mkpdbatom(pdbwriter, atomnm, resid, xyz):
    ''' write an PDB ATOM entry '''
    elem = atomnm.strip()[0]
    atomid = pdbwriter[1]
    pdbwriter[0] += "ATOM  %5d  %-3s ALA  %4d    %8.3f%8.3f%8.3f  1.00  1.00          %2s  \n" % (
        atomid, atomnm, resid, xyz[0], xyz[1], xyz[2], elem)



def mkhelix(nres, fn):
    ''' alpha helix '''

    pos = [None] * (nres + 3)

    dz = 1.46262
    dth = 100
    for i in range(nres + 2):
        # each amino acid is a dictionary
        aa = {}

        aa["xn"]  = mkxyz(i, 1.53605, -0.073, dz,  5.163, dth)
        aa["xhn"] = mkxyz(i, 1.52663, -1.041, dz, 14.478, dth)
        aa["xca"] = mkxyz(i, 2.27759,  0.812, dz, 30.711, dth)
        aa["xha"] = mkxyz(i, 2.99128,  1.311, dz, 18.007, dth)
        aa["xcb"] = mkxyz(i, 3.31126, -0.012, dz, 48.562, dth)
        aa["xc"]  = mkxyz(i, 1.65845,  1.865, dz, 59.869, dth)
        aa["xo"]  = mkxyz(i, 1.91447,  3.054, dz, 53.955, dth)

        pos[i] = aa

    pdbwriter = ["", 1]
    for i in range(1, nres + 1):
        aa = pos[i]
        mkpdbatom(pdbwriter, "N",  i, aa["xn"])
        mkpdbatom(pdbwriter, "NH", i, aa["xhn"])
        mkpdbatom(pdbwriter, "CA", i, aa["xca"])
        mkpdbatom(pdbwriter, "HA", i, aa["xha"])
        mkpdbatom(pdbwriter, "CB", i, aa["xcb"])
        mkpdbatom(pdbwriter, "C",  i, aa["xc"])
        mkpdbatom(pdbwriter, "O",  i, aa["xo"])

    print pdbwriter[0]



if __name__ == "__main__":
    doargs()
    mkhelix(nres, fnout)

#!/usr/bin/env python3



'''
make a PDB file for an ideal alpha helix
'''


import os, sys, getopt
from math import *


# number of residues
nres = 10
fnout = None
nter = 0  # add an ACE residue at the N-terminal
cter = 0  # add an NH2 residue at the C-terminal
verbose = 0



def help():
    ''' print usage and die '''
    print("%s [Options]" % sys.argv[0])
    print("  -n :             specify the number of residues")
    print("  -o :             specify the output file, if none, print(the output on screen")
    print("  --nter, -N:      add a N-terminal residue ACE")
    print("  --cter, -C:      add a C-terminal residue NH2")
    print("  --ter, -T:       add both terminal residues")
    sys.exit(1)


def do_args():
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "hn:o:TNcCv:",
             ["help", "output=", "ter", "nter", "cter", "verbose="])
    except getopt.GetoptError as err:
        # print(help information and exit:
        print(str(err)) # will print(something like "option -a not recognized"
        help()

    global nres, fnout, nter, cter, verbose

    for o, a in opts:
        if o in ("-n",):
            nres = int(a)
        elif o in ("-o", "--output"):
            fnout = a
        elif o in ("-N", "--nter"):
            nter = 1
        elif o in ("-C", "-c", "--cter"):
            cter = 1
        elif o in ("-T", "--ter"):
            cter = nter = 1
        elif o in ("-v",):
            verbose += 1
        elif o in ("--verbose",):
            verbose = int(a)
        elif o in ("-h", "--help"):
            help()
        else:
            print("unknown option '%s'" % o)
            help()



def mkxyz(i, r, z0, dz, th0, dth):
    ''' return the coordinates from helix parameters '''
    th = (th0 + i * dth) * pi / 180
    z = z0 + i * dz
    return [ r * cos(th), r * sin(th), z]



def mkpdbatom(pdbwriter, atomnm, resid, xyz, resname = "ALA"):
    ''' write an PDB ATOM entry '''
    elem = atomnm.strip()[0]
    atomid = pdbwriter[1]
    pdbwriter[0] += [ "ATOM  %5d  %-3s %3s  %4d    %8.3f%8.3f%8.3f  1.00  1.00          %2s  \n" % (
        atomid, atomnm, resname, resid, xyz[0], xyz[1], xyz[2], elem), ]



def make_helix(nres, fn):
    ''' alpha helix '''

    pos = [None] * (nres + 3)

    dz = 1.46262
    dth = 100
    for i in range(nres + 2):
        # each amino acid is a dictionary
        aa = {}

        aa["xn"]  = mkxyz(i, 1.53605, -0.073, dz,  5.163, dth)
        aa["xh"]  = mkxyz(i, 1.52663, -1.041, dz, 14.478, dth)
        aa["xca"] = mkxyz(i, 2.27759,  0.812, dz, 30.711, dth)
        aa["xha"] = mkxyz(i, 2.99128,  1.311, dz, 18.007, dth)
        aa["xcb"] = mkxyz(i, 3.31126, -0.012, dz, 48.562, dth)
        aa["xc"]  = mkxyz(i, 1.65845,  1.865, dz, 59.869, dth)
        aa["xo"]  = mkxyz(i, 1.91447,  3.054, dz, 53.955, dth)

        pos[i] = aa

    pdbwriter = [[], 1]

    resid = 1

    # add the N-terminal
    if nter:
        aa = pos[0]
        mkpdbatom(pdbwriter, "CH3", resid, aa["xca"], resname = "ACE")
        mkpdbatom(pdbwriter, "C",   resid, aa["xc"],  resname = "ACE")
        mkpdbatom(pdbwriter, "O",   resid, aa["xo"],  resname = "ACE")
        resid += 1

    # `i` is the internal residue index
    # `resid` is the output residue index
    for i in range(1, nres + 1):
        aa = pos[i]
        mkpdbatom(pdbwriter, "N",  resid, aa["xn"])
        if i > 1 or nter:
          mkpdbatom(pdbwriter, "H", resid, aa["xh"])
        mkpdbatom(pdbwriter, "CA", resid, aa["xca"])
        mkpdbatom(pdbwriter, "HA", resid, aa["xha"])
        mkpdbatom(pdbwriter, "CB", resid, aa["xcb"])
        mkpdbatom(pdbwriter, "C",  resid, aa["xc"])
        mkpdbatom(pdbwriter, "O",  resid, aa["xo"])
        resid += 1

    # add the C-terminal
    if cter:
        aa = pos[nres + 1]
        mkpdbatom(pdbwriter, "N",  nres + 1, aa["xn"], resname = "NH2")

    # add TER and END
    s = pdbwriter[0][-1]
    s = "TER   " + s[6:26] + " " * 54 + "\n"
    pdbwriter[0] += [s, "END" + " " * 77 + "\n"]

    # print(the output or save it to file
    if not fnout:
        print(''.join( pdbwriter[0] ))
    else:
        open(fnout, "w").writelines( pdbwriter[0] )



if __name__ == "__main__":
    do_args()
    make_helix(nres, fnout)

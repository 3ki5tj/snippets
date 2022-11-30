#!/usr/bin/env python

import os, sys, glob, re

fns = glob.glob("walk*.log")
fnout = "walk_Svsq.dat"
tval = 1000

dat = {}
qs = []
for fn in fns:
    m = re.match("walk_q([0-9]+e?[0-9]*)(_g.*)?", os.path.splitext(fn)[0])
    qv = float(m.group(1))
    for ln in open(fn).readlines():
        arr = ln.split()
        if int(arr[0]) == tval:
            dat[qv] = ln
            break
    else:
        print "cannot find time point %s for %s" % (tval, fn)
        continue
    qs += [qv,]

s = ""
for qv in sorted(qs):
    s += "\t".join( ["%g" % qv,] + dat[qv].split()[1:] ) + "\n"
print "saving data to %s" % fnout
open(fnout, "w").writelines(s)




#!/usr/bin/env python

''' remove line breaks for large numbers in IsingDOS*.txt '''

import sys, os, glob

def inttrim(fn):
    s = open(fn).read()
    s1 = s.replace("\\\n \n>    ", "")
    if s != s1:
        open(fn, "w").write(s1)

for fn in glob.glob("Ising*.txt"):
    inttrim(fn) 

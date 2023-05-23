#!/usr/bin/env python

from math import *


KB        = 0.0019872041
AVOGADRO  = 6.022140857e23
HBAR      = (1.054571800e-34*AVOGADRO*1e12/4184)
H = HBAR * (2*pi)

kT = KB*300
m0 = 12/418.4
b = 3.79
theta = 113*pi/180
vol = 1;

# 3D correction factor for Cartesian PCA on RMSD fit structure
n = 3
smass = m0 * n
detI = (m0*b*b)**3*2./9*(2-cos(theta))*(sin(theta)**2)
fac = 8 * pi**2 * vol * smass**1.5 * (2*pi*kT)**3 * detI**0.5 / H**6
print "D = 3", fac, log(fac)

# 2D correction factor for Cartesian PCA on RMSD fit structure
smass = m0 * n
detI = (m0*b*b)*2./3*(2-cos(theta))
fac = 2 * pi * vol * smass * (2*pi*kT)**1.5 * detI**0.5 / H**3
print "D = 2", fac, log(fac)

#!/usr/bin/env python3

""" verify the velocity inverse matrix """

import numpy as np
from math import pi, cos, sin

mA, mB, mC = 12.011, 1.008, 1.008
invmA, invmB, invmC = 1.0/mA, 1.0/mB, 1.0/mC

ang_hoh = 109.28*pi/180
ang_ohh = (pi - ang_hoh)*0.5
cos_hoh = cos(ang_hoh)
cos_ohh = cos(ang_ohh)
sin_ohh = sin(ang_ohh)

dist_oh = 1.0
dist_on = dist_oh * sin_ohh
dist_hn = dist_oh * cos_ohh
dist_hh = dist_hn * 2

dist_AB = dist_BA = dist_oh
dist_AC = dist_CA = dist_oh
dist_BC = dist_CB = dist_hh

cosA, cosB, cosC = cos_hoh, cos_ohh, cos_ohh

#mat = np.array([
#    [ma+mb, ma*cosB, mb*cosA],
#    [mc*cosB, mb+mc, mb*cosC],
#    [mc*cosA, ma*cosC, ma+mc],
#])
mat = np.array([
    [
        (invmB + invmC)*dist_BC*dist_BC,  # lambda_BC xBC
        invmC*dist_AC*dist_BC*cosC,
        invmB*dist_BA*dist_BC*cosB,
    ],
    [
        invmC*dist_BC*dist_CA*cosC,
        (invmC + invmA)*dist_CA*dist_CA,
        invmA*dist_BA*dist_CA*cosA,
    ],
    [
        invmB*dist_AB*dist_CB*cosB,
        invmA*dist_AB*dist_AC*cosA,
        (invmA + invmB)*dist_AB*dist_AB,   # lambda_AB  x_AB.x_AB
    ],
])
mat = -mat


print("Velocity matrix:", mat, sep="\n", end="\n\n")
invmat = np.linalg.inv(mat)
print("Inverse velocity matrix:", invmat, sep="\n", end="\n\n")


b = np.array([1., 1., 1.])
x = np.linalg.solve(mat, b)
#print(x)


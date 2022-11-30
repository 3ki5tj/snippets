

from math import *


def vadd(a, b):
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]



def vdiff(a, b):
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]



def vnadd(a, b):
  return [-a[0] - b[0], -a[1] - b[1], -a[2] - b[2]]



def vdot(a, b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]



def vsqr(a):
  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]



def vdist(a, b):
  return sqrt( vsqr( vdiff(a, b) ) )


def vsmul(a, s):
  return [a[0] * s, a[1] * s, a[2] * s]



def vcross(x, y):
  return [x[1]*y[2] - x[2]*y[1],
          x[2]*y[0] - x[0]*y[2],
          x[0]*y[1] - x[1]*y[0]]



def vang(a, b, c):
  ab = vdist(a, b)
  bc = vdist(b, c)
  ac = vdist(a, c)
  return acos((ab*ab + bc*bc - ac*ac)/(2*ab*bc))



def vdih(xi, xj, xk, xl):
  ''' return the dihedral angle of i-j-k-l,
      which is computed as the angle between the jik and jlk
      the dihedral angle is 0 in the cis configuration
      and pi in the trans configuration '''
  xij = vdiff(xi, xj)
  xkj = vdiff(xk, xj)
  xkl = vdiff(xk, xl)
  nxkj2 = vsqr(xkj)
  nxkj = sqrt(nxkj2)
  tol = nxkj2 * 1e-16

  m = vcross(xij, xkj)
  m2 = vsqr(m)
  n = vcross(xkj, xkl)
  n2 = vsqr(n)
  if m2 > tol and n2 > tol:
    cosphi = vdot(m, n)
    cosphi /= sqrt(m2 * n2)
    cosphi = max(min(cosphi, 1), -1)
  phi = acos(cosphi)
  if vdot(n, xij) < 0.0:
    phi = -phi
  return phi


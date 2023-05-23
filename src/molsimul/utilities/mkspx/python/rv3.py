#!/usr/bin/env python3

''' 3D real (float) vector '''

import math

def neg(v): return [-v[0], -v[1], -v[2]]

def lincomb2(a, b, u, v):
  return [ a[0]*u + b[0]*v, a[1]*u + b[1]*v, a[2]*u + b[2]*v ]

def sadd(a, b, v): return lincomb2(a, b, 1, v)

def add(a, b): return sadd(a, b, 1)

def diff(a, b): return sadd(a, b, -1)

def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def sqr(a): return dot(a, a)

def norm(a): return math.sqrt(sqr(a))

def normalize(a):
  s = norm(a)
  if s <= 0: s = 1
  return [a[0]/s, a[1]/s, a[2]/s]

def cross(a, b):
  return [ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] ]

def perpen(a, b):
  ''' perpendicular component of `a' to `b', normalized '''
  return normalize( sadd(a, b, -dot(a,b)/sqr(b)) )

def vmin(a, b):
  return [min(a[i], b[i]) for i in range(3)]

def vround(a, n = 0):
  return [round(x, n) for x in a]


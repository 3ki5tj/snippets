/* vector routines in two dimensions */



"use strict";



function vzero(x)
{
  x[0] = 0;
  x[1] = 0;
  return x;
}



function vcopy(x, y)
{
  x[0] = y[0];
  x[1] = y[1];
  return x;
}



function vinc(x, dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  return x;
}



function vdec(x, dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  return x;
}



function vsinc(x, dx, s)
{
  x[0] += dx[0] * s;
  x[1] += dx[1] * s;
  return x;
}



function vdiff(c, a, b)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  return c;
}



function vsmul(x, s)
{
  x[0] *= s;
  x[1] *= s;
  return x;
}



function vsqr(x)
{
  return x[0] * x[0] + x[1] * x[1];
}



function vdot(x, y)
{
  return x[0] * y[0] + x[1] * y[1];
}



function vcross(x, y)
{
  return x[0]*y[1] - x[1]*y[0];
}



function vpbc(v, l, invl)
{
  v[0] -= (Math.floor(v[0]*invl + 1000.5) - 1000.0) * l;
  v[1] -= (Math.floor(v[1]*invl + 1000.5) - 1000.0) * l;
  return v;
}



function vwrap(v, l)
{
  while ( v[0] < 0 ) v[0] += l; while ( v[0] > l ) v[0] -= l;
  while ( v[1] < 0 ) v[1] += l; while ( v[1] > l ) v[1] -= l;
}

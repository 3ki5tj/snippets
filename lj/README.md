# Simulation of Lennard-Jones fluid #

## Compilation ##

The code support both 2D and 3D versions.
Use `md.c` as an example.

For the 3D version
```
gcc md.c -lm
```

For the 2D version
```
gcc -DD=2 md.c -lm
```

## ljeos.h

The default call is
```
U = ljeos3d_get(rho, T, &P, &Fex, &muex);
```

This function computes the average energy (return value),
 pressure, Helmholtz free energy, and the chemical potential
 of a Lennard-Jones fluid at a given number density, rho and temperature, T.
The last two quantities exclude the corresponding ideal-gas
(kinetic energy) parts.  If we want to check the logarithmic
partition function,  then (Fex/T) is the quantity to go.

There are several alternative equation of states in the literature.
The code contains three equations, as shown below.
The last one is the default, which is mapped to
the above function `ljeos3d_get()`.
```
ljeos3d_MBWRJZG(rho, T, P, Fex, muex);
ljeos3d_MBWRKN(rho, T, P, Fex, muex);
ljeos3d_PVEhBH(rho, T, P, Fex, muex);
```

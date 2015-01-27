# Simulation of toy nucleic acids #

## Compilation ##

Type
```
make
```


## References ##

### Three-interaction-site (TIS) model ###

Coarse-grained model for predicting RNA folding thermodynamics.
Natalia A. Denesyuk and D. Thirumalai,
The Journal of Physical Chemistry, 117, 4901-4911 (2013).

### PDB format ###

http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

### Dielectric constant ###

Dielectric constant of water.
\[
\epsilon = 87.740 - 0.40008 t + 9.398\times 10^{-4} t^2 - 1.410 \times 10^{-8} t^3,
\]
where $t$ is the temperature in Celsius ($t = T -273.15\mathrm{K}$).

Source:
Dielectric Constant of Water from 0 to 100$^\circ$C
by C. G. Malmberg and A. A. Maryott
Journal of Research of the National Bureau of Standards
Vol. 56, No. 1, January 1956, Research Paper 2641.


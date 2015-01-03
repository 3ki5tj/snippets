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

File           | Description
---------------|--------------------
is2dos0.ma     | exact density of states (Mathematica)
is2exact.c     | exact density of states and free energy, energy, heat capacity (C)


# Exact results

The density of states can be computed in two ways.
`is2lndosnxm.dat` is the output of C code, `is2exact.c`
`is2lndosnxm.txt` is the output of the Mathematica script `is2dos0.ma`.
Note the difference in the extension `.dat` vs. `.txt`,
but the contents should be the same within numerical error.


# Simulation of the Ising model

## Compilation

```
gcc mc.c -lm
```

## Exact density of states

There is a Mathematica program `is2dos0.ma` of computing
the density of states
To compute the density of states of the n x m Ising model,
where n and m are two integers, type
```
    math < is2dos0.ma n m
```
Here `math` is the command-line interface of Mathematica.
It will produce the two text files

IsingDOSnxm.txt, for integer density of states.

is2lndosnxm.txt, for logarithm density of states.

The density of states will contain `n * m + 1` lines.
Each contains (the logarithm of) the density of states
of an energy level.  This is a discrete model,
so the energy levels are integers, and they are given by
```
-2*n*m,
-2*n*m + 4,
-2*n*m + 8,
...
2*n*m - 4,
2*n*m
```

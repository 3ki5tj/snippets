# Snippets: a collection of independent modules relevant to molecular simulations

## Overview

This project contains several modules mainly relevant to molecular simulations.

The modules within this project tend to be independent
with little or none interdependency on other modules.

The modules are usually written in C or JavaScript (for building static web pages).

The codes are not polished too much, as we want to expose
the essential parts of the algorithms as much as possible.
Each module often contains only a handleful files, and the reader can
readily locate the essential parts of the methods in a few
functions, and use them in their own projects.
The downside of this is that the codes may are not fully optimized
and may contain errors (sorry about that).

## Featured modules

* [basic modules in C](src/common/c) and [their JavaScript translations](src/js/modules/clib),
  such as [random number generators](src/common/c/rand),
  [solving the eigenvalues and eigenvetors of a real symmetric matrix](src/common/c/linalge/eig.h).
  [solving a set of linear equations](src/common/c/linalge/linsolve.h).

* [modules relevant to molecular simulations](src/molsimul).

  * Utility: building an extended spiral configuration
    from an amino-acid sequence
    and save it as a PDB file,
    in [C](src/molsimul/utilities/mkspx/c/mkspx.c),
    in [Python](src/molsimul/utilities/mkspx/python/mkspx.py),
    or as [a web app](src/molsimul/utilities/mkspx/web/mkspx.html).

  * Algorithms in molecular dynamics,
    including an implementation for
    [SETTLE: constraint algorithm for 3-point water molecule](src/molsimul/methods/settle/c),
    [calculating the Madelung constant for NaCl-like cubic lattices using Ewald sum](src/molsimul/methods/ewald/c/madelung.c).

## Github link

[Github link](https://github.com/3ki5tj/snippets)

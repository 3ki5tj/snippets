# Snippets: small code pieces

This is a collection of small code pieces.

## Overview

### Design principles

1. Being self-contained

    The module contains self-contained code.
    Basic modules are self-contained than those in zcom.h (common).
    The code is written in such a way that you can copy and paste
    some of the routines to a new project.

2. Written in C

    Aimed for rapid prototyping.
    C++ is not suitable, too much to take care of the language features.

    Cons:

    * Stick to double, unable to use template to shift between floating point types.

## Contents

* [basic](basic) Basic modules
* [molsimul](molsimul) Molecular simulations
* [misc](misc) Miscellaneous codes and applications

## Source code

[Github link](https://github.com/3ki5tj/snippets)

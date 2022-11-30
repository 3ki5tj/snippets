# Correction to the Langevin equation

## Usage

```sh
make
```

### Without correction

```sh
./langevin 0
```

### With correction

```sh
./langevin 1
```

```gnuplot
plot [][0:] "hist0.dat" u 1:3 w l, "hist1.dat" u 1:3 w l
```

"hist0.dat" should not be flat, "hist1.dat" should be flat.

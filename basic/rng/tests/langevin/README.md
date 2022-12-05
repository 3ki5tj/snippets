# Correction to the Langevin equation

## Langevin Equation

$$
dx = x \, dW
$$

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
plot [][0:] "hist0.dat" u 1:3 w l, "hist1.dat" u 1:3 w l, 2/x**2
```

The histogram "hist0.dat" from the uncorrected Langevin equation
should not be flat, but p(x) = 2/x**2,
while that "hist1.dat" from the corrected version should be flat.

## Theory

For the uncorrected Langevin equation,

$$ dx = x \, dW $$

This translate to a Fokker-Planck equation with $B(x) = x^2$.

$$
\frac{\partial}{\partial t} P(x, t)
=
\frac{1}{2}
\frac{\partial^2}{\partial x^2}
\left(
x^2 \,P(x)
\right).
$$

So the stationary distribution is

$$
P(x) \propto \frac{1}{x^2}.
$$

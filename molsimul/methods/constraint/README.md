
## Two-dimensional model `xy.c`

The infinite stiffness limt can depend on
a parameter alpha, depending on the interaction
between the two degrees of freedom.

### Theory

The potential is
$$
V(x, y) = (1/2) x^2  + (k/2) y^2 (1 + a x^2).
$$
Here, $k$ would tend to infinity ultimately.

With the scaling is $z = \sqrt{k} y$, we get
$$
V(x, z) = (1/2) (x^2  + z^2  + a x^2 z^2).
$$
which is symmetric with respect to $x$ and $z$.

Now by using the virial theorem, we have
\begin{align}
T
&= \langle x dV/dx \rangle \\
&= \langle x^2 \rangle  + a \langle x^2 z^2 \rangle \\
&\approx \langle x^2 \rangle \left( 1 + a \langle z^2 \rangle \right).
\end{align}

Let $t = \langle x^2 \rangle$, we have
$$
t = \frac{ 2 T }
{ 1 + \sqrt{1 + 4 a T } }.
$$

This means the average $\langle x^2 \rangle$
depends only on the coupling coefficient, $a$,
but not on the stiffness $k$.
So as the stiffness tends to the infinity,
$\langle x^2 \rangle$ is not uniquely defined.



## Spring model `lt.c`

Even without thermostat,
the infinite stiffness limit is
not the same as the result of constraints.


### Comparison in the canonical ensemble

With thermotat,
the constrained ensemble and unconstrained ensemble
produced distributions that differ by a factor of $l$.

```
plot [:1.3][:] "dat/l_c_nt.his" u 1:(log($2)) w l, "dat/l_nt.his" u 1:(log($2)) w l
```

#### Data

Compare `l_c.his` vs `l.his`

```
plot [:2][:] "dat/l_c.his" u 1:(log($2)) w l, "dat/l.his" u 1:(log($2/$1)) w l
```


## Trimer model `trimer.c`

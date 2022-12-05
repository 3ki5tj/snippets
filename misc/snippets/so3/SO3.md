# Uniform sampling of SO(3)

Idea: A normalized quaternion, $q = w + x i + y j + z k$,
with $x^2 + y^2 + z^2 + w^2 = 1$,
can be mapped to a rotation matrix:

$$
Q = \left(\begin{array}{ccc}
  1 - 2 y^2 - 2 z^2 & 2 x y - 2 z w    & 2 x z + 2 y w \\
  2 x y + 2 z w     & 1 - 2 x^2 - 2z^2 & 2 y z - 2 x w \\
  2 x z - 2 y w     & 2 y z + 2 x w    & 1 - 2 x^2 - 2 y^2
  \end{array}
\right)
$$

The quaternion can be generated from three $\mathcal U(0,1)$ uniform random numbers, $u_1$, $u_2$, $u_3$ as

$$q = \left(
    \sqrt{1-u_1} \sin(2\pi u_2),
    \sqrt{1-u_1} \cos(2\pi u_2),
    \sqrt{u_1} \sin(2\pi u_3),
    \sqrt{u_1} \cos(2\pi u_3)
\right).$$

The mapping from $SU(2)$ or quaternions to $SO(3)$ is a 2:1 surjection (from wikipedia).

References:

[1] [Map from unit quaternions to $SO(3)$](https://math.stackexchange.com/questions/1587309/map-from-unit-quaternions-to-so3)

[2] [The Quaternions and the Spaces $S^3$, $SU(2)$, $SO(3)$ and $RP^3$](https://www.cis.upenn.edu/~cis6100/geombchap8.pdf)

[3] [Random sampling on $SO(3)$](https://y7k4.github.io/2020/10/16/random-sampling-on-so-3.html)

[4] [HealPix](https://healpix.sourceforge.io/)

[5] [Healpy](https://github.com/healpy/healpy) [Healpy tutorial](https://healpy.readthedocs.io/en/latest/tutorial.html)

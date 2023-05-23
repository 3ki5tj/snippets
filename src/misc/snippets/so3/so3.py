#!/usr/bin/env python3

import numpy as np
import math

def make_rotation_matrix(u1, cphi, sphi, cpsi, spsi):
    ''' build a rotation matrix from u1, phi and psi

        Idea: A normalized quaternion, $q = w + x i + y j + z k$,
        with $x^2 + y^2 + z^2 + w^2 = 1$,
        can be mapped to a rotation matrix:

        $$
        Q = \left(
          \begin{array}{ccc}
          1 - 2*y^2 - 2*z^2 & 2xy - 2zw         & 2xz + 2yw       \\
          2xy + 2zw         & 1 - 2*x^2 - 2*z^2 & 2yz - 2xw       \\
          2xz - 2yw         & 2yz + 2xw         & 1 - 2*x^2 - 2*y^2
          \end{array}
        \right)
        $$

        The quaternion can be generated from three $\mathcal U(0,1)$ uniform random numbers, $u_1$, $u_2$, $u_3$ as

        $$
        q =
        \left(
            \sqrt{1-u_1} \sin(2 \pi u_2),
            \sqrt{1-u_1} \cos(2 \pi u_2),
            \sqrt{ u_1 } \sin(2 \pi u_3),
            \sqrt{ u_1 } \cos(2 \pi u_3)
        \right).
        $$

        The mapping from $SU(2)$ or quaternions to $SO(3)$ is a 2:1 surjection (from wikipedia).

        References:

        [1] https://math.stackexchange.com/questions/1587309/map-from-unit-quaternions-to-so3

        [2] https://www.cis.upenn.edu/~cis6100/geombchap8.pdf

        [3] https://y7k4.github.io/2020/10/16/random-sampling-on-so-3.html
    '''

    c1 = (1-u1)**0.5
    s1 = u1**0.5
    w = c1*cphi
    z = c1*sphi
    x = s1*cpsi
    y = s1*spsi
    return np.array([
        [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w,     2*x*z + 2*y*w    ],
        [2*x*y + 2*z*w,     1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w    ],
        [2*x*z - 2*y*w,     2*y*z + 2*x*w,     1 - 2*x*x - 2*y*y],
    ])



def so3_sampling_cubic(n, m):
    ''' making n**3 * (aspect_ratio) rotation matrices that are supposed to evenly
        cover SO(3)

        The ideal ratio of n:m is around 8:3
    '''

    inc1, inc2, inc3 = 0.5, 0.5, 0.5
    #inc1, inc2, inc3 = 2**0.5/2, 3**0.5/3, 5**0.5/5

    for phi in (np.arange(n)+inc1)*2*math.pi/n:
        cphi = math.cos(phi)
        sphi = math.sin(phi)
        for psi in (np.arange(n)+inc2)*2*math.pi/n:
            cpsi = math.cos(psi)
            spsi = math.sin(psi)
            for u1 in (np.arange(m)+inc3)/m:
                yield make_rotation_matrix(u1, cphi, sphi, cpsi, spsi)



def so3_sampling_fib(n):
    ''' making n**2 rotation matrices that are supposed to evenly
        cover the SO(3) manifold

        The u1 dimension is generated from the Fibonacci sequence
        '''

    fac = (5**0.5 - 1) / 2 # golden ratio
    #inc1, inc2, inc3 = 0.5, 0.5, 0.0
    inc1, inc2, inc3 = 2**0.5/2, 3**0.5/3, 5**0.5/5

    i = 0
    for phi in (np.arange(n)+inc1)*2*math.pi/n:
        cphi = math.cos(phi)
        sphi = math.sin(phi)
        for psi in (np.arange(n)+inc2)*2*math.pi/n:
            cpsi = math.cos(psi)
            spsi = math.sin(psi)
            i += 1
            u1 = math.fmod(i*fac + inc3, 1)
            yield make_rotation_matrix(u1, cphi, sphi, cpsi, spsi)



def make_random_vec():
    ''' make a random unit vector '''
    v = np.random.randn(3)
    # return the normalized vector
    return v / np.sum(v**2)**0.5


def random_vecs_test(rots, times=3, bins=10, title=""):
    ''' testing if the bilinear projection to two random vectors
        yield a uniform distribution in (0, 1) '''

    import matplotlib.pyplot as plt

    dx = (2.0/times/bins) * 0.8 # the last factor produces some margin between bins

    for t in range(times):
        if t == 0:
            u = v = np.array([1.0, 0.0, 0.0])
        elif t == 1:
            u = v = np.array([0.0, 1.0, 0.0])
        elif t == 2:
            u = v = np.array([0.0, 0.0, 1.0])
        else:
            u = make_random_vec()
            v = make_random_vec()

        xs = np.array([u @ rot @ v for rot in rots])
        plt.hist(xs+t*dx, bins=bins, range=(-1+t*dx,1+t*dx), width=dx, alpha=0.7)

    plt.title(title)
    plt.show()



def sample_and_test(sampling_method, name, **params):

    # generate rotation matrices
    rots = [rot for rot in sampling_method(**params)]
    print("generated", len(rots), "rotation matrices (" + name + ")")

    # test uniformity
    random_vecs_test(rots, title=name)



sample_and_test(so3_sampling_fib, "Fibonacci + square lattice", n = 100)

#sample_and_test(so3_sampling_cubic, "cubic lattice", n = 50, m = 20)


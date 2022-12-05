#!/usr/bin/env python3

import numpy as np
import healpy as hp
import math

def quaternion_to_rotation_matrix(q):
    ''' build a rotation matrix from a normalized quaternion,
        $q = w + x i + y j + z k$,
        with $x^2 + y^2 + z^2 + w^2 = 1$:

        $$
        Q = \left(
          \begin{array}{ccc}
          1 - 2*y^2 - 2*z^2 & 2xy - 2zw         & 2xz + 2yw       \\
          2xy + 2zw         & 1 - 2*x^2 - 2*z^2 & 2yz - 2xw       \\
          2xz - 2yw         & 2yz + 2xw         & 1 - 2*x^2 - 2*y^2
          \end{array}
        \right)
        $$
    '''

    x, y, z, w = q
    
    return np.array([
        [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w,     2*x*z + 2*y*w    ],
        [2*x*y + 2*z*w,     1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w    ],
        [2*x*z - 2*y*w,     2*y*z + 2*x*w,     1 - 2*x*x - 2*y*y],
    ]).transpose([2, 0, 1])



def angles_to_rotation_matrix(psi, theta, phi):
    ''' build a rotation matrix from the three angles, psi, theta, phi '''
    
    # build quaternions from the three angles, psi, theta, phi
    # Eq. (4) of Yershova Jain, LaValle, Mitchell
    # Generating Uniform Incremental Grids on SO(3) Using the Hopf Fibration
    # Int J Rob Res. (2010)
    x1 = np.cos(theta/2) * np.cos(psi/2)
    x2 = np.cos(theta/2) * np.sin(psi/2)
    x3 = np.sin(theta/2) * np.cos(phi + psi/2)
    x4 = np.sin(theta/2) * np.sin(phi + psi/2)

    return quaternion_to_rotation_matrix([x1, x2, x3, x4])



def make_rotation_matrices_fib(nside):
    ''' create 12 * nside**2 rotation matrices that evenly cover SO(3) '''

    npix = hp.nside2npix(nside) # 12 * nside**2
    
    # hp.pix2ang() returns theta, phi, for the S^2 dimension
    # theta goes from 0 to pi
    # phi goes from 0 to 2*pi
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    
    # psi represent the S^1 dimension
    golden_ratio = (5**0.5 - 1)/2
    psi = np.mod(np.arange(1, npix+1) * golden_ratio, 1) * (math.pi*2)

    return angles_to_rotation_matrix(psi, theta, phi)


#rots = make_rotation_matrices_fib(16)
#print(rots, len(rots))


def make_rotation_matrices_lattice(nside, m):
    ''' create 12 * nside**2 * m rotation matrices that evenly cover SO(3) '''

    npix = hp.nside2npix(nside) # 12 * nside**2
    
    # pix2ang returns theta, phi, for the S^2 dimension
    # theta goes from 0 to pi
    # phi goes from 0 to 2*pi
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    theta = np.repeat(theta, m)
    phi = np.repeat(phi, m)
 
    # psi represents the S^1 dimension
    psi = ((np.arange(m) + 0.5) / m) * (math.pi*2)
    psi = np.repeat(psi, npix)

    return angles_to_rotation_matrix(psi, theta, phi)



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



def sample_and_test(gridding_method, name, **params):

    # generate rotation matrices
    rots = gridding_method(**params)
    print("generated", len(rots), "rotation matrices (" + name + ")")

    # test uniformity
    random_vecs_test(rots, title=name)


sample_and_test(make_rotation_matrices_fib, "HealPix + Fibonacci", nside=32)

#sample_and_test(make_rotation_matrices_lattice, "HealPix lattice", nside=8, m=32)

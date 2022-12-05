#!/usr/bin/env python3

import numpy as np
import math



def make_x_rotation_matrix(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [ 1,  0,  0],
        [ 0,  c, -s],
        [ 0,  s,  c],
    ])



def make_y_rotation_matrix(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [ c,  0,  s],
        [ 0,  1,  0],
        [-s,  0,  c],
    ])



def make_z_rotation_matrix(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [ c, -s,  0],
        [ s,  c,  0],
        [ 0,  0,  1],
    ])



def euler_sampling_cubic(n):
    ''' making n**3 rotation matrices over SO(3) '''

    #inc1, inc2, inc3 = 0.5, 0.5, 0.5
    inc1, inc2, inc3 = 2**0.5/2, 3**0.5/3, 5**0.5/5

    for phi in (np.arange(n)+inc1)*2*math.pi/n:
        m1 = make_y_rotation_matrix(phi)
        for psi in (np.arange(n)+inc2)*2*math.pi/n:
            m2 = make_x_rotation_matrix(psi)
            for theta in (np.arange(n)+inc3)*2*math.pi/n:
                m3 = make_z_rotation_matrix(theta)
                yield m1 @ m2 @ m3



def euler_sampling_fib(n):
    ''' making n**2 rotation matrices over SO(3) '''

    fac = (5**0.5 - 1)/2

    #inc1, inc2, inc3 = 0.5, 0.5, 0.5
    inc1, inc2 = 2**0.5/2, 3**0.5/3

    i = 0
    for phi in (np.arange(n)+inc1)*2*math.pi/n:
        m1 = make_y_rotation_matrix(phi)
        for psi in (np.arange(n)+inc2)*2*math.pi/n:
            m2 = make_x_rotation_matrix(psi)
            i += 1
            theta = math.fmod(i * fac, 1)*2*math.pi
            m3 = make_z_rotation_matrix(theta)
            yield m1 @ m2 @ m3



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



def sample_and_test(sampling_method, param, name):

    # generate rotation matrices
    rots = [rot for rot in sampling_method(param)]
    print("generated", len(rots), "rotation matrices (" + name + ")")

    # test uniformity
    random_vecs_test(rots, title=name)



sample_and_test(euler_sampling_cubic, 36, "Euler-angles sampling, cubic lattice")

sample_and_test(euler_sampling_fib, 216, "Euler-angles sampling, Fib + square lattice")

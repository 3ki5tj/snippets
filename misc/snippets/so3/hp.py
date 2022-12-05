#!/usr/bin/env python3

''' healpy documentation: https://healpy.readthedocs.io/en/latest/tutorial.html

    pip3 install healpy

'''

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def sphere(nside):
    ''' @param nside should be a power of 2 '''
  
    print(
        "Approximate resolution at nside {} is {:.2} deg".format(
            nside, hp.nside2resol(nside, arcmin=True) / 60
        )
    )
    
    npix = hp.nside2npix(nside)
    print(npix)
    
    # pix2ang returns theta, phi
    # theta goes from 0 to pi
    # phi goes from 0 to 2*pi
    print (hp.pix2ang(nside, np.arange(npix)))
    print (hp.pix2vec(nside, 1))

    #m = np.arange(npix)
    #hp.mollview(m, title="Mollview image RING")
    #hp.graticule()
    #plt.show()
    
    
sphere(8)
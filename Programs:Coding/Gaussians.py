#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 09:36:28 2020

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from scipy.optimize import curve_fit

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-0266-51630-0053.fits')
redshift = test_data[2].data['Z']
lam = 10 ** (test_data[1].data['loglam'] / (1 + redshift))
flux = test_data[1].data['flux'] / (1 + redshift)
ivar = test_data[1].data['ivar']

#MgII properties:
MgII_bounds = (lam > 2700)&(lam < 2900)
MgII_flux = flux[MgII_bounds]
MgII_lam = lam[MgII_bounds]
MgII_ivar = ivar[MgII_bounds]
sig = 1 / np.sqrt(abs(MgII_ivar))

#C4 properties
C4_bounds = (lam > 1500) & (lam < 1600)
C4_flux = flux[C4_bounds]
C4_lam = lam[C4_bounds]
C4_ivar = ivar[C4_bounds]
sigma = 1 / np.sqrt(abs(C4_ivar))

def gaussian(x, m, sigma):
    g = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-.5 * ((x - m) / sigma)**2)
    return g

# m (mu) is the mean expectation (can be mode or median too)
# sigma is the standard deviation, and you can find it from FWHM
# x is the position of where you are
    
# in the Shein2011 paper, they found that for MgII, the FWHM is 1200 km s^-1

gauss_fit = gaussian(MgII_lam, 0 ,sig)




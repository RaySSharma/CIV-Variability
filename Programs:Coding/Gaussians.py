#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 09:36:28 2020

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7294-56739-0440.fits')
redshift = test_data[2].data['Z']
lam = (10 ** (test_data[1].data['loglam'])) / (1 + redshift)
flux = test_data[1].data['flux'] / (1 + redshift)
ivar = test_data[1].data['ivar']

#MgII properties:
MgII_bounds = (lam > 2700) & (lam < 2900)
MgII_flux = flux[MgII_bounds]
MgII_lam = lam[MgII_bounds]
MgII_ivar = ivar[MgII_bounds]
#sig = 1 / np.sqrt(abs(MgII_ivar))
sig = 1200 / 2.35482 #km/s

#C4 properties
C4_bounds = (lam > 1500) & (lam < 1600)
C4_flux = flux[C4_bounds]
C4_lam = lam[C4_bounds]
C4_ivar = ivar[C4_bounds]
#sig_C4 = 1 / np.sqrt(abs(C4_ivar))
sig_C4 = 24500 / 2.35482 #km/s

def gaussian(x, k, m, sigma):
    g = k * np.exp(-.5 * ((x - m) / sigma)**2)
    return g

def gaussian3(x, m1, sigma1, k1, m2, sigma2, k2, m3, sigma3, k3):
    gauss = (k1 * np.exp(-.5 * ((x - m1) / sigma1)**2)) + (k2 * np.exp(-.5 * ((x - m2) / sigma2)**2)) + (k3 * np.exp(-.5 * ((x - m3) / sigma3)**2))
    return gauss

# m (mu) is the mean expectation (can be mode or median too)
# sigma is the standard deviation, and you can find it from FWHM
# x is the position of where you are
    
# in the Shein2011 paper, they found that for MgII, the FWHM is 1200 km s^-1

MgII_conditions = (1000, 2798, sig) #this is for the MgII line
gauss_fit, pcov = curve_fit(gaussian, MgII_lam, MgII_flux, p0 = MgII_conditions)

CIV_conditions = (1549, sig_C4, 1000, 2549, sig_C4 / 2, 500, 3549, 2 * sig_C4, 2000)
C4_gauss_fit, pcov2 = curve_fit(gaussian3, C4_lam, C4_flux, p0 = CIV_conditions)

plt.figure()
plt.plot(MgII_lam, gaussian(MgII_lam, *gauss_fit), 'b', label = 'Gaussian Curve')
plt.plot(lam, flux, 'r', label = 'Continuum')
plt.legend(loc = 'best')
plt.title('MgII Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(2600, 3000)
plt.ylim(-2, 2)

plt.figure()
plt.plot(C4_lam, gaussian3(C4_lam, *C4_gauss_fit), 'k', label = 'Gaussian Curve')
plt.plot(lam, flux, 'y', label = 'Continuum')
plt.legend(loc = 'best')
plt.title('CIV Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(1400, 1700)
plt.ylim(-1, 2.5)





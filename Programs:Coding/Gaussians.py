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
sig = 1200 / 2.35482 #A/s,

#C4 properties
C4_bounds = (lam > 1500) & (lam < 1600)
C4_flux = flux[C4_bounds]
C4_lam = lam[C4_bounds]
C4_ivar = ivar[C4_bounds]
sig_C4 = 245000 / 2.35482 #A/s, adjust these to wavelength (Angstroms)

def gaussian(x, k, m, sigma):
    g = k * np.exp(-.5 * ((x - m) / sigma)**2)
    return g

def gaussian3(x, m, sigma1, k1, sigma2, k2, sigma3, k3):
    gauss = (k1 * np.exp(-.5 * ((x - m) / sigma1)**2)) + (k2 * np.exp(-.5 * ((x - m) / sigma2)**2)) + (k3 * np.exp(-.5 * ((x - m) / sigma3)**2))
    return gauss

# m (mu) is the mean expectation (can be mode or median too)
# sigma is the standard deviation, and you can find it from FWHM
# x is the position of where you are
    
# in the Shein2011 paper, they found that for MgII, the FWHM is 1200 km s^-1

MgII_conditions = (1000, 2798, sig) #this is for the MgII line
gauss_fit, pcov = curve_fit(gaussian, MgII_lam, MgII_flux, p0 = MgII_conditions)

CIV_condition1 = (1549, sig_C4 / 5000, 100, sig_C4, 0, sig_C4, 0)
C4_gauss_fit1, pcov2 = curve_fit(gaussian3, C4_lam, C4_flux, p0 = CIV_condition1)
CIV_condition2 = (1549, sig_C4, 0, sig_C4 / 50000 , 100, sig_C4, 0)
C4_gauss_fit2, pcov3 = curve_fit(gaussian3, C4_lam, C4_flux, p0 = CIV_condition2)
CIV_condition3 = (1549, sig_C4, 0, sig_C4, 0, sig_C4 / 10000, 1000)
C4_gauss_fit3, pcov4 = curve_fit(gaussian3, C4_lam, C4_flux, p0 = CIV_condition3)

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
plt.plot(C4_lam, gaussian3(C4_lam, *C4_gauss_fit1), 'k', label = 'Gaussian Curve 1')
plt.plot(lam, flux, 'y')
plt.legend(loc = 'best')
plt.title('CIV Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(1400, 1700)
plt.ylim(-0.5, 3)

plt.plot(C4_lam, gaussian3(C4_lam, *C4_gauss_fit2), 'm', label = 'Gaussian Curve 2')
plt.plot(lam, flux, 'y')
plt.legend(loc = 'best')
plt.title('CIV Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(1400, 1700)
plt.ylim(-0.5, 3)

plt.plot(C4_lam, gaussian3(C4_lam, *C4_gauss_fit3), 'c', label = 'Gaussian Curve 3')
plt.plot(lam, flux, 'y', label = 'Continuum')
plt.legend(loc = 'best')
plt.title('CIV Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(1400, 1700)
plt.ylim(-0.5, 3)


# Change the sigma to delta(lambda), this will probably change the different 
# heights of the curves. 
# put bounds on k: 0 - 1000
# put bounds on sigma: depending on C4 or MgII, MgII: 5000 km/s, C4: 15000 as 
# upper bounds. Remember all of these must be positive
# if after all these corrections you find that the fit is still looking strange,
# you can always adjust mu again to see if it can go outside the boundaries.


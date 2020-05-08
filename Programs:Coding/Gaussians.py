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

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7665-57328-0452.fits')
redshift = test_data[2].data['Z']
lam = (10 ** (test_data[1].data['loglam'])) / (1 + redshift)
flux = test_data[1].data['flux'] / (1 + redshift)
ivar = test_data[1].data['ivar']
c = 3 * 10**5

#MgII properties:
MgII_bounds = (lam > 2700) & (lam < 2900)
MgII_flux = flux[MgII_bounds]
MgII_lam = lam[MgII_bounds]
MgII_ivar = ivar[MgII_bounds]
sig = (1200 * 2790) / (c) #A

#C4 properties
C4_bounds = (lam > 1500) & (lam < 1600)
C4_flux = flux[C4_bounds] 
C4_lam = lam[C4_bounds]
C4_ivar = ivar[C4_bounds]
sig_C4 = (24500 * 1549) / (c) #A

def gaussian(x, m, sigma, k):
    sigma = (sigma / c) * m
    g = k * np.exp(-.5 * ((x - m) / sigma)**2)
    return g

def gaussian3(x, m, sigma1, k1, sigma2, k2, sigma3, k3):
    gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
    return gauss
    
# in the Shein2011 paper, they found that for MgII, the FWHM is 1200 km s^-1

MgII_conditions = (2798, sig, 100)
MgII_boundaries = [[2700, 0, 0], [2900, 5000, 1000]]
gauss_fit, pcov = curve_fit(gaussian, MgII_lam, MgII_flux, p0 = MgII_conditions, bounds = MgII_boundaries)

CIV_condition = (1549, 500, 100, 1000, 100, 10000, 100)
C4_boundaries = [[1500, 0, 0, 0, 0, 0, 0], [1600, np.inf, 1000, np.inf, 1000, np.inf, 1000]]
C4_gauss_fit, pcov2 = curve_fit(gaussian3, C4_lam, C4_flux, p0 = CIV_condition, bounds = C4_boundaries)


print(pcov, pcov2)
print(gauss_fit, C4_gauss_fit)

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
plt.plot(C4_lam, gaussian3(C4_lam, *C4_gauss_fit), label = 'CIV Full')
plt.plot(C4_lam, gaussian(C4_lam, *C4_gauss_fit[[0, 1, 2]]), 'k', label = 'Gaussian Curve 1')
plt.plot(C4_lam, gaussian(C4_lam, *C4_gauss_fit[[0, 3, 4]]), 'b', label = 'Gaussian Curve 2')
plt.plot(C4_lam, gaussian(C4_lam, *C4_gauss_fit[[0, 5, 6]]), 'm', label = 'Gaussian Curve 3')
plt.plot(lam, flux, 'y', alpha = 0.5, label = 'Continuum')
plt.legend(loc = 'best')
plt.title('CIV Gaussian Fit')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux')
plt.xlim(1400, 1700)
plt.ylim(-0.5, 3)

# put bounds on k: 0 - 1000
# put bounds on sigma: depending on C4 or MgII, MgII: 5000 km/s, C4: 15000 as 
# upper bounds. Remember all of these must be positive
# if after all these corrections you find that the fit is still looking strange,
# you can always adjust mu again to see if it can go outside the boundaries.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:40:22 2020

@author: RachelCampo
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.io import fits
from astropy import units as u
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15 as p15
import astropy.constants as const
from scipy.interpolate import UnivariateSpline

import FESubtraction
import Gaussians2
import MasterCode

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-5407-55926-0636.fits')
redshift = test_data[2].data['Z']
wavelength = 10 ** test_data[1].data['loglam'] / (1 + redshift)
flux = test_data[1].data['flux']
c = 3 * 10 ** 5

# pulling out the pf values for the gaussian CIV from table:

pf = final_list[1:,15]

def gaussian(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k * np.exp(-.5 * ((x - m) / sigma)**2)
        return g

def gaussian3(x, m, sigma1, k1, sigma2, k2, sigma3, k3):
        gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
        return gauss

C4_wav = np.linspace(1500, 1600, 1000)
C4_flux = gaussian3(C4_wav, *pf)

d = p15.luminosity_distance(redshift).to('cm')
C4_flux = C4_flux * 10 ** -17 * u.erg / u.s / u.cm / u.cm / u.Angstrom
C4_wav = C4_wav * u.Angstrom

def FWHM(x, y):
    spline = UnivariateSpline(x, y - y.max() / 2, s = 0)
    a, b = spline.roots()
    return abs(a - b), a, b

def A_to_kms(fwhm, m):
    return const.c * fwhm / m

fwhm, left, right = FWHM(C4_wav, C4_flux)
fwhm = fwhm * u.Angstrom
fwhm_kms = A_to_kms(fwhm, pf[0] * u.Angstrom)
print(fwhm, fwhm_kms.to('km/s'))

C4_luminosity = C4_flux * (4 * np.pi * d ** 2)

def C4_lum(wav, lum):
    return np.trapz(lum, wav)

C4_L = C4_lum(C4_wav, C4_luminosity)
print(C4_L)

def L1350(C4_L):
    return 10 ** (7.66 + 0.863 * np.log10(C4_L / (u.erg / u.s))) * u.erg / u.s

L_1350 =L1350(C4_L)
print(L_1350)

def mass_bh(lum, fwhm, a = 0.660, b = 0.53):
    return 10 ** (a + b * np.log10(lum / (1e44 * u.erg / u.s))
                  + 2 * np.log10(fwhm / (u.km / u.s))) * u.solMass
                  

#will most likely put this into a for loop/function in order to iterate
#over every quantity
lam = 1549
delt_lam = np.abs(final_list[1:,16] - lam)
v = (delt_lam * c) / lam

mass_BH = mass_bh(L_1350, v)
print(np.log10(mass_BH / u.solMass), mass_BH)






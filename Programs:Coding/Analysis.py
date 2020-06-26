#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:40:22 2020

@author: RachelCampo
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15 as p15
import astropy.constants as const
from scipy.interpolate import UnivariateSpline
import pandas as pd
import pdb

import FESubtraction
import Gaussians2
import MasterCode

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7665-57328-0452.fits')
redshift = test_data[2].data['Z']
wavelength = 10 ** test_data[1].data['loglam'] / (1 + redshift)
flux = test_data[1].data['flux']
c = 3 * 10 ** 5

final_list = pd.read_csv('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/final_list.csv', sep = ',', index_col = False)
Q = len(final_list)

#C4pf = final_list.loc[:, ['CIV Gaussian Fit']]
#Mg2pf = final_list.loc[:, ['MgII Gaussian Fit']]
mu = final_list.loc[:, 'CIV Mu Value from Gaussian Fitting'].values
sig1 = final_list.loc[:, 'CIV Sigma 1 Value from Gaussian Fitting'].values
c4k1 = final_list.loc[:, 'CIV K1 Value from Gaussian Fitting'].values
sig2 = final_list.loc[:, 'CIV Sigma 2 Value from Gaussian Fitting'].values
c4k2 = final_list.loc[:, 'CIV K2 Value from Gaussian Fitting'].values
sig3 = final_list.loc[:, 'CIV Sigma 3 Value from Gaussian Fitting'].values
c4k3 = final_list.loc[:, 'CIV K3 Value from Gaussian Fitting'].values

#print(mu, sig1, c4k1, sig2, c4k2, sig3, c4k3)
def gaussian(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k * np.exp(-.5 * ((np.ones((Q, x)) - m.reshape(Q,1)) / sigma)**2)
        return g

def gaussian3(x, m, sigma1, k1, sigma2, k2, sigma3, k3):
        gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
        return gauss

C4_wav = np.linspace(1500, 1600, 1000)
C4_flux = gaussian3(C4_wav, mu, sig1, c4k1, sig2, c4k2, sig3, c4k3)

d = p15.luminosity_distance(redshift).to('cm')
C4_flux = C4_flux * 10 ** -17 * u.erg / u.s / u.cm / u.cm / u.Angstrom
C4_wav = C4_wav * u.Angstrom

def FWHM(x, y):
    spline = UnivariateSpline(x, y - y.max() / 2, s = 0)
    a, b = spline.roots()
    return abs(a - b), a, b
#abs(a-b) is actually delt_lam, it rotational velocity

def A_to_kms(fwhm, m):
    return const.c * fwhm / m

fwhm, left, right = FWHM(C4_wav, C4_flux)
fwhm = fwhm * u.Angstrom
fwhm_kms = A_to_kms(fwhm, mu * u.Angstrom)
print(fwhm, fwhm_kms.to('km/s'))

C4_luminosity = C4_flux * (4 * np.pi * d ** 2)

def C4_lum(wav, lum):
    return np.trapz(lum, wav)

C4_L = C4_lum(C4_wav, C4_luminosity)
print(C4_L)

def L1350(C4_L):
    return 10 ** (7.66 + 0.863 * np.log10(C4_L / (u.erg / u.s))) * u.erg / u.s

L_1350 = L1350(C4_L)
print(L_1350)

def mass_bh(lum, fwhm, a = 0.660, b = 0.53):
    return 10 ** (a + b * np.log10(lum / (1e44 * u.erg / u.s))
                  + 2 * np.log10(fwhm / (u.km / u.s))) * u.solMass
                  

#will most likely put this into a for loop/function in order to iterate
#over every quantity
lam = 1549
mu = final_list.loc[:, 'CIV Mu Value from Gaussian Fitting']
delt_lam = np.abs(mu - lam) # this may not be true, this could
# be in a different position
# use the abs(a-b) for rotational, what you calculated here is the translational
v = (delt_lam * c) / lam

mass_BH = mass_bh(L_1350, v)
print(np.log10(mass_BH / u.solMass), mass_BH)






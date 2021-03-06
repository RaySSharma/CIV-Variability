#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:40:22 2020

@author: RachelCampo
"""

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as p15
import astropy.constants as const
from scipy.interpolate import UnivariateSpline
import pandas as pd


#test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0105.fits')
#redshift = test_data[2].data['Z']
#wavelength = 10 ** test_data[1].data['loglam'] / (1 + redshift)
#flux = test_data[1].data['flux']
c = 3 * 10 ** 5

final_list = pd.read_csv('data2/rlc186/QuasarData/final_list.csv', sep = ',', index_col = False)
Q = len(final_list)

mu = final_list.loc[:, 'CIV Mu Value from Gaussian Fitting'].values
mu_err = final_list.loc[:, 'Error of CIV Mu from Gaussian Fitting'].values
sig1 = final_list.loc[:, 'CIV Sigma 1 Value from Gaussian Fitting'].values
sig1_err = final_list.loc[:, 'Error of CIV Sigma 1 from Gaussian Fitting'].values
c4k1 = final_list.loc[:, 'CIV K1 Value from Gaussian Fitting'].values
c4k1_err = final_list.loc[:, 'Error of CIV K1 from Gaussian Fitting'].values
sig2 = final_list.loc[:, 'CIV Sigma 2 Value from Gaussian Fitting'].values
sig2_err = final_list.loc[:, 'Error of CIV Sigma 2 from Gaussian Fitting'].values
c4k2 = final_list.loc[:, 'CIV K2 Value from Gaussian Fitting'].values
s4k2_err = final_list.loc[:, 'Error of CIV K2 from Gaussian Fitting'].values
sig3 = final_list.loc[:, 'CIV Sigma 3 Value from Gaussian Fitting'].values
sig3_err = final_list.loc[:, 'Error of CIV Sigma 3 from Gaussian Fitting'].values
c4k3 = final_list.loc[:, 'CIV K3 Value from Gaussian Fitting'].values
c4k3_err = final_list.loc[:, 'Error of CIV K3 from Gaussian Fitting'].values
redshift = final_list.loc[:, 'Redshift'].values

#print(mu, sig1, c4k1, sig2, c4k2, sig3, c4k3, redshift)

def gaussian(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k.reshape(Q, 1) * np.exp(-.5 * ((x * np.ones((Q, len(x))) - m.reshape(Q,1)) / sigma.reshape(Q,1))**2)
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
    #need loop through C4_flux and take out each individual row and feed that 
    #through FWHM
    a, b = spline.roots()
    return abs(a - b), a, b
#abs(a-b) is actually delt_lam, it rotational velocity

def A_to_kms(fwhm, m):
    return const.c * fwhm / m

for x in C4_flux:
    list = x

fwhm, left, right = FWHM(C4_wav, list) #may not actually need left and right value
fwhm = fwhm * u.Angstrom
fwhm_kms = A_to_kms(fwhm, mu * u.Angstrom)
print(fwhm, fwhm_kms.to('km/s'))

C4_luminosity = (C4_flux.T * (4 * np.pi * d ** 2)).T # transpose C4_flux

#may have to loop through again
def C4_lum(wav, lum):
    return np.trapz(lum, wav)

C4_L = C4_lum(C4_wav, C4_luminosity)
print(C4_L)

def L1350(C4_L):
    return 10 ** (7.66 + 0.863 * np.log10(C4_L / (u.erg / u.s))) * u.erg / u.s

L_1350 = L1350(C4_L)
print(L_1350)

#may not need to loop through anything for this one
def mass_bh(lum, fwhm, a = 0.660, b = 0.53):
    return 10 ** (a + b * np.log10(lum / (1e44 * u.erg / u.s))
                  + 2 * np.log10(fwhm / (u.km / u.s))) * u.solMass
                  
#def other_mass_bh(lum, fwhm):
#    lam_0 = np.trapz(lam * P) / (np.trapz(P)) #P is value of flux at that lambda
#    #we got this from Gaussians, the linspace
#    sigma_line = (np.trapz((lam - lam_0)**2 * P) / np.trapz(P))**(1/2)
#    return sigma_line
#will most likely put this into a for loop/function in order to iterate
#over every quantity
lam = 1549
mu = final_list.loc[:, 'CIV Mu Value from Gaussian Fitting']
delt_lam = np.abs(mu - lam) # this may not be true, this could
# be in a different position
# use the abs(a-b) for rotational, what you calculated here is the translational
v = (delt_lam * c) / lam

mass_BH = mass_bh(L_1350, fwhm_kms)
print(np.log10(mass_BH / u.solMass), mass_BH)






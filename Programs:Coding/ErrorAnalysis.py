#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:53:09 2020

@author: RachelCampo
"""

#this is going to be the definition of the error analysis, the propogation of
#uncertainty

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.integrate as spi
from uncertainties import unumpy
import pandas as pd
import Analysis2
import Gaussians
from astropy import units as u
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15 as p15
import astropy.constants as const
from scipy.interpolate import UnivariateSpline


#for the error, you need to calculate lum, fwhm, a, b.

test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0105.fits')
redshift = test_data[2].data['Z']
wavelength = 10 ** test_data[1].data['loglam'] / (1 + redshift)
flux = test_data[1].data['flux']
c = 3 * 10 ** 5



lum, fwhm, bhm = Analysis2.analysis()
final_list = pd.read_csv('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/final_list.csv', sep = ',', index_col = False)
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
c4k2_err = final_list.loc[:, 'Error of CIV K2 from Gaussian Fitting'].values
sig3 = final_list.loc[:, 'CIV Sigma 3 Value from Gaussian Fitting'].values
sig3_err = final_list.loc[:, 'Error of CIV Sigma 3 from Gaussian Fitting'].values
c4k3 = final_list.loc[:, 'CIV K3 Value from Gaussian Fitting'].values
c4k3_err = final_list.loc[:, 'Error of CIV K3 from Gaussian Fitting'].values


#We are also going to assume that a,b have 0 uncertainty (although unlikely)

#lum_err = unumpy.uarray([lum, lum], [0.009, 0.41]) #first is (m,b) and (m_err, b_err)

#fwhm_err = unumpy.uarray([fwhm], [C4_ivar]) #uses the flux variance from Gaussians2



def gaussian(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k.reshape(Q, 1) * unumpy.exp(-.5 * ((x * np.ones((Q, len(x))) - m.reshape(Q,1)) / sigma.reshape(Q,1))**2)
        return g

def gaussian3(x, unp_array):
        m, sigma1, k1, sigma2, k2, sigma3, k3 = unp_array
        gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
        return gauss



C4_wav = np.linspace(1500, 1600, 1000)
unp_array = unumpy.uarray([mu, sig1, c4k1, sig2, c4k2, sig3, c4k3],
                         [mu_err, sig1_err, c4k1_err, sig2_err, c4k2_err, sig3_err, c4k3_err])

C4_flux = gaussian3(C4_wav, unp_array)
d = p15.luminosity_distance(redshift).to('cm')
#C4_flux = C4_flux * 10 ** -17 * u.erg / u.s / u.cm / u.cm / u.Angstrom
#C4_wav = C4_wav * u.Angstrom

def FWHM(x, y):
    spline = UnivariateSpline(x, y - y.max() / 2, s = 0)
    a, b = spline.roots()
    return abs(a - b), a, b

def A_to_kms(fwhm, m):
    return const.c * fwhm / m

for x in C4_flux:
    list = x

C4_luminosity = C4_flux * (4 * np.pi * d ** 2)

def C4_lum(wav, lum):
    return np.trapz(lum, wav)

C4_L = C4_lum(C4_wav, C4_luminosity) #luminosity error

#doing FWHM manually:

mu_new = np.random.normal(loc = mu[i], scale = mu_err[i], size = 1000)
sig1_new = np.random.normal(loc = sig1, scale = sig1_err, size = 1000)
sig2_new = np.random.normal(loc = sig2, scale = sig2_err, size = 1000)
sig3_new = np.random.normal(loc = sig3, scale = sig3_err, size = 1000)
c4k1_new = np.random.normal(loc = c4k1, scale = c4k1_err, size = 1000)
c4k2_new = np.random.normal(loc = c4k2, scale = c4k2_err, size = 1000)
c4k3_new = np.random.normal(loc = c4k3, scale = c4k3_err, size = 1000)

flux_new = [gaussian3(C4_wav, [mu_new[i], sig1_new[i], c4k1_new[i], sig2_new[i],
            c4k2_new[i], sig3_new[i], c4k3_new[i]]) for i in range(1000)]
    
fwhm_new = [FWHM(C4_wav, x) for x in flux_new] #fwhm err

#finding black hole mass error:

def mass_bh(lum, fwhm, a = 0.660, b = 0.53):
    return 10 ** (a + b * np.log10(lum / (1e44 * u.erg / u.s))
                  + 2 * np.log10(fwhm / (u.km / u.s))) * u.solMass

bhm_err = mass_bh(C4_L, fwhm_new) #bhm_err


#hdr = ['Name', 'MJD', 'Fiber ID', 'Plate', 
#       'Error of Luminosity Using Bisector: Slope (m)', 'Error of Luminosity Using Bisector: (b)',
#       'Error of Luminosity Using Y|X: Slope (m)', 'Error of Luminosity Using Y|X: b',
#       'Error of FWHM',
#       'Error of Black Hole Mass: Bisector', 'Error of Black Hole Mass: Y|X']
#
#error_list = [quasar[2].data['THING_ID'], quasar[2].data['MJD'], 
#              quasar[2].data['FIBERID'], quasar[2].data['PLATE'],
#              lum_err_m, lum_err_b, lum_err_m2, lum_err_b2, fwhm_err,
#              bhm_err, bhm_err2
#
#data_frame = pd.DataFrame(error_list, columns = hdr, dtype = str)
#data_frame.to_csv('error_list.csv', index = False)


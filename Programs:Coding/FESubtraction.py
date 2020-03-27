#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 20:34:48 2020

@author: RachelCampo
"""
# this code will be the example code to see if it works on one quasar before
# applying it to the quasar catalog

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from scipy.optimize import curve_fit

# first I will import from the SDSS the platelist file, which containts all the
# plates used on the telescope.
# then I will import 
platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/platelist.fits')
specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')
test_data = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7572-56944-0012.fits')
FE_Template = pd.read_csv('/Users/RachelCampo/Desktop/Research/CIV-Variability/CIV-Variability-master/Fe_UVtemplt_A.dat', delim_whitespace = True)

# for the iron template, from the Vestergaard paper, the units are the same
# as in the SDSS

#properties from quasar
wavelength = test_data[1].data['loglam']
redshift = test_data[2].data['Z']
flux_rf = test_data[1].data['flux'] / (1 + redshift) #no need to worry about rf the flux
wavelength_rf = (10**wavelength) / (1 + redshift)
ivar = test_data[1].data['ivar']

#properties from iron template
FE_wavelength = FE_Template['wavelength'].values
FE_flux = FE_Template['flux'].values

#C4 properties
C4_bounds = (wavelength_rf > 1435) & (wavelength_rf < 1710)
C4_flux = flux_rf[C4_bounds]
C4_wavelength = wavelength_rf[C4_bounds]
C4_ivar = ivar[C4_bounds]
sigma = 1 / np.sqrt(abs(C4_ivar))


#the three parameters for the widening of the iron plate, rebinning, and fitting
#the function
def gauss(x, m, sig):
    sigma_conv = np.sqrt(sig**2 - 900**2) / (2 * np.sqrt(2 * np.log(2)))
    broadened_sigma = np.exp(- (x - m)**2 / (2 * (sigma_conv)**2))
    return broadened_sigma

def rebin_log(x, y):
    log_x = np.log10(x)
    new_x = np.logspace(log_x[1], log_x[-2], len(x))
    return new_x, interp1d(x, y)

def rebin_lin(x, y):
    x_new = np.linspace(x[1], x[-2], len(x))
    return x_new, interp1d(x, y)

log_FE_wavelength, log_FE_spline = rebin_log(FE_wavelength, FE_flux)

def fit_func(lam, A, k, B, m, sigma):
    global log_wavelength
    global ix
    FE_convolution = np.convolve(log_FE_spline(log_wavelength), gauss(log_wavelength, m, sigma), mode = 'same')
    return (A * lam**k) + (10**B * FE_convolution[ix])
#the fit_func is fitting both the continuum AND the FE plate! That's why we are adding
    #both the Alambda^K and the FE.

boundaries = [[0.1, -2, 10, 1500, 905], [1e3, 2, 20, 1600, 10000]]
p0 = [10, 0, 14, 1550, 1000]

log_wavelength, log_flux = rebin_log(C4_wavelength, C4_flux)
ix = ((log_wavelength > 1435)&(log_wavelength < 1465)) | ((log_wavelength > 1690)&(log_wavelength < 1710))
pf, covariances = curve_fit(fit_func, log_wavelength[ix], log_flux(log_wavelength[ix]), sigma = sigma[ix], bounds = boundaries, p0 = p0)


C4_cutoffs = (log_wavelength > 1465) & (log_wavelength < 1710) # did not use the pf variables, possible reason
# why it was not subtracted off correctly.
continuum_flux = fit_func(log_wavelength[ix], *pf) #put *pf in for the numbers

#this will be in logspace
plt.plot(log_wavelength[ix], log_flux(log_wavelength[ix]), label = 'Original')
plt.plot(log_wavelength[ix], continuum_flux, label = 'Continuum + Iron')
#plt.plot(log_wavelength, log_flux(log_wavelength) - continuum_flux, label = 'Subtracted Continuum')
plt.legend(loc = 'best')

plt.plot(log_wavelength, log_flux(log_wavelength), label = 'Original')
plt.plot(log_wavelength, continuum_flux, label = 'Continuum + Iron')
plt.plot(log_wavelength, log_flux(log_wavelength) - continuum_flux, label = 'Subtracted')
plt.xlabel('wavelength')
plt.ylabel('flux')
plt.legend(loc = 'best')

#converting back to linear space
subt_log_flux = log_flux(log_wavelength) - continuum_flux
linear_wavelength, linear_flux = rebin_lin(log_wavelength, subt_log_flux)

plt.plot(log_wavelength, log_flux(log_wavelength) - continuum_flux, label = 'Subtracted Continuum in log space')
plt.plot(linear_wavelength, linear_flux(log_wavelength), label = 'Subtracted Continuum in lin space')
plt.xlabel('flux')
plt.ylabel('wavelength')
plt.legend(loc = 'best')






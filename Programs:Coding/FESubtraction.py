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
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.optimize import curve_fit

platelist = fits.open('/home/ray/Downloads/platelist.fits')
specdatalist = fits.open('/home/ray/Downloads/DR14Q_v4_4.fits')
test_data = fits.open('/home/ray/Downloads/spec-7572-56944-0012.fits')
FE_Template = pd.read_csv('../Fe_UVtemplt_A.dat', delim_whitespace = True)

# for the iron template, from the Vestergaard paper, the units are the same
# as in the SDSS

#properties from quasar
wavelength = test_data[1].data['loglam']
redshift = test_data[2].data['Z']
flux_rf = test_data[1].data['flux'] / (1 + redshift) #no need to worry about rf the flux
wavelength_rf = (10**wavelength) / (1 + redshift)
ivar = test_data[1].data['ivar']
sigma = 1 / np.sqrt(abs(ivar))

#properties from iron template
FE_wavelength = FE_Template['wavelength'].values
FE_flux = FE_Template['flux'].values


#the three parameters for the widening of the iron plate
def gauss(x, m, sigma):
    sigma_conv = np.sqrt(sigma**2 - 900**2) / (2 * np.sqrt(2 * np.log(2)))
    broadened_sigma = np.exp(- (x - m)**2 / (2 * (sigma_conv)**2))
    return broadened_sigma

def rebin_log(x, y):
    log_x = np.log10(x)
    new_x = np.logspace(log_x[1], log_x[-2], len(x))
    return new_x, IUS(x, y)

log_FE_wavelength, log_FE_spline = rebin_log(FE_wavelength, FE_flux)

def fitting_function(lam, A, k, B, mu, sigma):
    global log_wavelength
    global cutoffs
    iron_convolution = np.convolve(log_FE_spline(log_wavelength), gauss(log_wavelength, mu, sigma), mode='same')
    return (A * lam**k) + (10**B * iron_convolution[cutoffs])

#the fit_func is fitting both the continuum AND the FE plate! That's why we are adding
    #both the Alambda^K and the FE.

cutoffs = ((wavelength_rf > 1435)&(wavelength_rf < 1465)) | ((wavelength_rf > 1690)&(wavelength_rf < 1710))
boundaries = [[0, -10, 0, 1000, 10], [100, 10, 20, 2000, 10000]]
p0 = [10, -2, 13, 1000, 1000]
log_wavelength, log_flux = rebin_log(wavelength_rf, flux_rf)
pf, covariances = curve_fit(fitting_function, log_wavelength[cutoffs], log_flux(log_wavelength[cutoffs]), sigma = sigma[cutoffs], bounds = boundaries, p0 = p0)


cutoffs = (log_wavelength > 1465) & (log_wavelength < 1710) # did not use the pf variables, possible reason
# why it was not subtracted off correctly.
continuum_flux = fitting_function(log_wavelength[cutoffs], *pf) #put *pf in for the numbers

plt.plot(wavelength_rf[cutoffs], flux_rf[cutoffs], label = 'Original')
plt.plot(wavelength_rf[cutoffs], continuum_flux, label = 'Continuum + Iron')
plt.plot(wavelength_rf[cutoffs], flux_rf[cutoffs] - continuum_flux, label = 'Subtracted Continuum')
plt.legend(loc = 'best')



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
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from scipy.optimize import curve_fit

def FE_sub(lam, f, redshift, ivar):

    FE_Template = pd.read_csv('/data2/rlc186/QuasarData/Fe_UVtemplt_A.dat', delim_whitespace = True)

    #properties from quasar
    wavelength = lam
    flux_rf = f #no need to worry about rf the flux
    wavelength_rf = (wavelength) / (1 + redshift)
    var = 1 / ivar

    #properties from iron template
    FE_wavelength = FE_Template['wavelength'].values
    FE_flux = FE_Template['flux'].values

    #C4 properties
    C4_bounds = (wavelength_rf > 1445) & (wavelength_rf < 1705)
    C4_flux = flux_rf[C4_bounds]
    C4_wavelength = wavelength_rf[C4_bounds]
    C4_var = var[C4_bounds]
    sigma = np.sqrt(abs(C4_var))
    
    #MgII properties
    MgII_bounds = (wavelength_rf > 2200) & (wavelength_rf < 3090)
    MgII_flux = flux_rf[MgII_bounds]
    MgII_wavelength = wavelength_rf[MgII_bounds]
    MgII_var = var[MgII_bounds]
    MgII_sigma = np.sqrt(abs(MgII_var))


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

    
    def fit_func_mgii(lam, A, k, B, m, sigma):
        nonlocal log_wavelength
        nonlocal ix
        FE_convolution = np.convolve(log_FE_spline(log_wavelength), gauss(log_wavelength, m, sigma), mode = 'same')
        if ix is not None:
            return (A * lam**k) + (10**B * FE_convolution[ix])
        else:
            return (A * lam**k) + (10**B * FE_convolution)
    
    def fit_func_civ(lam, A, k):
        return (A * lam**k)

    log_FE_wavelength, log_FE_spline = rebin_log(FE_wavelength, FE_flux)

 
    #the fit_func is fitting both the continuum AND the FE plate! That's why we are adding
        #both the Alambda^K and the FE.

    boundaries = [[0, -2], [1000, 2]]
    p0 = [10, 0]
    ix = ((C4_wavelength > 1445)&(C4_wavelength < 1465)) | ((C4_wavelength > 1700)&(C4_wavelength < 1705))
    pf, covariances = curve_fit(fit_func_civ, C4_wavelength[ix], C4_flux[ix], sigma = sigma[ix], bounds = boundaries, p0 = p0, max_nfev=1e3)
    
    ix = None
    continuum_flux = fit_func_civ(C4_wavelength, *pf)

    #converting back to linear space
    subt_C4_flux = C4_flux - continuum_flux
    
    MgII_boundaries = [[0, -2, 10, 2700, 0], [100, 2, 20, 2850, 5000]] #mu = 2798
    MgII_p0 = [10, 0, 15, 2798, 1000]
    log_wavelength, MgIIlog_flux = rebin_log(MgII_wavelength, MgII_flux)
    ix = ((log_wavelength > 2200)&(log_wavelength < 2700)) | ((log_wavelength > 2900)&(log_wavelength < 3090))
    MgII_pf, covar = curve_fit(fit_func_mgii, log_wavelength[ix], MgIIlog_flux(log_wavelength[ix]), sigma = MgII_sigma[ix], bounds = MgII_boundaries, p0 = MgII_p0, max_nfev=1e4)
    
    ix = None

    #C4_cutoffs = (log_wavelength > 1465) & (log_wavelength < 1710) 
  
    MgII_continuum_flux = fit_func_mgii(log_wavelength, *MgII_pf)
    flux_sub = MgIIlog_flux(MgII_wavelength) - MgII_continuum_flux
    lin_lam, lin_flux = rebin_lin(log_wavelength, flux_sub)

    return C4_wavelength, subt_C4_flux, pf, np.diag(covariances), lin_lam, lin_flux(lin_lam), MgII_pf, np.diag(covar)



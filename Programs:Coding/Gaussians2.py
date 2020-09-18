#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 09:36:28 2020

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit

def gauss_fit(quasar, full_wav, c4_wav, c4_flux, full_ivar, mg2_wav, mg2_flux):
    redshift = quasar[2].data['Z']
    full_wav_rf = full_wav / (1 + redshift)
    c = 3 * 10**5
    mg2_min, mg2_max, mg2_cen = 2700, 2900, 2798
    c4_min, c4_max, c4_cen = 1500, 1600, 1549


    #MgII properties:
    ivar_mg2 = full_ivar[[abs(w - full_wav_rf).argmin() for w in mg2_wav]]
    mg2_std = 1 / np.sqrt(np.abs(ivar_mg2))

    #C4 properties
    ivar_c4 = full_ivar[[abs(w - full_wav_rf).argmin() for w in c4_wav]]
    c4_std = 1 / np.sqrt(np.abs(ivar_c4))

    def gaussian(x, m, sigma, k):
        sigma_A = (sigma / c) * m  # convert km/s to delta lam in angstroms
        g = k * np.exp(-.5 * ((x - m) / sigma_A)**2)
        return g

    def gaussian3(x, m, sigma1, k1, sigma2, k2, sigma3, k3):
        gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
        return gauss

    mg2_width = 1200  # km/s
    mg2_prior = (mg2_cen, mg2_width, 1)
    mg2_bounds = [[mg2_min, 0, 0], [mg2_max, 5000, 1000]]
    mg2_fit, mg2_cov = curve_fit(gaussian, mg2_wav, mg2_flux, p0 = mg2_prior, bounds = mg2_bounds, sigma = mg2_std, max_nfev=1e3)
    plt.plot(mg2_wav, mg2_flux, 'C0-')
    plt.plot(mg2_wav, gaussian(mg2_wav, *mg2_prior), 'C1-')
    plt.plot(mg2_wav, gaussian(mg2_wav, *mg2_fit), 'C2-')

    c4_width = 5e3  # km/s
    c4_prior = (c4_cen, c4_width, 1, c4_width/2, 1, c4_width/4, 1)
    c4_bounds = [[c4_min, 0, 0, 0, 0, 0, 0], [c4_max, 5e5, 1000, 5e5, 1000, 5e5, 1000]]
    c4_fit, c4_cov = curve_fit(gaussian3, c4_wav, c4_flux, p0 = c4_prior, bounds = c4_bounds, sigma = c4_std, max_nfev=1e3) 
    plt.plot(c4_wav, c4_flux, 'C0-')
    plt.plot(c4_wav, gaussian3(c4_wav, *c4_prior), 'C1-')
    plt.plot(c4_wav, gaussian3(c4_wav, *c4_fit), 'C2-')
    return mg2_fit, np.diag(mg2_cov), c4_fit, np.diag(c4_cov)


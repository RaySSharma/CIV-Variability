#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 21:08:08 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

spectra = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7572-56944-0042.fits')

redshift = spectra[2].data['Z']
flux = (10 ** 17 * spectra[1].data['flux'])
lam = 10 ** (spectra[1].data['loglam']) / (1 + redshift)

var = 1 / spectra[1].data['ivar'] 
var2 = spectra[1].data['ivar']
sigma = np.sqrt(abs(var))

blue_flux = (lam > 1435)&(lam < 1465)
red_flux = (lam > 1690)&(lam < 1710)

bflux = lam[blue_flux]
rflux = lam[red_flux]

flux_fit = np.concatenate([flux[blue_flux], flux[red_flux]])
lam_fit = np.concatenate([lam[blue_flux], lam[red_flux]])
var_fit = np.concatenate([var2[blue_flux], var2[red_flux]])
line_fit = np.poly1d(np.polyfit(lam_fit, flux_fit, deg = 1, w = var_fit))

subtraction = flux - line_fit[flux]
print(subtraction)

#plt.plot(lam, flux, 'b')
#plt.plot(bflux, flux[blue_flux], 'r')
#plt.plot(rflux, flux[red_flux], 'r')
#line = np.linspace(1400, 1800, 100)
#plt.plot(line, line_fit(line), 'k')
#plt.xlim(1400, 1800)
#plt.xlabel('Wavelength (Angstroms)', fontsize = 10)
#plt.ylabel('Flux', fontsize = 10)
#plt.title('Spectra 7572-56944-0042', fontsize = 10)



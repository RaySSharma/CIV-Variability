#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:32:29 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

spectra = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7572-56944-0005.fits')

redshift = spectra[2].data['Z']
flux = (10 ** 17 * spectra[1].data['flux'])
lam = 10 ** (spectra[1].data['loglam']) / (1 + redshift)

var = 1 / spectra[1].data['ivar'] 
var2 = spectra[1].data['ivar']
sigma = np.sqrt(abs(var))

blue_flux = (lam > 2200)&(lam < 2700)
red_flux = (lam > 2900)&(lam < 3090)

bflux = lam[blue_flux]
rflux = lam[red_flux]

flux_fit = np.concatenate([flux[blue_flux], flux[red_flux]])
lam_fit = np.concatenate([lam[blue_flux], lam[red_flux]])
var_fit = np.concatenate([var2[blue_flux], var2[red_flux]])
line_fit = np.poly1d(np.polyfit(lam_fit, flux_fit, deg = 1, w = var_fit))

subtraction = flux - line_fit(lam)
print(subtraction)

plt.plot(lam, subtraction, 'b')
#plt.plot(bflux, flux[blue_flux], 'r')
#plt.plot(rflux, flux[red_flux], 'r')
line = np.linspace(2000, 3100, 100)
plt.plot(line, line_fit(line), 'k')
plt.xlim(2000, 3100)
plt.xlabel('Wavelength (Angstroms)', fontsize = 10)
plt.ylabel('Flux', fontsize = 10)
plt.title('Spectra 7572-56944-0005', fontsize = 10)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:28:37 2020

@author: RachelCampo
"""
import dustmaps.sfd
dustmaps.sfd.fetch()
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord as SC
from astropy import units as u
from dust_extinction.parameter_averages import F99

def dust_cor(RA, Dec, lam, f, SFD):
    wavelength = lam
    flux = f

    coordinates = SC(RA, Dec, frame = 'icrs', unit = u.deg)

    sfd = SFD
    ebv = sfd(coordinates)
    R_v = 3.1

#print('EVB for test quasar:'.format(ebv))

    ext = F99(R_v)
    flux = flux / ext.extinguish(wavelength*u.AA, ebv)
    return flux, ebv


#plt.semilogy(wavelength, flux, 'b', label = "Test Quasar's Emission")
#plt.semilogy(wavelength, flux/ext.extinguish(wavelength*u.AA, Ebv = ebv), 'y', label = 'Dereddened Quasar Emission')
#plt.xlabel('Wavelength (Angstrom)')
#plt.ylabel('Flux (ergs/cm^2/s)')
#plt.title('Dereddened Emission')
#plt.legend(loc = 'best')

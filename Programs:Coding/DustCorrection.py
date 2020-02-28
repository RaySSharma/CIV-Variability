#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:28:37 2020

@author: RachelCampo
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord as SC
from astropy import units as u
from dust_extinction.parameter_averages import F99


platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/platelist.fits')
test_quasar = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-0266-51630-0053.fits')


#Plate_RA = platelist[1].data['RACEN']
#Plate_Dec = platelist[1].data['DECCEN']
#plate = platelist[1].data['PLATE']
RA = test_quasar[2].data['PLUG_RA']
Dec = test_quasar[2].data['PLUG_DEC']
wavelength = 10**test_quasar[1].data['loglam']
flux = test_quasar[1].data['flux']

#plt.plot(wavelength, flux, 'b')

coordinates = SC(RA, Dec, frame = 'icrs', unit = u.deg)
#print(coordinates)

sfd = SFDQuery()
ebv = sfd(coordinates)
R_v = 3.1

#print('EVB for test quasar:'.format(ebv))

ext = F99(R_v)
print(ext.extinguish(wavelength*u.AA, ebv))

plt.semilogy(wavelength, flux, 'b', label = "Test Quasar's Emission")
plt.semilogy(wavelength, flux/ext.extinguish(wavelength*u.AA, Ebv = ebv), 'y', label = 'Dereddened Quasar Emission')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (ergs/cm^2/s)')
plt.title('Dereddened Emission')
plt.legend(loc = 'best')

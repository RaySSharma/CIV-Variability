#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:28:37 2020

@author: RachelCampo
"""

# This code will be the beginnings of dust corrections. I will focus on trying
# to get used to using the dust maps and subtracting them from a sample quasar

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd


platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/platelist.fits')
specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')
dustmap1 = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/SFD_dust_4096_ngp.fits')
dustmap2 = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/SFD_dust_4096_sgp.fits')
dustmap3 = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/SFD_i100_4096_ngp.fits')
dustmap4 = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/SFD_i100_4096_sgp.fits')
# should we be using SFD of should we use Schlafly & Finkbeiner?
#are there headers in the file?

Plate_RA = platelist[1].data['RACEN']
Plate_Dec = platelist[1].data['DECCEN']
plate = platelist[1].data['PLATE']

RA = specdatalist[1].data['RA']
Dec = specdatalist[1].data['DEC']
plate_number = specdatalist[1].data['PLATE']

dec_list = []
ra_list = []


#Maybe the first thing I have to do is look at both plate_ra/dec and the 
#quasar's ra/dec and figure out how much  of a difference there is?
#Possibly only have to look at the plate_ra/dec because that what was observed?

#I think I'll start with the latter. But first we have to define some functions

#coloring function
def color_function(B, V, B_0, V_0):
    dust = (B_0 - V_0) - (B - V) # this will later be plugged into E
    # also 'dust' is the dust map quantity
    #I may have to rearrange this function since the dust map give us E(B-V)
    return dust

#dust extinction
def extinction_function(F_l, F_l0):
    dust_extinction = -2.5*np.log10(F_l/F_l0)
    return dust_extinction

#extinction law
def dust_law(E, B, V):
    A_v = 3.1*E
    return A_v

for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    plateRA = Plate_RA[ix]
    ra_list.append(plateRA)
    
for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    plateDEC = Plate_Dec[ix]
    dec_list.append(plateDEC)
    

    
#print(dec_list)
#print(len(dec_list))
#print(ra_list)
#print(len(ra_list))
    
print(dustmap1[0].header)




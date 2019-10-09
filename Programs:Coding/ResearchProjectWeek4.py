#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 21:33:15 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/platelist.fits')
specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')

wavelength = specdatalist[1].data['loglam']
redshift = specdatalist[2].data['Z']
SN_Ratio = specdatalist[1].data['W1SNR']
BAL_Indicator = specdatalist[1].data['BI_CIV']
RA = platelist.data['RACEN']
Dec = platelist.data['DECCEN']
plate_quality = platelist.data['PLATEQUALITY']

#pulling the specific properties from data

data_list = np.array([])
null_list = np.array([])

for x in Dec:
    condition = specdatalist[x]
    if x > 1:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = condition
        data_list = np.concatenate([data_list, [success]])
        

for x in RA:
    condition = specdatalist[x]
    if x > 1:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = condition
        data_list = np.concatenate([data_list, [success]])


for x in BAL_Indicator:
    condition = specdatalist[x]
    if x > 0:
        failure = x
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = condition
        data_list = np.concatenate([data_list, [success]])


for x in SN_Ratio:
    condition = specdatalist[x]
    if x >= 2:
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = x
        null_list = np.concatenate([null_list, [failure]])


for x in plate_quality:
    condition = specdatalist[x]
    if x == 'good':
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])


for x in redshift:
    condition = specdatalist[x]
    if x > 1.7 and x < 2.6:
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = x
        null_list = np.concatenate([null_list, [failure]])
    
#this searches through the data and sorts out what spectra will be considered
#when graphing the data.
        
for i in data_list:
    flux = np.log10(1e17 * specdatalist[1].data['flux'])
    rest_frame_Z = 10**wavelength / (1+i)
    
    plt.plot(rest_frame_Z, flux, 'g')
    plt.ylabel('log(flux [erg/s/cm^2/Angstrom])')
    plt.xlabel('Wavelength [Angstrom]')
    
#here, I created a for loop for graphing each redshift in the data_list
#array. I adjusted the redshift to the rest frame and then it takes the rest
#frame redshift and creates a graph for each piece of data with respect to 
#flux.
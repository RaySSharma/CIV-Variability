#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 22:57:58 2019

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
plate_quality = platelist.data['PLATEQUALITY']

Plate_RA = platelist.data['RACEN']
Plate_Dec = platelist.data['DECCEN']
RA = specdatalist[1].data['RA']
Dec = specdatalist[1].data['DEC']


success_list = np.array([])
data_list = np.array([])
null_list = np.array([])

dec_list = np.array([])
ra_list = np.array([])


for x in enumerate(specdatalist):
    RAdiff = abs(Plate_RA[x] - RA[x])
    ra_list = np.concatenate([ra_list, [RAdiff]])
    
   
for x in enumerate(specdatalist):
    DECdiff = abs(Plate_Dec[x] - Dec[x])
    dec_list = np.concatenate([dec_list, [DECdiff]])

  
for i, x in enumerate(dec_list):
    condition = i
    if x > 1:
        failure = condition
        null_list = np.concatenate([null_list , [failure]])
    else:
        success = condition
        success_list = np.concatenate([success_list, [success]])
        
for i, x in enumerate(ra_list):
    condition = i
    if x > 1:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = condition
        success_list = np.concatenate([success_list, [success]])


for i, x in enumerate(BAL_Indicator):
    condition = i
    if x > 0:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = condition
        data_list = np.concatenate([data_list, [success]])


for i, x in enumerate(SN_Ratio):
    condition = i
    if x >= 2:
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])


for i, x in enumerate(plate_quality):
    condition = i
    if x == 'good':
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])


for i, x in enumerate(redshift):
    condition = i
    if x > 1.7 and x < 2.6:
        success = condition
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = condition
        null_list = np.concatenate([null_list, [failure]])
        

for x in success_list:
    if dec_list[x] <= 1 and ra_list[x] is <= 1 and BAL_Indicator[x] == 0 and SN_Ratio[x] >= 2 and plate_quality[x] == 'good' and redshift[x] > 1.7 and redshift[x] < 2.6:
        success = x
        data_list = np.concatenate([data_list, [success]]) 


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
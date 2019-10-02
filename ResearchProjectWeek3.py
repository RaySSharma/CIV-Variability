#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 07:54:17 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')

wavelength = specdatalist[1].data['loglam']
redshift = specdatalist[2].data['Z']
SN_Ratio = specdatalist[1].data['W1SNR']
BAL_Indicator = specdatalist[1].data['BI_CIV']
PlateBoreSight = specdatalist[1].data['']

#pulling the specific proerties from data

data_list = [0]
null_list = [0]

#creating a list for the data that we want, and then creating a list of
#data that we do not want to bother looking at. 

for x in PlateBoreSight:
    if x > 0.01:
        failure = x
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = x
        data_list = np.concatenate([data_list, [success]])


for x in BAL_Indicator:
    if x > 0:
        failure = x
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = x
        data_list = np.concatenate([data_list, [success]])


for x in SN_Ratio:
    if x >= 2:
        success = x
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = x
        null_list = np.concatenate([null_list, [failure]])


for x in redshift:
    if x > 1.7 and x < 2.6:
        success = x
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = x
        null_list = np.concatenate([null_list, [failure]])
    
#this searches through the data and sorts out what spectra fit in our viewing 
#window.
        
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
    
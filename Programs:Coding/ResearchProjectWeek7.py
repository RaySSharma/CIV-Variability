#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 12:44:00 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits


platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/platelist.fits')
specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')


redshift = specdatalist[1].data['Z']
SN_Ratio = platelist[1].data['PLATESN2']
BAL_Indicator = specdatalist[1].data['BI_CIV'] 
plate_quality = platelist[1].data['PLATEQUALITY']
plate_number = specdatalist[1].data['PLATE']
plate = platelist[1].data['PLATE']

Plate_RA = platelist[1].data['RACEN']
Plate_Dec = platelist[1].data['DECCEN']
RA = specdatalist[1].data['RA']
Dec = specdatalist[1].data['DEC']


success_list = np.array([])
data_list = np.array([])
null_list = np.array([])

dec_list = np.array([])
ra_list = np.array([])


for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    plateRA = Plate_RA[ix]
    ra_list = np.concatenate([ra_list, [plateRA]])
    
   
for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    plateDec = Plate_Dec[ix]
    dec_list = np.concatenate([dec_list, [plateDec]])

  
for i, x in enumerate(dec_list):
    condition = x
    if x > 1:
        failure = i
        null_list = np.concatenate([null_list , [failure]])
    else:
        success = i
        success_list = np.concatenate([success_list, [success]])
        
for i, x in enumerate(ra_list):
    condition = x
    if x > 1:
        failure = i
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = i
        success_list = np.concatenate([success_list, [success]])


for i, x in enumerate(plate_number):
    ix = plate == plate_number[i]
    SN = SN_Ratio[ix]
    if len(SN) == 1:
        if x >= 2:
            success = i
            success_list = np.concatenate([success_list, [success]])
        else:
            failure = i
            null_list = np.concatenate([null_list, [failure]])


for i, x in enumerate(plate_number):
    ix = plate == plate_number[i]
    quality = plate_quality[ix]
    if len(quality) == 1:
        if quality[0] == 'good':
            success = i
            success_list = np.concatenate([success_list, [success]])
        else:
            failure = i
            null_list = np.concatenate([null_list, [failure]])


for i, x in enumerate(BAL_Indicator):
    condition = x
    if x > 0:
        failure = i
        null_list = np.concatenate([null_list, [failure]])
    else:
        success = i
        success_list = np.concatenate([success_list, [success]])
        

for i, x in enumerate(redshift):
    condition = x
    if x > 1.7 and x < 2.6:
        success = i
        success_list = np.concatenate([success_list, [success]])
    else:
        failure = i
        null_list = np.concatenate([null_list, [failure]])
        

for x in success_list:
    if (dec_list[x] <= 1) and (ra_list[x] <= 1) and (BAL_Indicator[x] == 0) and (SN_Ratio[x] >= 2) and (plate_quality[x] == 'good') and (redshift[x] > 1.7) and (redshift[x] < 2.6):
        success = x
        data_list = np.concatenate([data_list, [success]])
    else:
        failure = x
        null_list = np.concatenate([null_list, [failure]])

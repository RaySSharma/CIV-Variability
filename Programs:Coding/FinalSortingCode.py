#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:52:52 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits

platelist = fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/plates-dr16.fits')
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


success_list = []
data_list = []
null_list = []

dec_list = []
ra_list = []


for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    RADiff = abs(RA[x] - Plate_RA[ix])
    ra_list.append(RADiff)
    
#when running RADiff, what we noticed is that some of the elements were empty
#and this could be due to possibly a quasar not having a plate number that is 
#in the platelist file. This problem may crop up later on, so make sure to keep
#an eye out for this.
   
for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    DecDiff = abs(Dec[y] - Plate_Dec[ix])
    dec_list.append(DecDiff)

  
for i, x in enumerate(dec_list):
    condition = x
    if any(x > 1):
        failure = i
        null_list.append(failure)
    else:
        success = i
        success_list.append(success)
        
for i, x in enumerate(ra_list):
    condition = x
    if any(x > 1):
        failure = i
        null_list.append(failure)
    else:
        success = i
        success_list.append(success)


for i, x in enumerate(plate_number):
    ix = plate == plate_number[i]
    SN = SN_Ratio[ix]
    if len(SN) == 1:
        if SN >= 2:
            success = i
            success_list.append(success)
        else:
            failure = i
            null_list.append(failure)


for i, x in enumerate(plate_number):
    ix = plate == plate_number[i]
    quality = plate_quality[ix]
    if len(quality) == 1:
        if quality[0] == 'good':
            success = i
            success_list.append(success)
        else:
            failure = i
            null_list.append(failure)


for i, x in enumerate(BAL_Indicator):
    condition = x
    if x > 0:
        failure = i
        null_list.append(failure)
    else:
        success = i
        success_list.append(success)
        

for i, x in enumerate(redshift):
    condition = x
    if x > 1.7 and x < 2.6:
        success = i
        success_list.append(success)
    else:
        failure = i
        null_list.append(failure)
        

nums, counts = np.unique(success_list, return_counts = True)

for nums, count in zip(nums, counts):
    if count == 6:
        data_list.append(nums)

print(len(data_list))
print(data_list)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:50:06 2019

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


success_list = []
data_list = []
null_list = []

dec_list = []
ra_list = []


for x, y in enumerate(plate_number):
    ix = plate == plate_number[x]
    RADiff = abs(RA[y] - Plate_RA[ix])
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
    if x.any() > 1:
        failure = i
        null_list.append(failure)
    else:
        success = i
        success_list.append(success)
        
for i, x in enumerate(ra_list):
    condition = x
    if x.any() > 1:
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

## Here is the data download section for this code ##

import astropy.io.fits as fits
from subprocess import call

spectra_output_dir = '/data2/rlc186/QuasarData/'  # Directory to download spectra to
spectra_url_root = 'https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/lite/'  # Spectra directory URL

quasar_list = data_list  # List of quasar indices
quasar_catalog = fits.open('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/DR14Q_v4_4.fits')  # Quasar catalog file

quasar_plate = quasar_catalog[1].data['PLATE']  # Quasar catalog plates
quasar_fiber = quasar_catalog[1].data['FIBERID']  # Quasar catalog fibers
quasar_mjd = quasar_catalog[1].data['MJD']  # Quasar catalog MJDs

for i, qso_idx in enumerate(quasar_list):
    print('Downloading spectra:', i + 1, '/', len(quasar_list))
    plate = int(quasar_plate[qso_idx])  # Corresponding quasar plate
    fiber = int(quasar_fiber[qso_idx])  # Corresponding quasar fiber
    mjd = int(quasar_mjd[qso_idx])  # Corresponding quasar MJD

    plate, mjd, fiber = str(plate), str(mjd), str(fiber).zfill(4)  # Format the three quantities

    spectra_filename = spectra_url_root + plate + '/' + 'spec-' + plate + '-' + mjd + '-' + fiber + '.fits'  # Full URL to spec


    call(['curl', '-o', spectra_output_dir + 'spec-' + plate + '-' + mjd + '-' + fiber, spectra_filename])  # Tell bash to wget the file


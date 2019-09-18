#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 22:31:11 2019

@author: RachelCampo
"""

import numpy as np
from astropy.io import fits


#this code is going to be the part where it runs through my data and figures
#out what the wavelength is after the redshift, uses the doppler shift to find
#what the wavelength was before the redshift, and then determines what 
#part of the spectrum is was from.

specdatalist = fits.open('/Users/RachelCampo/Desktop/Research/Data/spec-0266-51630-0053.fits')

wavelength = specdatalist[1].data['loglam']
redshift = specdatalist[1].data['Z']

wave_list = np.array([])
result_list = np.array([])

radio_list = ['Radio']
micro_list = ['Microwave']
infra_list = ['Infrared']
visi_list = ['Visible']
ultra_list = ['Ultraviolet']
xray_list = ['X-Ray']
gamma_list = ['Gamma Ray']
error_list = ['Error. Wavelength not found in EM spectrum.']


for x in wavelength:
    ds_eq = (10**wavelength) / (1 + redshift)
    originwav = ds_eq[0]
    wave_list = np.concatenate([wave_list, [originwav]])

#above, I took the data array and searched through it and calculated what the 
#wavelength originally looked like when first emitted by using the Doppler 
#Shift. Then after getting that result, I put it in a list called wave_list

for x in wave_list:
        
    if x >= 10**9:
        result_list = np.concatenate([result_list, [radio_list]])
    elif x < 10**9 and x >= 10**6:
       result_list = np.concatenate([result_list, [micro_list]])
    elif x < 10**6 and x >= 7000:
       result_list = np.concatenate([result_list, [infra_list]])
    elif x < 7000 and x >= 4000:
       result_list = np.concatenate([result_list, [visi_list]])
    elif x < 4000 and x >= 10:
       result_list = np.concatenate([result_list, [ultra_list]])
    elif x <10 and x >= 0.1:
       result_list = np.concatenate([result_list, [xray_list]])
    elif x < 0.1:
       result_list = np.concatenate([result_list, [gamma_list]])
    else:
       result_list = np.concatenate([result_list, [error_list]])
       
#above, I searched through each element in wave_list and looked at it's 
#wavelength and checked which part of the EM spectrum it was emitted in.
#then, for each time one of the statements were true, it would add to the list
#what corresponding wavelength it was in. The units are angstrom here.
       

print(result_list)

#the SDSS can see in a range from 3600 A to 10,400 A. The CIV line emmits in the
#1549 A range and the MgII line emmits in the 2798 A range. This is in the UV
#range.

    
    

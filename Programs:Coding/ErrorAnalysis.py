#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:53:09 2020

@author: RachelCampo
"""

#this is going to be the definition of the error analysis, the propogation of
#uncertainty

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.integrate as spi
from uncertainties import unumpy
import pandas as pd
import Analysis2
import Gaussians


#for the error, you need to calculate lum, fwhm, a, b.
#the standard deviations for these are in shen 2011, 3.7
#in Shen2011, the sd for luminosity for mgII is 0.09
#for the luminosity for cIV, 0.27

quasar = fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0105.fits')

lum, fwhm, bhm = Analysis2.analysis()


#We are also going to assume that a,b have 0 uncertainty (although unlikely)

#starting with error of lum_CIV (bisector):

lum_err_m = unumpy.uarray([lum], [0.006]) #lum will be imported from analysis
lum_err_b = unumpy.uarray([lum], [0.27])

#now with error from the Y|X lum_CIV:

lum_err_m2 = unumpy.uarray([lum], [0.009])
lum_err_b2 = unumpy.uarray([lum], [0.41])

#onto fwhm uncertainties:

fwhm_err = unumpy.uarray([fwhm], [C4_ivar]) #uses the flux variance from Gaussians2


bhm_err = unumpy.uarray([bhm], [lum_err_m + lum_err_b + fwhm_err]) #error of
#bhm with everything included (using the bisector lum_CIV error)

bhm_err2 = unumpy.uarray([bhm], [lum_err_m2 + lum_err_b2 + fwhm_err]) #error of
#bhm with everything included (using the Y|X lum_CIV)

hdr = ['Name', 'MJD', 'Fiber ID', 'Plate', 
       'Error of Luminosity Using Bisector: Slope (m)', 'Error of Luminosity Using Bisector: (b)',
       'Error of Luminosity Using Y|X: Slope (m)', 'Error of Luminosity Using Y|X: b',
       'Error of FWHM',
       'Error of Black Hole Mass: Bisector', 'Error of Black Hole Mass: Y|X']

error_list = [quasar[2].data['THING_ID'], quasar[2].data['MJD'], 
              quasar[2].data['FIBERID'], quasar[2].data['PLATE'],
              lum_err_m, lum_err_b, lum_err_m2, lum_err_b2, fwhm_err,
              bhm_err, bhm_err2]

data_frame = pd.DataFrame(error_list, columns = hdr, dtype = str)
data_frame.to_csv('error_list.csv', index = False)


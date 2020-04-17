#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 09:16:12 2020

@author: RachelCampo
"""

#this will be the master code for the entire project

#import pristine_to_mock as ptm
from astropy.io import fits
import DustCorrection
import FESubtraction
import Gaussians2
import glob

quasar_list = glob.glob('data2/rlc186/QuasarData/spec-*')
final_list = [['Name', 'MJD', 'Final Wavelength (A)', 'Final Flux', 'MgII Gaussian Fit'
               'CIV Gaussian Fit', 'EBV']]

for i in quasar_list:
    quasar = fits.open(i)
    
    wavelength = 10**quasar[1].data['loglam']
    flux = quasar[1].data['flux']
    #starting to go through each code and extracting the properties from it.
    extinguished_flux = DustCorrection.dust_cor(quasar, wavelength, flux) #returns extinguished 
    #flux and ebv
    FE_sub = FESubtraction(quasar, wavelength, extinguished_flux[0]) # returns
    #wavelength, flux
    gauss_fit = Gaussians2(quasar, FE_sub[0], FE_sub[1]) # this returns the fits
    #for both the MgII gaussian and the C4 gaussian.
    final_list.append([quasar[0].data['NAME'], quasar[0].data['MJD'], FE_sub[0],
    FE_sub[1], extinguished_flux[1]])
    
    
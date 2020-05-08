#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 09:16:12 2020

@author: RachelCampo
"""

#this will be the master code for the entire project

from astropy.io import fits
import DustCorrection
import FESubtraction
import Gaussians2
import glob

quasar_list = [fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7665-57328-0452.fits')]
#quasar_list = glob.glob('data2/rlc186/QuasarData/spec-*')
final_list = [['Name', 'MJD', 'Fiber ID', 'Plate', 'EBV', 'Fe pf Values', 
               'Fe Value of A', 'Fe Value of k', 'Fe Value of B', 'Fe Value of Mu',
               'FE Value of Sigma', 'MgII Gaussian Fit', 'MgII Mu Value', 'MgII Sigma Value',
               'MgII K Value', 'CIV Gaussian Fit', 'CIV Mu Value', 'CIV Sigma 1 Value',
               'CIV K1 Value', 'CIV Sigma 2 Value', 'CIV K2 Value', 'CIV Sigma 3 Value',
               'CIV K3 Value']]

for i in quasar_list:
    quasar = i
    
    wavelength = 10**quasar[1].data['loglam']
    flux = quasar[1].data['flux']
    var = quasar[1].data['ivar']
    #starting to go through each code and extracting the properties from it.
    extinguished_flux, ebv = DustCorrection.dust_cor(quasar, wavelength, flux) #returns extinguished 
    #flux and ebv
    FE_wav, FE_flux, pf, pcov = FESubtraction.FE_sub(quasar, wavelength, extinguished_flux) # returns
    #wavelength, flux, pf, diagonal pcov
    mg2gauss, mg2pcov, c4gauss, c4pcov  = Gaussians2.gauss_fit(quasar, FE_wav, FE_flux, var) # this returns the fits
    #for both the MgII gaussian and the C4 gaussian along with respective diagonal pcov.
    
    A = pcov[0]
    k = pcov[1]
    B = pcov[2]
    mu = pcov[3]
    sigma = pcov[4]
    #pulling out each number from the diagonal from Fe pcov
    
    MgII_mu = mg2pcov[0]
    MgII_sigma = mg2pcov[1]
    MgII_k = mg2pcov[2]
    #pulling out each number from MgII pcov
    
    CIV_mu = c4pcov[0]
    CIV_sig1 = c4pcov[1]
    CIV_k1 = c4pcov[2]
    CIV_sig2 = c4pcov[3]
    CIV_k2 = c4pcov[4]
    CIV_sig3 = c4pcov[5]
    CIV_k3 = c4pcov[6]
    
    final_list.append([quasar[0].data['THING_ID'], quasar[0].data['MJD'], 
                       quasar[0].data['FIBER_ID'], quasar[0].data['PLATE'], ebv, 
                       pf, A, k, B, mu, sigma, mg2gauss, MgII_mu, MgII_sigma, MgII_k, 
                       *c4gauss, CIV_mu, CIV_sig1, CIV_k1, CIV_sig2, CIV_k2, 
                       CIV_sig3, CIV_k3])
    
  # take diagonals of the covariance matricies and add these into the list.
print(final_list)
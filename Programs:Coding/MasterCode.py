#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 09:16:12 2020

@author: RachelCampo
"""

#this will be the master code for the entire project

from astropy.io import fits
import numpy as np
import DustCorrection
import FESubtraction
import Gaussians2
import glob

quasar_list = [fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7665-57328-0452.fits'), fits.open('/Users/RachelCampo/Desktop/Research/Data/Other Spectra/spec-7128-56567-0284.fits')]
#quasar_list = glob.glob('data2/rlc186/QuasarData/spec-*')
final_list = [['Name', 'MJD', 'Fiber ID', 'Plate', 
               'EBV', 
               'CIV pf Values from FE Subtraction', 'MgII pf Values from FE Subtraction'
               'CIV Value of A from FE Subtraction', 'CIV Value of k from FE Subtraction', 
               'CIV Value of B from FE Subtraction', 'CIV Value of Mu from FE Subtraction',
               'CIV Value of Sigma from FE Subtraction', 'MgII Value of A from FE Subtraction',
               'MgII Value of k from FE Subtraction', 'MgII Value of B from FE Subtraction',
               'MgII Value of Mu from FE Subtraction', 'MgII Value of Sigma from Fe Subtraction',
               'MgII Gaussian Fit', 'MgII Mu Value from Gaussian Fitting', 
               'MgII Sigma Value from Gaussian Fitting',
               'MgII K Value from Gaussian Fitting', 'CIV Gaussian Fit', 
               'CIV Mu Value from Gaussian Fitting', 'CIV Sigma 1 Value from Gaussian Fitting',
               'CIV K1 Value from Gaussian Fitting', 'CIV Sigma 2 Value from Gaussian Fitting', 
               'CIV K2 Value from Gaussian Fitting', 'CIV Sigma 3 Value from Gaussian Fitting',
               'CIV K3 Value from Gaussian Fitting']]

for i in quasar_list:
    
    quasar = i
    
    wavelength = 10**quasar[1].data['loglam']
    flux = quasar[1].data['flux']
    var = quasar[1].data['ivar']
    #starting to go through each code and extracting the properties from it.
    extinguished_flux, ebv = DustCorrection.dust_cor(quasar, wavelength, flux) #returns extinguished 
    #flux and ebv
    C4_wav, C4_flux, C4pf, C4pcov, MgII_wav, MgII_flux, MgIIpf, MgIIpcov = FESubtraction.FE_sub(quasar, wavelength, extinguished_flux) # returns
    #wavelength, flux, pf, diagonal pcov
    mg2gauss, mg2pcov, c4gauss, c4pcov  = Gaussians2.gauss_fit(quasar, extinguished_flux, C4_wav, C4_flux, var, MgII_wav, MgII_flux) # this returns the fits
    #for both the MgII gaussian and the C4 gaussian along with respective diagonal pcov.
    
    A = C4pcov[0]
    k = C4pcov[1]
    B = C4pcov[2]
    mu = C4pcov[3]
    sigma = C4pcov[4]
    #pulling out each number from the diagonal from Fe pcov
    
    A_mg = MgIIpcov[0]
    k_mg = MgIIpcov[1]
    B_mg = MgIIpcov[2]
    mu_mg = MgIIpcov[3]
    sigma_mg = MgIIpcov[4]
    #pulling out covariances of mgII from fe_sub
    
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
    
    final_list.append([quasar[2].data['THING_ID'], quasar[2].data['MJD'], 
                       quasar[2].data['FIBERID'], quasar[2].data['PLATE'], ebv, 
                       C4pf, MgIIpf, A, k, B, mu, sigma, A_mg, k_mg, B_mg, mu_mg,
                       sigma_mg, mg2gauss, MgII_mu, MgII_sigma, MgII_k, 
                       c4gauss, CIV_mu, CIV_sig1, CIV_k1, CIV_sig2, CIV_k2, 
                       CIV_sig3, CIV_k3])
    
  # take diagonals of the covariance matricies and add these into the list.
np.savetxt('final_list.csv', final_list, fmt = '%s')
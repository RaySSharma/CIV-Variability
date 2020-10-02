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
import pandas as pd
import tqdm

spectra_dir = '/data2/rlc186/QuasarData/'
spectra = glob.glob('/data2/rlc186/QuasarData/spec-*')
final_list = []

hdr = ['Name', 'MJD', 'Fiber ID', 'Plate', 'Redshift',
               'EBV',
               'CIV Value of A from FE Subtraction', 'Error of CIV A', 'CIV Value of k from FE Subtraction', 
               'Error of CIV K', 'MgII Value of A from FE Subtraction',
               'Error of MgII A', 'MgII Value of k from FE Subtraction', 'Error of MgII k', 'MgII Value of B from FE Subtraction', 'Error of MgII B',
               'MgII Value of Mu from FE Subtraction', 'Error of MgII Mu', 'MgII Value of Sigma from Fe Subtraction', 'Error of MgII Sigma',
               'MgII Mu Value from Gaussian Fitting', 'Error of MgII Mu from Gaussian Fitting',
               'MgII Sigma Value from Gaussian Fitting', 'Error of MgII Sigma from Gaussian Fitting',
               'MgII K Value from Gaussian Fitting', 'Error of MgII K from Gaussian Fitting',
               'CIV Mu Value from Gaussian Fitting', 'Error of CIV Mu from Gaussian Fitting',
               'CIV Sigma 1 Value from Gaussian Fitting', 'Error of CIV Sigma 1 from Gaussian Fitting',
               'CIV K1 Value from Gaussian Fitting', 'Error of CIV K1 from Gaussian Fitting',
               'CIV Sigma 2 Value from Gaussian Fitting', 'Error of CIV Sigma 2 from Gaussian Fitting',
               'CIV K2 Value from Gaussian Fitting', 'Error of CIV K2 from Gaussian Fitting',
               'CIV Sigma 3 Value from Gaussian Fitting', 'Error of CIV Sigma 3 from Gaussian Fitting',
               'CIV K3 Value from Gaussian Fitting', 'Error of CIV K3 from Gaussian Fitting']

def get_safe_pixels(mask, ivar):
    import numpy as np
    NOPLUG = mask & 2**0
    BADTRACE = mask & 2**1
    BADFLAT = mask & 2**2
    BADARC = mask & 2**3
    MANYBADCOLUMNS = mask & 2**4
    MANYREJECTED = mask & 2**5
    NEARBADPIXEL = mask & 2**16
    LOWFLAT = mask & 2**17
    FULLREJECT = mask & 2**18
    SCATTEREDLIGHT = mask & 2**20
    NOSKY = mask & 2**22
    BRIGHTSKY = mask & 2**23
    COMBINERJ = mask & 2**25
    REDMONSTER = mask & 2**28
    
    block = (NOPLUG + BADTRACE + BADFLAT + BADARC + MANYBADCOLUMNS + MANYREJECTED + NEARBADPIXEL + LOWFLAT + FULLREJECT + SCATTEREDLIGHT + NOSKY + BRIGHTSKY + COMBINERJ + REDMONSTER)

    safe_pixels = (block == 0) & (ivar > 0)
    return safe_pixels

for i in spectra:
    
    quasar = fits.open(i, ignore_missing_end = True)
    
    safe_pixels = get_safe_pixels(quasar[1].data['and_mask'], quasar[1].data['ivar'])
    wavelength = 10**quasar[1].data['loglam'][safe_pixels]
    flux = quasar[1].data['flux'][safe_pixels]
    ivar = quasar[1].data['ivar'][safe_pixels]
    z = quasar[2].data['Z']
    #starting to go through each code and extracting the properties from it.
    try:
        extinguished_flux, ebv = DustCorrection.dust_cor(quasar, wavelength, flux) #returns extinguished 
    #flux and ebv
    except:
        print('Quasar ' + str(quasar[2].data['THING_ID'] + quasar[2].data['MJD'] +
              quasar[2].data['FIBERID'] + quasar[2].data['PLATE']) + 
              ' has failed the Dust Correction Code')
        continue
    
    try:
        C4_wav, C4_flux, C4pf, C4pcov, MgII_wav, MgII_flux, MgIIpf, MgIIpcov = FESubtraction.FE_sub(wavelength, extinguished_flux, z, ivar) # returns
    #wavelength, flux, pf, diagonal pcov
    except:
        print('Quasar ' + str(quasar[2].data['THING_ID'] + quasar[2].data['MJD'] +
              quasar[2].data['FIBERID'] + quasar[2].data['PLATE']) + 
              ' has failed the Iron Subtraction Code')
        continue
        
    try:    
        mg2gauss, mg2pcov, c4gauss, c4pcov  = Gaussians2.gauss_fit(quasar, wavelength, C4_wav, C4_flux, ivar, MgII_wav, MgII_flux) # this returns the fits
    #for both the MgII gaussian and the C4 gaussian along with respective diagonal pcov.
    except:
        print('Quasar ' + str(quasar[2].data['THING_ID'] + quasar[2].data['MJD'] +
              quasar[2].data['FIBERID'] + quasar[2].data['PLATE']) + 
              ' has failed the Gaussian Fitting Code')
        continue
    
    A = C4pf[0]
    A_err = C4pcov[0]
    k = C4pf[1]
    k_err = C4pcov[1]
    #pulling out each number from the diagonal from Fe pcov
    #make sure to change to the pf values, not the pcov values!! pcov are the variances
    
    A_mg = MgIIpf[0]
    A_mg_err = MgIIpcov[0]
    k_mg = MgIIpf[1]
    k_mg_err = MgIIpcov[1]
    B_mg = MgIIpf[2]
    B_mg_err = MgIIpcov[2]
    mu_mg = MgIIpf[3]
    mu_mg_err = MgIIpcov[3]
    sigma_mg = MgIIpf[4]
    sigma_mg_err = MgIIpcov[4]
    #pulling out covariances of mgII from fe_sub
    
    MgII_mu = mg2gauss[0]
    MgII_mu_err = mg2pcov[0]
    MgII_sigma = mg2gauss[1]
    MgII_sigma_err = mg2pcov[1]
    MgII_k = mg2gauss[2]
    MgII_k_err = mg2pcov[2]
    #pulling out each number from MgII pcov
    
    CIV_mu = c4gauss[0]
    CIV_mu_err = c4pcov[0]
    CIV_sig1 = c4gauss[1]
    CIV_sig1_err = c4pcov[1]
    CIV_k1 = c4gauss[2]
    CIV_k1_err = c4pcov[2]
    CIV_sig2 = c4gauss[3]
    CIV_sig2_err = c4pcov[3]
    CIV_k2 = c4gauss[4]
    CIV_k2_err = c4pcov[4]
    CIV_sig3 = c4gauss[5]
    CIV_sig3_err = c4pcov[5]
    CIV_k3 = c4gauss[6]
    CIV_k3_err = c4pcov[6]
    
    final_list.append([quasar[2].data['THING_ID'], quasar[2].data['MJD'], 
                       quasar[2].data['FIBERID'], quasar[2].data['PLATE'], z, ebv, 
                       A, A_err, k, k_err, 
                       A_mg, A_mg_err, k_mg, k_mg_err, B_mg, B_mg_err, mu_mg, mu_mg_err,
                       sigma_mg, sigma_mg_err, MgII_mu, MgII_mu_err, MgII_sigma, 
                       MgII_sigma_err, MgII_k, MgII_k_err,
                       CIV_mu, CIV_mu_err, CIV_sig1, CIV_sig1_err, CIV_k1, CIV_k1_err,
                       CIV_sig2, CIV_sig2_err, CIV_k2, CIV_k2_err,
                       CIV_sig3, CIV_sig3_err, CIV_k3, CIV_k3_err])

    
    quasar.close()


data_frame = pd.DataFrame(final_list, columns = hdr, dtype = float)

data_frame.to_csv('final_list.csv', index = False)

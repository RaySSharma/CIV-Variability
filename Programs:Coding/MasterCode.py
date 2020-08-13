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

quasar_list = [fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0180.fits'), fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0105.fits')]
                #its.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0106.fits'), fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0109.fits'), 
               # fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0112.fits'), fits.open('/Users/RachelCampo/Desktop/Research/Data/Official Data/spec-5202-55824-0152.fits')]
#quasar_list = glob.glob('data2/rlc186/QuasarData/spec-*')
final_list = []
hdr = ['Name', 'MJD', 'Fiber ID', 'Plate', 'Redshift',
               'EBV',
               'CIV Value of A from FE Subtraction', 'Error of CIV A', 'CIV Value of k from FE Subtraction', 
               'Error of CIV K', 'CIV Value of B from FE Subtraction', 'Error of CIV B', 'CIV Value of Mu from FE Subtraction',
               'Error of CIV Mu', 'CIV Value of Sigma from FE Subtraction', 'Error of CIV Sigma', 'MgII Value of A from FE Subtraction',
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

for i in quasar_list:
    
    quasar = i
    
    wavelength = 10**quasar[1].data['loglam']
    flux = quasar[1].data['flux']
    var = quasar[1].data['ivar']
    z = quasar[2].data['Z']
    #starting to go through each code and extracting the properties from it.
    extinguished_flux, ebv = DustCorrection.dust_cor(quasar, wavelength, flux) #returns extinguished 
    #flux and ebv
    C4_wav, C4_flux, C4pf, C4pcov, MgII_wav, MgII_flux, MgIIpf, MgIIpcov = FESubtraction.FE_sub(quasar, wavelength, extinguished_flux) # returns
    #wavelength, flux, pf, diagonal pcov
    mg2gauss, mg2pcov, c4gauss, c4pcov  = Gaussians2.gauss_fit(quasar, extinguished_flux, C4_wav, C4_flux, var, MgII_wav, MgII_flux) # this returns the fits
    #for both the MgII gaussian and the C4 gaussian along with respective diagonal pcov.
    
    A = C4pf[0]
    A_err = C4pcov[0]
    k = C4pf[1]
    k_err = C4pcov[1]
    B = C4pf[2]
    B_err = C4pcov[2]
    mu = C4pf[3]
    mu_err = C4pcov[3]
    sigma = C4pf[4]
    sigma_err = C4pcov[4]
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
                       A, A_err, k, k_err, B, B_err, mu, mu_err, sigma, sigma_err, 
                       A_mg, A_mg_err, k_mg, k_mg_err, B_mg, B_mg_err, mu_mg, mu_mg_err,
                       sigma_mg, sigma_mg_err, MgII_mu, MgII_mu_err, MgII_sigma, 
                       MgII_sigma_err, MgII_k, MgII_k_err,
                       CIV_mu, CIV_mu_err, CIV_sig1, CIV_sig1_err, CIV_k1, CIV_k1_err,
                       CIV_sig2, CIV_sig2_err, CIV_k2, CIV_k2_err,
                       CIV_sig3, CIV_sig3_err, CIV_k3, CIV_k3_err])
    #data_frame = pd.DataFrame({hdr[0]: [final_list[0][0]], hdr[1]: [final_list[0][1]],
#                               hdr[2]: [final_list[0][2]], hdr[3]: [final_list[0][3]],
#                               hdr[4]: [final_list[0][4]], hdr[5]: [final_list[0][5]],
#                               hdr[6]: [final_list[0][6]], hdr[7]: [final_list[0][7]],
#                               hdr[8]: [final_list[0][8]], hdr[9]: [final_list[0][9]],
#                               hdr[10]: [final_list[0][10]], hdr[11]: [final_list[0][11]],
#                               hdr[12]: [final_list[0][12]], hdr[13]: [final_list[0][13]],
#                               hdr[14]: [final_list[0][14]], hdr[15]: [final_list[0][15]],
#                               hdr[16]: [final_list[0][16]], hdr[17]: [final_list[0][17]],
#                               hdr[18]: [final_list[0][18]], hdr[19]: [final_list[0][19]],
#                               hdr[20]: [final_list[0][20]], hdr[21]: [final_list[0][21]],
#                               hdr[22]: [final_list[0][22]], hdr[23]: [final_list[0][23]],
#                               hdr[24]: [final_list[0][24]]})
    
    data_frame = pd.DataFrame(final_list, columns = hdr, dtype = float)


data_frame.to_csv('final_list.csv', index = False)
#comp = dict(method = 'zip', archive_name = 'final_list.csv')
#data_frame.to_csv('final_list.zip', index = False, compression = comp)
#np.savetxt('final_list.csv', data_frame, fmt = '%s', delimiter = ',')

#make sure you close each quasar
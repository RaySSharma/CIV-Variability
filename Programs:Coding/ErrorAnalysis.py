#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:53:09 2020

@author: RachelCampo
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.integrate as spi
from uncertainties import unumpy
import uncertainties
import pandas as pd
import Analysis2
import Gaussians
from astropy import units as u
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15 as p15
import astropy.constants as const
from scipy.interpolate import UnivariateSpline
import pdb

#for the error, you need to calculate lum, fwhm, a, b.
c = 3 * 10 ** 5
lum, fwhm, bhm = Analysis2.analysis()
final_list = pd.read_csv('/Users/RachelCampo/Desktop/Research/CIV-Variability/Programs/final_list.csv', sep = ',', index_col = False)
Q = len(final_list)

mu = final_list.loc[:, 'CIV Mu Value from Gaussian Fitting'].values
mu_err = final_list.loc[:, 'Error of CIV Mu from Gaussian Fitting'].values
sig1 = final_list.loc[:, 'CIV Sigma 1 Value from Gaussian Fitting'].values
sig1_err = final_list.loc[:, 'Error of CIV Sigma 1 from Gaussian Fitting'].values
c4k1 = final_list.loc[:, 'CIV K1 Value from Gaussian Fitting'].values
c4k1_err = final_list.loc[:, 'Error of CIV K1 from Gaussian Fitting'].values
sig2 = final_list.loc[:, 'CIV Sigma 2 Value from Gaussian Fitting'].values
sig2_err = final_list.loc[:, 'Error of CIV Sigma 2 from Gaussian Fitting'].values
c4k2 = final_list.loc[:, 'CIV K2 Value from Gaussian Fitting'].values
c4k2_err = final_list.loc[:, 'Error of CIV K2 from Gaussian Fitting'].values
sig3 = final_list.loc[:, 'CIV Sigma 3 Value from Gaussian Fitting'].values
sig3_err = final_list.loc[:, 'Error of CIV Sigma 3 from Gaussian Fitting'].values
c4k3 = final_list.loc[:, 'CIV K3 Value from Gaussian Fitting'].values
c4k3_err = final_list.loc[:, 'Error of CIV K3 from Gaussian Fitting'].values
redshift = final_list.loc[:, 'Redshift'].values

mu_mg = final_list.loc[:, 'MgII Mu Value from Gaussian Fitting'].values
mu_mg_er = final_list.loc[:, 'Error of MgII Mu from Gaussian Fitting'].values
mg_sig = final_list.loc[:, 'MgII Sigma Value from Gaussian Fitting'].values
mg_sig_er = final_list.loc[:, 'Error of MgII Sigma from Gaussian Fitting'].values
mg_k = final_list.loc[:, 'MgII K Value from Gaussian Fitting'].values
mg_k_er = final_list.loc[:, 'Error of MgII K from Gaussian Fitting'].values

name = (final_list.loc[:, 'Name'].values)
mjd = (final_list.loc[:, 'MJD'].values)
fiberid = (final_list.loc[:, 'Fiber ID'].values)
plate = (final_list.loc[:, 'Plate'].values)

#We are also going to assume that a,b have 0 uncertainty (although unlikely)

def gaussian(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k.reshape(Q, 1) * unumpy.exp(-.5 * ((x * np.ones((Q, len(x))) - m.reshape(Q, 1)) / sigma.reshape(Q,1))**2)
        return g

def gaussian3(x, unp_array):
        m, sigma1, k1, sigma2, k2, sigma3, k3 = unp_array
        gauss = gaussian(x, m, sigma1, k1) + gaussian(x, m, sigma2, k2) + gaussian(x, m, sigma3, k3)
        return gauss

def Gaussian(x, uarray):
        m, sigma, k = uarray
        sigma = (sigma / c) * m
        g = k.reshape(Q, 1) * unumpy.exp(-.5 * ((x * np.ones((Q, len(x))) - m.reshape(Q, 1)) / sigma.reshape(Q,1))**2)
        return g

C4_wav = np.linspace(1500, 1600, 1000)
Mg_wav = np.linspace(2750, 2850, 1000)
unp_array = unumpy.uarray([mu, sig1, c4k1, sig2, c4k2, sig3, c4k3],
                         [mu_err, sig1_err, c4k1_err, sig2_err, c4k2_err, sig3_err, c4k3_err])
uarray = unumpy.uarray([mu_mg, mg_sig, mg_k], [mu_mg_er, mg_sig_er, mg_k_er])

C4_flux = gaussian3(C4_wav, unp_array)
Mg_flux = Gaussian(Mg_wav, uarray)
d = p15.luminosity_distance(redshift).to('cm')
C4_flux = C4_flux * 10 ** -17
mg_flux = Mg_flux * 10 ** -17
#C4_wav = C4_wav * u.Angstrom

def FWHM(x, y):
    y_new = y - y.min()
    spline = UnivariateSpline(x, y_new - y_new.max() / 2, s = 0)
    try:
        a, b = spline.roots()
    except:
        return np.nan
    return abs(a - b)

def A_to_kms(fwhm, m):
    return const.c * fwhm / m

C4_luminosity = (C4_flux.T * (4 * np.pi * d ** 2)).T #ergs/s/A
Mg_luminosity = (mg_flux.T * (4 * np.pi * d ** 2)).T

def C4_lum(wav, lum):
    return np.trapz(lum, wav)

C4_L = C4_lum(C4_wav, C4_luminosity) #luminosity error
Mg_L = C4_lum(Mg_wav, Mg_luminosity)

#doing FWHM manually:
mean_list = []
std_list = []

mg_mean_list = []
mg_std_list = []

for j in range(len(final_list)):
    mu_new = np.random.normal(loc = mu[j], scale = mu_err[j], size = 1000)
    sig1_new = np.random.normal(loc = sig1[j], scale = sig1_err[j], size = 1000)
    sig2_new = np.random.normal(loc = sig2[j], scale = sig2_err[j], size = 1000)
    sig3_new = np.random.normal(loc = sig3[j], scale = sig3_err[j], size = 1000)
    c4k1_new = np.random.normal(loc = c4k1[j], scale = c4k1_err[j], size = 1000)
    c4k2_new = np.random.normal(loc = c4k2[j], scale = c4k2_err[j], size = 1000)
    c4k3_new = np.random.normal(loc = c4k3[j], scale = c4k3_err[j], size = 1000)
    
    mg_mu_new = np.random.normal(loc = mu_mg[j], scale = mu_mg_er[j], size = 1000)
    mg_sig_new = np.random.normal(loc = mg_sig[j], scale = mg_sig_er[j], size = 1000)
    mg_k_new = np.random.normal(loc = mg_k[j], scale = mg_k_er[j], size = 1000)

    def gauss(x, m, sigma, k):
        sigma = (sigma / c) * m
        g = k * np.exp(-.5 * ((x - m) / sigma)**2)
        return g

    def gauss3(x, mu_new, sig1_new, sig2_new, sig3_new, c4k1_new, c4k2_new, c4k3_new):
        g = gauss(x, mu_new, sig1_new, c4k1_new) + gauss(x, mu_new, sig2_new, c4k2_new) + gauss(x, mu_new, sig3_new, c4k3_new)
        return g

    flux_new = [gauss3(C4_wav, mu_new[i], sig1_new[i], sig2_new[i], sig3_new[i],
                c4k1_new[i], c4k2_new[i], c4k3_new[i]) for i in range(1000)]
    mg_flux_new = [gauss(Mg_wav, mg_mu_new[u], mg_sig_new[u], mg_k_new[u]) for u in range(1000)]

    fwhm_new = [FWHM(C4_wav, x) for x in flux_new] #fwhm err
    fwhm_mg = [FWHM(Mg_wav, y) for y in mg_flux_new]
    fwhm_mg_std = np.nanstd(fwhm_mg)
    fwhm_mg_mean = np.nanmean(fwhm_mg)
    fwhm_err_std = np.nanstd(fwhm_new)
    fwhm_err_mean = np.nanmean(fwhm_new)

   
    fwhm_final_mean = A_to_kms(fwhm_err_mean * u.Angstrom, mu[j] * u.Angstrom).to('km/s').value
    fwhm_final_std = A_to_kms(fwhm_err_std * u.Angstrom, mu[j] * u.Angstrom).to('km/s').value #you want to make a list of this for each quasar
    #convert these into a unumpy for uncertainties
    
    fwhm_final_mean_mg = A_to_kms(fwhm_mg_mean * u.Angstrom, mu_mg[j] * u.Angstrom).to('km/s').value
    fwhm_final_std_mg = A_to_kms(fwhm_mg_std * u.Angstrom, mu_mg[j] * u.Angstrom).to('km/s').value
    
    mean_list.append(fwhm_final_mean)
    std_list.append(fwhm_final_std)
    
    mg_mean_list.append(fwhm_final_mean_mg)
    mg_std_list.append(fwhm_final_std_mg)

fwhm_unumpy = unumpy.uarray(mean_list, std_list)
fwhm_unumpy_mg = unumpy.uarray(mg_mean_list, mg_std_list)


#finding black hole mass error:
def mass_bh(lum, fwhm, a = 0.660, b = 0.53):
    return 10 ** (a + b * unumpy.log10(lum / 1e44)
                  + 2 * unumpy.log10(fwhm))

bhm_err = mass_bh(C4_L.value, fwhm_unumpy) #bhm_err
bhm_err_mg = mass_bh(Mg_L.value, fwhm_unumpy_mg)
print(bhm_err_mg)

bhm_n = unumpy.nominal_values(bhm_err)
bhm_std = unumpy.std_devs(bhm_err)

C4_L_values = unumpy.nominal_values(C4_L.value)
C4_L_std = unumpy.std_devs(C4_L.value)

fwhm_unp = unumpy.nominal_values(fwhm_unumpy)
fwhm_std = unumpy.std_devs(fwhm_unumpy)

line_shift = (1550 - mu)
mg_line_shift = (2200 - mu_mg)

fig, ax = plt.subplots()

ax.errorbar(bhm_n, C4_L_values, xerr = bhm_std, yerr = C4_L_std, fmt = 'o')
ax.set_ylabel('Luminosity w/ Error (erg/s)')
ax.set_xlabel('Black Hole Mass w/ Error (Solar Masses)')
ax.set_title('Luminiosity vs. Black Hole Mass')
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()

#
#hdr = ['Name', 'MJD', 'Fiber ID', 'Plate']
#       'Black Hole Mass Using CIV (Solar Mass)', 'Black Hole Mass Using CIV Error',
#       'Luminosity Using CIV (ergs/s)', 'Luminosity Using CIV Error',
#       'Full Width Half Max Using CIV (km/s)', 'Full Width Half Max Using CIV Error']
#       'Black Hole Mass Using MgII (Solar Mass)', 'Black Hole Mass Using MgII Error',
#       'Luminosity Using MgII (erg/s)', 'Luminosity Using MgII Error',
#       'Full Width Half Max Using MgII (km/s)', 'Full Width Half Max Using MgII Error'
#       'Line Shift of CIV (km/s)', 'Line Shift of MgII (km/s)']
#
#error_list = [name, mjd, fiberid, plate]
#              bhm_n, bhm_std,
#              C4_L_values, C4_L_std,
#              fwhm_unp, fwhm_std]
#
#data_frame = pd.DataFrame(data = error_list, columns = hdr, dtype = str)
#data_frame.to_csv('error_list.csv', index = False)

#make sure to track the units!
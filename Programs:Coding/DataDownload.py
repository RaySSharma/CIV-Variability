#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 17:33:51 2019

@author: RachelCampo
"""

import astropy.io.fits as fits
from subprocess import call

spectra_output_dir = '/data/rlc186/QuasarData/'  # Directory to download spectra to
spectra_url_root = 'https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/lite/'  # Spectra directory URL

quasar_list = [0, 1, 2, 3, 4, 5]  # List of quasar indices
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


    call(['curl', '-O', spectra_output_dir, spectra_filename])  # Tell bash to wget the file
import astropy.io.fits as fits
from subprocess import call

spectra_output_dir = '/data/ramonsharma/spectra/'  # Directory to download spectra to
spectra_url_root = 'https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/lite/'  # Spectra directory URL

quasar_list = [0, 1, 2, 3, 4, 5]  # List of quasar indices
quasar_catalog = fits.open(***)  # Quasar catalog file

quasar_plate = quasar_catalog[1].data['PLATE']  # Quasar catalog plates
quasar_fiber = quasar_catalog[1].data['FIBERID']  # Quasar catalog fibers
quasar_mjd = quasar_catalog[1].data['MJD']  # Quasar catalog MJDs

for i, qso_idx in enumerate(quasar_list):
    print('Downloading spectra:', i + 1, '/', len(quasar_list))
    plate = int(quasar_plate[qso_idx])  # Corresponding quasar plate
    fiber = int(***)  # Corresponding quasar fiber
    mjd = int(***)  # Corresponding quasar MJD

    plate, mjd, fiber = str(plate), str(mjd), str(fiber).zfill(4)  # Format the three quantities

    spectra_filename = spectra_url_root + plate + '/' + 'spec-' + plate + '-' + mjd + '-' + fiber + '.fits'  # Full URL to spec


    call(['wget', '-P', spectra_output_dir, spectra_filename])  # Tell bash to wget the file
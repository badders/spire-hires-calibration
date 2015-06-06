from common import *

hdulist = fits.open('calibration/hires_spectra.fits')

params = fitgaussian2D(hdulist[1].data)
print(sqrt(params))

from pylab import *
from astropy.io import fits
from astropy.units import arcmin, degree
normal = 'plots/m74_PLW_10_0.fits'
hires = 'plots/m74_PLW_10_0_hires.fits'
truth = 'plots/TRUTH_m74_PLW.fits'

normal_fits = fits.open(normal)
normal_data = normal_fits[1].data
normal_data *= zeros_like(normal_data)
normal_header = normal_fits[1].header
normal_pix_size = abs(normal_header['CDELT1'])

hires_fits = fits.open(hires)
hires_data = hires_fits[1].data
hires_data *= zeros_like(hires_data)
hires_header = hires_fits[1].header
hires_pix_size = abs(hires_header['CDELT1'])

truth_fits = fits.open(truth)
truth_data = truth_fits[1].data
truth_data *= zeros_like(truth_data)
truth_header = truth_fits[1].header
truth_pix_size = abs(truth_header['CDELT1'])

def gauss2D(x, y, cx, cy, width, amplitude=1):
    return amplitude * e**(-((x - cx)**2 / (2*width**2) + (y - cy)**2 / (2*width**2)))

# Make gaussian

# Arcminutes
width = arcmin / (normal_pix_size * degree)
w, h = normal_data.shape
cx = w / 2
cy = h / 2

print(normal_pix_size)
print(degree)
print(arcmin)
print('Creating normal gaussian ...')
print('1 arcminute = ', width,' pixels')

for i in arange(w):
    for j in arange(h):
        normal_data[i][j] = gauss2D(i, j, cx, cy, width)

normal_fits.writeto('calibration/normal.fits', clobber=True)

width = arcmin / (hires_pix_size * degree)
w, h = hires_data.shape
cx = w / 2
cy = h / 2

print('Creating hires gaussian ...')
print('1 arcminute = ', width,' pixels')

for i in arange(w):
    for j in arange(h):
        hires_data[i][j] = gauss2D(i, j, cx, cy, width)

hires_fits.writeto('calibration/hires.fits', clobber=True)

width = arcmin / (truth_pix_size * degree)
w, h = truth_data.shape
cx = w / 2
cy = h / 2

print('Creating truth gaussian ...')
print('1 arcminute = ', width,' pixels')

for i in arange(w):
    for j in arange(h):
        truth_data[i][j] = gauss2D(i, j, cx, cy, width)

truth_fits.writeto('calibration/truth.fits', clobber=True)

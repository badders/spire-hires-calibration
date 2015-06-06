from common import *

factors = [0, 0.1, 0.2, 0.5, 0.8, 1, 1.3, 1.6,  2, 2.8, 3.5, 5, 7.5, 10, 20]

snr = []
beam_size = []

for f in factors:
    hdulist2 = fits.open('plots/m74_PLW_{}_0_beam.fits'.format(f))
    data = hdulist2[1].data
    w, h = data.shape
    cx, cy = w / 2, h / 2
    snr.append(SNR('plots/m74_PLW_{}_0.fits'.format(f)))
    params = fitgaussian2D(data[cy-10:cy+10,cx-10:cx+10])
    beam_size.append(params[3])
    if f==0 and False:
        imshow(data[cy-10:cy+10,cx-10:cx+10])
        figure()

plot(snr, beam_size)
print(beam_size)

show()

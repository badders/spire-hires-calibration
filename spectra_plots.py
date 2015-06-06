from pylab import *
import glob
from astropy.io import fits

spectra_files = glob.glob('plots/*spectra*fits')

for f in spectra_files:
    hdulist = fits.open(f)
    data = hdulist[1].data
    w, h = data.shape
    cx, cy = w // 2, h //2
    rs = zeros(w * h)
    vs = zeros(w * h)

    c = 0
    for i in range(0, w):
        for j in range(0, h):
            r = sqrt((cx - i)**2 + (cy - j)**2)
            v = data[i][j]
            rs[c] = r
            vs[c] = v
            c +=1

    binned_vs = zeros(int(rs.max() + 1))
    bin_counts = zeros_like(binned_vs)

    for i in range(len(rs)):
        r = rs[i]
        v = vs[i]
        binned_vs[int(r)] += v
        bin_counts[int(r)] += 1

    binned_vs = binned_vs / bin_counts

    shrunk_vs = zeros(int(len(binned_vs) / 2))

    shrunk_vs = binned_vs[::2]
    i = 0
    for j in binned_vs[1::2]:
        shrunk_vs[i] += j
        i += 1

    shrunk_vs /= 2
    shrunk_bins = arange(int(ceil(len(binned_vs) / 2)), step=2)

    plot(binned_vs, 'x')
    ylim((0,50))
    xlim((0,40))

xlabel('r')

show()

from common import *

factors = [0.2, 0.5, 0.8, 1.3, 2, 2.8, 5, 10, 20]

figure(figsize=(14, 10))

def spectra(spec_file, smooth=True):
    hdulist = fits.open(spec_file)
    data = hdulist[1].data
    h, w = data.shape
    cx, cy = w / 2, h / 2
    rs = zeros(w * h)
    vs = zeros(w * h)
    aspect = w / h

    c = 0
    for i in range(0, h):
        for j in range(0, w):
            r = sqrt((cx - j)**2 + (cy - i)**2)
            if r < cx:
                v = data[i][j]
                rs[c] = r
                vs[c] = v
                c +=1

    binned_vs = zeros(int(rs.max() + 1))
    bin_counts = zeros_like(binned_vs)

    for i in range(len(rs)):
        r = rs[i]
        v = vs[i]
        binned_vs[int(floor(r))] += v
        bin_counts[int(floor(r))] += 1

    out =  binned_vs / bin_counts

    if smooth:
        win = bartlett(3)
        win = win / win.sum()
        out = convolve(out, win, 'same')

    return  out

# File names : obs_band_bgfactor_bgn_type.fits
render = loglog
#render = plot

truth = spectra('plots/TRUTH_SPECTRA_m74_PLW.fits')
background = spectra('plots/m74_PLW_0_0_spectra.fits')

count = 1
for factor in factors:
    subplot(3, 3, count)
    count += 1

    spec = spectra('plots/m74_PLW_{}_0_spectra.fits'.format(factor))
    spec_hires = spectra('plots/m74_PLW_{}_0_spectra_hires.fits'.format(factor))
    x = arange(0, len(truth))

    x = x * (1.789257 / 8.5714)
    #normed_hires = spec_hires * (spec[-5:].mean() / spec_hires[-5:].mean())
    #normed_truth = truth * (spec[-5:].mean() / truth[-5:].mean())

    normed_spec = spec# * (truth.max() / spec.max())
    normed_hires = spec_hires * spec.sum() / spec_hires.sum() # / 4 # * (truth.max() / spec_hires.max())
    normed_truth = truth * spec.sum() / truth.sum() # / (truth.mean() / normed_hires.mean()) #* (spec.max() / truth.max())
    x_spec = x[:len(normed_spec)]
    x_hires = x[:len(normed_hires)] #* cdelt_orig / cdelt_hires
    x_truth = x[:len(normed_truth)] #* cdelt_orig / cdelt_truth

    if True:
        render(x_spec, normed_spec, label='Original')
        render(x_spec, background, 'k--', label='Background')
        render(x_hires, normed_hires, label='HiRes')
        render(x_truth, normed_truth, label='Truth')
    else:
        render(x_spec, (normed_spec - background), label='Original')
        render(x_spec, (normed_hires - background), label='HiRes')
        render(x_spec, (normed_truth - background), label='Truth')

    snr = SNR('plots/m74_PLW_{}_0.fits'.format(factor))

    xlabel('Spatial Frequency / arcmin$^{-1}$')
    ylabel('Power')
    title('SNR: {:.1f}'.format(snr))
    legend(loc=3, prop={'size':9})
    xlim(0.2, 18)
    #ylim(0, 1000)

tight_layout()

savefig('doc/report/figures/power-spectra.pdf')
show()

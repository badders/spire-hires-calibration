from common import *
from astropy import wcs
from scipy import interpolate
import aplpy

N = 8
factors = [0.1, 0.5, 1, 2, 3.5, 10]

truth_hdulist = fits.open('plots/TRUTH_m74_PLW.fits')
truth_image = truth_hdulist[1].data
truth_wcs = wcs.WCS(truth_hdulist[1].header)

#normalise
truth_image = truth_image * (0.304 / 3.2)

mask_image = truth_image.copy()
w, h = mask_image.shape
cx, cy = 86, 95
r_max = 30

for i in range(h):
    for j in range(w):
        r = sqrt((cx - j)**2 + (cy - i)**2)
        if r < r_max:
            mask_image[i][j] = 1
        else:
            mask_image[i][j] = 0

truth_image = nan_to_num(truth_image)
truth_image *= mask_image
truth_peak = truth_image.max()

# imshow(truth_image)
# figure()
avg_rms = zeros_like(factors)
avg_hires_rms = zeros_like(avg_rms)

render = plot
snr = array([SNR('plots/m74_PLW_{}_0.fits'.format(f)) for f in factors])

differences = []
hires_differences = []

for n in range(N):
    difference = zeros_like(factors)
    hires_difference = zeros_like(factors)
    diff_maps = []
    hires_diff_maps = []

    for i, f in enumerate(factors):
        hdulist = fits.open('plots/m74_PLW_{}_{}_hires.fits'.format(f, n))
        hires_image = hdulist[1].data
        hdulist2 = fits.open('plots/m74_PLW_{}_{}.fits'.format(f, n))
        image = hdulist2[1].data

        image *= mask_image
        hires_image *= mask_image

        max_loc = truth_image.argmax()
        image *= (nansum(truth_image) / nansum(image))
        hires_image *= (nansum(truth_image) / nansum(hires_image))

        h, w = image.shape
        cx, cy = 86., 95.

        pix_diffs = []
        hires_pix_diffs = []

        diff_maps.append(truth_image - image)
        hires_diff_maps.append(truth_image - hires_image)

        diff = sqrt(nanfilter((truth_image - image)**2).sum() / mask_image.sum())
        hires_diff = sqrt(nanfilter((truth_image - hires_image)**2).sum() / mask_image.sum())
        difference[i] = diff
        hires_difference[i] = hires_diff

    difference = difference / truth_peak
    hires_difference = hires_difference / truth_peak

    differences.append(difference)
    hires_differences.append(hires_difference)

    avg_rms += difference
    avg_hires_rms += hires_difference

    if n == 0:
        fig = figure(figsize=(19,12))
        for i, s in enumerate(snr):
            print(i)
            f1 = aplpy.FITSFigure(fabs(diff_maps[i]), figure=fig , subplot=(3, 4, 2*i +1))
            title('SNR: {:.1f}'.format(s))
            f1.set_theme('publication')
            f1.hide_yaxis_label()
            f1.hide_xaxis_label()
            f1.hide_ytick_labels()
            f1.hide_xtick_labels()
            f1.ticks.hide()
            f1.show_colorscale(cmap='gist_heat', vmin=0., vmax=0.4)
            #xlim(0.4, 0.8)
            #ylim(0.4, 0.8)
            f1.recenter(85, 94, 35)
            f1.add_colorbar()
            f1.colorbar.set_font(size='small')
            f1.set_nan_color('black')

            f2 = aplpy.FITSFigure(fabs(hires_diff_maps[i]), figure=fig , subplot=(3, 4, 2*i + 2))
            title('HiRes')
            f2.set_theme('publication')
            f2.hide_yaxis_label()
            f2.hide_xaxis_label()
            f2.hide_ytick_labels()
            f2.hide_xtick_labels()
            f2.ticks.hide()
            f2.show_colorscale(cmap='gist_heat', vmin=0., vmax=0.4)
            #xlim(0.4, 0.8)
            #ylim(0.4, 0.8)

            f2.recenter(85, 94, 35)
            f2.add_colorbar()
            f2.colorbar.set_font(size='small')
            f2.set_nan_color('black')

        print('SAVING')
        tight_layout()
        savefig('doc/report/figures/diff_maps.pdf')

    import sys
    sys.exit()
    # render(snr, hires_difference, 'kx')
    # render(snr, difference, 'kx')

avg_rms /= N
avg_hires_rms /= N

xnew = linspace(snr[0], snr[-1], num=400)
tck = interpolate.splrep(snr, avg_rms)
hires_tck = interpolate.splrep(snr, avg_hires_rms)

figure()
render(snr, avg_rms, 'b',  label='Truth - Simulation')
render(snr, avg_hires_rms, 'g',  label='Truth - HiRes')

std_devs = zeros_like(factors)
hires_std_devs = zeros_like(factors)

for i in range(len(factors)):
    points = zeros(N)
    hires_points = zeros(N)
    for j in range(len(points)):
        points[j] = differences[j][i]
        hires_points[j] = hires_differences[j][i]

    std_devs[i] = std(points)
    hires_std_devs = std(hires_points)

errorbar(snr, avg_rms, fmt='b.', yerr=std_devs)
errorbar(snr, avg_hires_rms, fmt='g.', yerr=hires_std_devs)


# render(xnew, interpolate.splev(xnew, tck), 'b', label='Simulation')
# render(xnew, interpolate.splev(xnew, hires_tck), 'g', label='HiRes')

legend(loc='best')
xlabel('Peak SNR')
ylabel('RMS Image Pixel Difference (Fraction of Peak)')
grid()
ylim(0.02, 0.2)
xlim(5, snr.max())
savefig('doc/report/figures/differences.pdf')
show()

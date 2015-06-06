from common import *

fig = figure(figsize=(10,6))

factors = [0.2, 0.5, 0.8, 1, 1.3, 1.6,  2, 2.8, 3.5, 5, 7.5, 10, 20][:9]
factors = [0.2, 0.8, 1, 2, 3.5, 20]

def SNR_figure():
    for i in range(len(factors)):
        f = factors[i]
        snr = SNR('plots/m74_PLW_{}_0.fits'.format(f))

        #f1 = aplpy.FITSFigure('plots/m74_PLW_{}_0.fits'.format(f), figure=fig, subplot=(2, 3, i + 1))
        f1 = aplpy.FITSFigure('plots/m74_PLW_{}_0_hires.fits'.format(f), figure=fig, subplot=(2, 3, i + 1))

        title('SNR: {:.1f}'.format(snr))
        f1.set_theme('publication')
        f1.hide_yaxis_label()
        f1.hide_xaxis_label()
        f1.hide_ytick_labels()
        f1.hide_xtick_labels()
        f1.ticks.hide()

        f1.show_colorscale(cmap='gist_heat')
        f1.add_colorbar()
        f1.set_nan_color('black')
        
    tight_layout()
    savefig('doc/report/figures/simulated-observations.pdf')

SNR_figure()
show()

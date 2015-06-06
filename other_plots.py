from common import *
import aplpy

def plot_spectra():
    spec = aplpy.FITSFigure('plots/m74_PLW_20_0_spectra.fits')
    spec.set_theme('publication')
    spec.hide_yaxis_label()
    spec.hide_xaxis_label()
    spec.hide_ytick_labels()
    spec.hide_xtick_labels()
    spec.ticks.hide()
    spec.show_colorscale(cmap='gist_heat', stretch='sqrt')
    spec.add_colorbar()
    tight_layout()
    savefig('doc/report/figures/power-spec-sample.pdf')

def plot_beam():
    spec = aplpy.FITSFigure('fits/beams/0x5000241aL_PLW_pmcorr_1arcsec_norm_beam.fits')
    spec.set_theme('publication')
    spec.hide_yaxis_label()
    spec.hide_xaxis_label()
    spec.hide_ytick_labels()
    spec.hide_xtick_labels()
    spec.ticks.hide()
    spec.show_colorscale(cmap='gist_heat', pmin=0, stretch='power', exponent=1/3)
    spec.add_colorbar()
    tight_layout()
    savefig('doc/report/figures/beam.pdf')

def plot_compare():
    fig = figure(figsize=(5,8))
    spitzer = aplpy.FITSFigure('fits/spitzer/m74.fits', figure=fig, north=True, subplot=(2,1,1))
    spitzer.set_theme('publication')
    spitzer.hide_yaxis_label()
    spitzer.hide_xaxis_label()
    spitzer.hide_ytick_labels()
    spitzer.hide_xtick_labels()
    spitzer.ticks.hide()
    spitzer.show_colorscale(cmap='gist_heat')
    spitzer.add_colorbar()
    spitzer.recenter(24.1735, 15.783, 0.08)

    simulation = aplpy.FITSFigure('fits/herschel/m74.fits', north=True, figure=fig, subplot=(2,1,2))
    simulation.set_theme('publication')
    simulation.hide_yaxis_label()
    simulation.hide_xaxis_label()
    simulation.hide_ytick_labels()
    simulation.hide_xtick_labels()
    simulation.ticks.hide()
    simulation.show_colorscale(cmap='gist_heat')
    simulation.add_colorbar()
    simulation.recenter(24.1735, 15.783, 0.08)

    tight_layout()
    savefig('doc/report/figures/comparison.pdf')


#plot_spectra()
plot_beam()
#plot_compare()
show()

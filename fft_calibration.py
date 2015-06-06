loc = '/Users/tom/herschel/calibration/'

def fix_spectra(spectra):
    dimMap = spectra.dimensions
    specShift = SHIFT(spectra,dimMap[0]/2+1,0)
    specShift = SHIFT(specShift,dimMap[1]/2+1,1)
    imSpec = SimpleImage()
    imSpec.setImage(specShift)
    return imSpec

norm = fitsReader(loc + 'normal.fits')
hires = fitsReader(loc + 'hires.fits')
truth = fitsReader(loc + 'truth.fits')

_, spec_norm =fFT2d(norm)
_, spec_hires =fFT2d(hires)
_, spec_truth =fFT2d(truth)

simpleFitsWriter(fix_spectra(spec_norm), loc + 'normal_spectra.fits')
simpleFitsWriter(fix_spectra(spec_hires), loc + 'hires_spectra.fits')
simpleFitsWriter(fix_spectra(spec_truth), loc + 'truth_spectra.fits')

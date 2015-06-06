import urllib, os

# Settings
poolHsa = 'SPITZER_M74'
poolPath = '/Users/tom/spire/pool-local'
beamPath = '/Users/tom/herschel/fits/beams/'
spitzer24Src = '/Users/tom/herschel/fits/spitzer/m74.fits'
outDir     = "/Users/tom/herschel/plots/"
obsid = 1342189427
bgid = 1342186110
nanval = float('NAN')
bands = ['PLW', 'PMW', 'PSW']
bands = ['PLW']
imgfactors = [0, 0.1, 0.2, 0.5, 0.8, 1, 1.3, 1.6,  2, 2.8, 3.5, 5, 7.5, 10, 20]
#imgfactors = [20]

bg_n = 8
bg_shifts = [[0.0, 0.0], [0.1, 0.1], [0.2, 0.2], [0.0, 0.1], [0.0, 0.2], [0.1, 0.2], [0.2, 0.1], [0.2, 0.0], [0.1, 0.0]]
bg_shift = bg_shifts[bg_n]

# Get the observation from the archive if not already stored in the pool locally
def loadObservation(obsid):
    try:
        print 'Trying to get observation from pool ...'
        obsIn = getObservation(obsid=obsid, poolName=poolHsa, instrument='SPIRE', poolLocation=poolPath)
        obsIn.refs.remove('level0')
        obsIn.refs.remove('level0_5')
    except:
        print 'Failed to get observation from pool, downloading from HSA ...'
        obsIn = getObservation(obsid=obsid, useHsa=True, instrument='SPIRE')
        obsIn.refs.remove('level0')
        obsIn.refs.remove('level0_5')
        spiaSaveObs(obsIn, Pool=poolHsa, nameTag=str(obsid), PoolPath=poolPath)
    return obsIn
    
# Load in the spitzer image
def loadSpitzer():
    spitzer = fitsReader(spitzer24Src)
    # Convert spitzer image from MJ/sr to J/pixel
    sim = spitzer.getImage()
    sim *= 1e6
    sqdeg2sr = 3.04616e-4
    pix_size = abs(spitzer.wcs.cdelt1 * spitzer.wcs.cdelt2) * sqdeg2sr
    sim *= pix_size
    return spitzer

def fix_spectra(spectra):
    dimMap = spectra.dimensions
    specShift = SHIFT(spectra,dimMap[0]/2+1,0)
    specShift = SHIFT(specShift,dimMap[1]/2+1,1)
    imSpec = SimpleImage()
    imSpec.setImage(specShift)
    return imSpec

obsIn = loadObservation(obsid)
obsBG = loadObservation(bgid)
spitzer = loadSpitzer()

spitzerN = {}
spitzerT = {}

for band in bands:
    # Convolve spiter image with SPIRE beams
    # Beams found at https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/
    beamName = "0x5000241aL_%s_pmcorr_1arcsec_norm_beam.fits" % band
    beamFile= beamPath + beamName
    bcenter = {'PSW':(700,699), 'PMW':(700,700), 'PLW':(698,700)} # Positions of peak pixel
    beamfull = fitsReader(beamFile)
    beam = crop(beamfull, int(bcenter[band][0] - 100) , int(bcenter[band][1] - 100), int(bcenter[band][0] + 101), int(bcenter[band][1] + 101))
    
    assert beam.dimensions[0] % 2 == 1, "Beam dimensions are not odd"
    assert beam.dimensions[1] % 2 == 1, "Beam dimensions are not odd"
    assert beam.getIntensity(beam.dimensions[0] / 2, beam.dimensions[1] / 2) == MAX(NAN_FILTER(beam.image)), "Beam is not centred on central pixel"
    
    # Convolve spitzer image with beam
    print 'Convolving spitzer source with beam ...'
    
    print 'Beam cdelt1', beam.wcs.cdelt1
    print 'Spitzer cdelt1', spitzer.wcs.cdelt1
    print 'Beam cdelt2', beam.wcs.cdelt2
    print 'Spitzer cdelt2', spitzer.wcs.cdelt2
    
    # Generate Convultion beams
    nwcs = beam.wcs.copy()
    nwcs.cdelt1 = spitzer.wcs.cdelt1
    nwcs.cdelt2 = spitzer.wcs.cdelt2
    nwcs.naxis1 = spitzer.wcs.naxis1
    nwcs.naxis2 = spitzer.wcs.naxis2
    nwcs.crpix1 = (nwcs.naxis1 + 1) / 2
    nwcs.crpix2 = (nwcs.naxis2 + 1) / 2
    
    # Generate 'truth' beam
    twcs = beam.wcs.copy()
    twcs.cdelt1 = spitzer.wcs.cdelt1 * 2
    twcs.cdelt2 = spitzer.wcs.cdelt2 * 2
    twcs.naxis1 = spitzer.wcs.naxis1
    twcs.naxis2 = spitzer.wcs.naxis2
    twcs.crpix1 = (twcs.naxis1 + 1) / 2
    twcs.crpix2 = (twcs.naxis2 + 1) / 2
    
    # Regrid and fix NaNs
    nbeam = regrid(beam, wcs=nwcs)
    im = nbeam.getImage()
    im[im.where(IS_NAN(im))] = 0
    
    # Normalise back to 1.0 at peak 
    im *= 1/max(im)
    
    # Regrid truth beam
    tbeam = regrid(beam, wcs=twcs)
    tim = tbeam.getImage()
    tim[tim.where(IS_NAN(tim))] = 0
    tim *= 1/max(tim)
    
    convolvedSpitzer = imageConvolution(image=spitzer, kernel=nbeam)
    truthSpitzer = imageConvolution(image=spitzer, kernel=tbeam)
    
    # Remove background level
    im = convolvedSpitzer.getImage()
    im -= MEDIAN(NAN_FILTER(im))
    tim = truthSpitzer.getImage()
    tim -= MEDIAN(NAN_FILTER(tim))
    tim[tim.where(IS_NAN(tim))] = 0
    
    spitzerN[band] = convolvedSpitzer
    spitzerT[band] = truthSpitzer

newScans = {}

# Create new level 1 observation with data from spitzer
for band in bands:
    bgimg =  obsBG.refs["level2"].product.refs['psrc' + band].product
    bg_ra_start = bgimg.meta['ra'].getDouble() + bg_shift[0]
    bg_dec_start = bgimg.meta['dec'].getDouble() + bg_shift[1]
    
    src_ra_start = obsIn.meta['ra'].getDouble()
    src_dec_start = obsIn.meta['dec'].getDouble()
    scanlines = obsIn.level1
    
    newScan = {}
    for imgfactor in imgfactors:
        newScan[imgfactor] = Level1Context()
        newScan[imgfactor].meta = scanlines.meta

    # Process each scan in the existing observation and generate new one with spitzer data and BG Noise
    for prodID in range(scanlines.count):
        print 'Processing scanline %d' % prodID
        line = scanlines.getProduct(prodID)
        # Fetch bolometer names for band we are processing
        bolometers = [name for name in line.getChannelNames() if name.startswith(band)]
        lines = {}
        for imgfactor in imgfactors:
            lines[imgfactor] = scanlines.getProduct(prodID).copy()
           
        for b in bolometers:
            # Generate data for each bolometer
            ra = line.getRa(b)
            dec = line.getDec(b)
            sig = line.getSignal(b)
           
            # Loop over each value and get data from beam convolved spitzer image
            for i in range(len(sig)):
                try:
                    # Load signal from spitzer
                    intensity = spitzerN[band].getIntensityWorldCoordinates(ra[i], dec[i])
                except java.lang.IndexOutOfBoundsException:
                    instensity = nanval
    
                bg_ra = bg_ra_start + (ra[i] - src_ra_start)
                bg_dec = bg_dec_start + (dec[i] - src_dec_start)
                    
                bgintensity =  bgimg.getIntensityWorldCoordinates(bg_ra, bg_dec)
                    
                # Remove any NaNs from source data
                if IS_NAN(intensity):
                     for imgfactor in imgfactors:
                        lines[imgfactor].getSignal(b)[i] = 0
                else:
                    # Add noise factor from our background source
                    for imgfactor in imgfactors:
                        lines[imgfactor].getSignal(b)[i] = intensity * imgfactor + bgintensity
                   
        # Add scanline to new product
        for imgfactor in imgfactors:
            newScan[imgfactor].addProduct(lines[imgfactor])
    newScans[band] = newScan

# Create maps with new data
imagesize = [20, 20] # y-, x- dimensions in arcminutes
# Create hires map
print 'Generating HiRes maps ...'

maps_hires = {}
spectra_hires = {}
hist_hires = {}
beam_hires = {}

for band in bands:
    wcs = obsIn.level2.getProduct("psrc"+band).wcs.copy()
    # Prepare Wcs with half the pixel size of standard map
    wcs.cdelt1 /= 2.0
    wcs.cdelt2 /= 2.0
    wcs.naxis1 = int(imagesize[1]/60./abs(wcs.cdelt1) + 0.5)
    wcs.naxis2 = int(imagesize[0]/60./abs(wcs.cdelt2) + 0.5)
    wcs.crpix1 = (wcs.naxis1 + 1) / 2.
    wcs.crpix2 = (wcs.naxis2 + 1) / 2.

    newmaps = {}
    newspectra = {}
    newbeams = {}
    newhists = {}
    for imgfactor in imgfactors:
        newmaps[imgfactor], newbeams[imgfactor] = hiresMapper(newScans[band][imgfactor], array=band, beam=beam, wcs=wcs)
        _,  nspectra = fFT2d(newmaps[imgfactor])
        newspectra[imgfactor] = fix_spectra(nspectra)
        newhists[imgfactor] = imageHistogram(newmaps[imgfactor])
    
    maps_hires[band] = newmaps
    beam_hires[band] = newbeams
    spectra_hires[band] = newspectra
    
    hist_hires[band] = newhists

maps = {}
spectra = {}

print 'Generating maps ...'

for band in bands:
    wcs = obsIn.level2.getProduct("psrc"+band).wcs.copy()
    imagecenter = [wcs.crval1, wcs.crval2]
    
    wcs.naxis1 = int(imagesize[1]/60./abs(wcs.cdelt1) + 0.5)
    wcs.naxis2 = int(imagesize[0]/60./abs(wcs.cdelt2) + 0.5)
    wcs.crpix1 = (wcs.naxis1 + 1) / 2.
    wcs.crpix2 = (wcs.naxis2 + 1) / 2.
    print(band)
   
    newmaps = {}
    newspectra = {}
    for imgfactor in imgfactors:
        newmaps[imgfactor] = naiveScanMapper(newScans[band][imgfactor], array=band, wcs=wcs)
        newmaps[imgfactor] = regrid(source=newmaps[imgfactor], target=maps_hires[band][imgfactor])
        _,  nspectra = fFT2d(newmaps[imgfactor])
        newspectra[imgfactor] = fix_spectra(nspectra)
    
    maps[band] = newmaps
    spectra[band] = newspectra
    
print 'Saving output images ...'
for band in bands:
    spitzer_convolved = spitzerN[band]
    spitzer_truth = spitzerT[band]
    spitzer_convolved = regrid(source=spitzer_convolved, target=maps_hires[band][imgfactors[0]])
    simpleFitsWriter(spitzer_convolved, outDir + 'SPITZER_m74_' + band + '.fits')
    spitzer_truth = regrid(source=spitzer_truth, target=maps_hires[band][imgfactors[0]])
    simpleFitsWriter(spitzer_truth, outDir + 'TRUTH_m74_' + band + '.fits')
    tim_wn = spitzer_truth.getImage().copy()
    tim = spitzer_truth.getImage()
    tim[tim.where(IS_NAN(tim))] = 0
    _,  nspectra = fFT2d(spitzer_truth)
    simpleFitsWriter(fix_spectra(nspectra), outDir + 'TRUTH_SPECTRA_m74_' + band + '.fits')
        
    for imgfactor in imgfactors:
        map = maps[band][imgfactor]
        im = map.getImage()
        im[tim_wn.where(IS_NAN(tim_wn))] = nanval
        hires = maps_hires[band][imgfactor]
        him = hires.getImage()
        him[tim_wn.where(IS_NAN(tim_wn))] = nanval
        
        beam = beam_hires[band][imgfactor]
        hires_spectra = spectra_hires[band][imgfactor]
        spec = spectra[band][imgfactor]
        simpleFitsWriter(map, outDir + 'm74_'+ band + '_' + str(imgfactor) + '_' + str(bg_n) + '.fits')
        simpleFitsWriter(beam, outDir + 'm74_'+ band + '_' + str(imgfactor) + '_' + str(bg_n) + '_beam.fits')
        simpleFitsWriter(hires, outDir + 'm74_'+ band + '_' + str(imgfactor) + '_' + str(bg_n) + '_hires.fits')
        simpleFitsWriter(spec, outDir + 'm74_'+ band + '_' + str(imgfactor) + '_' + str(bg_n) + '_spectra.fits')
        simpleFitsWriter(hires_spectra, outDir + 'm74_'+ band + '_' + str(imgfactor) + '_' + str(bg_n) + '_spectra_hires.fits')
import urllib, os

poolHsa = 'SPITZER_M74'
poolPath = '/Users/tombadran/hcss/pool'

spitzer24Src = '/Users/tombadran/herschel/fits/spitzer/m74.fits'

cal = spireCal()

obsid = 1342189427
band = 'PLW'

# Get the observation from the archive if not already stored in the pool locally
def loadObservation(obsid):
    try:
        print 'Trying to get observation from pool ...'
        obsIn = getObservation(obsid=obsid, poolName=poolHsa, instrument='SPIRE', poolLocation=poolPath)
    except:
        print 'Failed to get observation from pool, downloading from HSA ...'
        obsIn = getObservation(obsid=obsid, useHsa=True, instrument='SPIRE')
        obsIn.refs.remove('level0')
        obsIn.refs.remove('level0_5')
        obsIn.refs.remove('level2')
        spiaSaveObs(obsIn, Pool=poolHsa, nameTag=str(obsid), PoolPath=poolPath)
    return obsIn
    
# Load in the spitzer image
def loadSpitzer():
    srcFits = fitsReader(spitzer24Src)
    smoothedSrc = gaussianSmoothing(image=srcFits, sigma=5)
    croppedImage = crop(image=smoothedSrc, row1=197, column1=800, row2=1448, column2=2189)
    return croppedImage

obsIn = loadObservation(obsid)
spitzer = loadSpitzer()

# Convolve spiter image with SPIRE beam
workDir = Configuration.getProperty('var.hcss.workdir')
beamName = "0x5000241aL_%s_pmcorr_1arcsec_norm_beam.fits" % band
urllib.urlretrieve ("https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/"+beamName, os.path.join(workDir,beamName))

bcenter = {'PSW':(700,699), 'PMW':(700,700), 'PLW':(698,700)} # Positions of peak pixel
beamfull = fitsReader(file = os.path.join(workDir,beamName))
beam = crop(beamfull, int(bcenter[band][0] - 100) , int(bcenter[band][1] - 100), int(bcenter[band][0] + 101), int(bcenter[band][1] + 101))

assert beam.dimensions[0] % 2 == 1, "Beam dimensions are not odd"
assert beam.dimensions[1] % 2 == 1, "Beam dimensions are not odd"
assert beam.getIntensity(beam.dimensions[0] / 2, beam.dimensions[1] / 2) == MAX(NAN_FILTER(beam.image)), "Beam is not centred on central pixel"

beamFT, _ = fFT2d(beam)
spitzerFT, _ = fFT2d(spitzer)
spitzerFt = beamFT * spitzerFT

spitzer, _ = fFT2d(spitzerFT)

# Overwrite herschel observation with data from spitzer
scanlines = obsIn.level1
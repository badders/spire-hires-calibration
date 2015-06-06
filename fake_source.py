###
##Add Gaussian source to timelines and maps, with width based on beam profile
###

import SpireHandbookBundle as hb
import sources_dev as srcMod
import os,sys,math

poolHsa='CheckApCorr_HSA'
poolSrc='CheckApCorr_Src'
poolPath='/Users/tombadran/hcss/pool'

bands = ['PSW','PMW','PLW']

cal=spireCal()

skipDone=False

#set observation properties
obsid=1342186110 #Lockman North (1deg wide relatively shallow deep field map)
raList=[-10.,0.,10.] #arcmin
decList=[-10.,0.,10.] #arcmin

#obsid=1342212430 #Alpha Ari (1)
#raList=[-1.,0.,1.] #arcmin
#decList=[-1.,0.,1.] #arcmin

locList=[]
for ra in raList:
    for dec in decList:
        locList.append([ra,dec])
alphaArr=Float1d.range(10)-4
   
alphaPip=-1

try:
    obsIn=getObservation(obsid=obsid,poolName=poolHsa,instrument='SPIRE',poolLocation=poolPath)
except:
    obsIn=getObservation(obsid=obsid,useHsa=True,instrument='SPIRE')
    obsIn.refs.remove('level0')
    obsIn.refs.remove('level0_5')
    obsIn.refs.remove('axiliary')
    spiaSaveObs(obsIn,Pool=poolHsa,nameTag=str(obsid),PoolPath=poolPath)
    
#set source peak flux density location(s) [wrt centre of map], and spectral index(es)
srcAmplitudes=[50.,100.]
for srcAmplitude in srcAmplitudes:
    
    areaPip={}
    k4P={}
    k4E={}
    for band in bands:
        areaPip[band]=cal.getPhot().getProduct('RadialCorrBeam').meta['beamPipeline%sArc'%band.capitalize()].value
        k4P[band]=cal.getPhot().getProduct('FluxConvList')[0].meta['k4P_%s'%band].value
        k4E[band]=cal.getPhot().getProduct('FluxConvList')[0].meta['k4E_%s'%band].value
    for alpha in alphaArr:
        kPsrc={}
        srcAmp={}
        for band in bands:
            kPsrc[band]=cal.getPhot().getProduct('ColorCorrKList').getProduct('point').getAlphaCorrection(alpha,band)
            srcAmp[band]=srcAmplitude/kPsrc[band]
        
        apRad={'PSW':22.,'PMW':30.,'PLW':45,}
        effBeam = hb.calcSpireEffBeam(alpha)
        effBeamPip = hb.calcSpireEffBeam(alphaPip)
        #apCorr no background
        apCorr = hb.calcApCorr(alpha)[0]
        #apCorr inc background subtraction
        apCorrBg = hb.calcApCorr(alpha)[1]
        
        effBeamSrc={}
        effBeamSrcPip={}
        effBeamMap={}
        effBeamSrcInt={}
        effBeamMapArea={}
        effBeamSrcArea={}
        effBeamSrcPipArea={}
        effBeamSrcFwhm={}
        effBeamSrcSigma={}
        effBeamSrcPipFwhm={}
        effBeamApPhot={}
        effBeamApPhotFlux={}
        effBeamApPhotFluxCorr={}
        extent={}
        
        for band in bands:
            effBeamSrc[band]=srcMod.SourceProfile(Double1d.range(len(effBeam[band])),effBeam[band])
            #normalise by area
            extent[band]=(effBeamSrc[band].maxRad-10.)/3600. # extent of source profile
            
            #calculate area and fwhm
            #effBeamMapArea[band]=sum(effBeamMap[band].image)
            #print '%.1f %s map area: %f'%(alpha,band,effBeamMapArea[band])
            effBeamSrcArea[band]=effBeamSrc[band].calcArea()
            print '%.1f %s source area: %f'%(alpha,band,effBeamSrcArea[band])
            effBeamSrcFwhm[band]=effBeamSrc[band].calcFwhm()
            print '%.1f %s source FWHM: %f'%(alpha,band,effBeamSrcFwhm[band])
            effBeamSrcSigma[band]=effBeamSrcFwhm[band]/SQRT(8.*LOG(2.))
                
        #get zero point offsets
        zPointOffset={}
        for band in bands:
            zPointOffset[band]=obsIn.getLevel2().getProduct('extd%s'%band).meta['zPointOffset']
        #
        for loc in locList:
            #construct tag of output
            dRa=loc[0] #arcmin
            dDec=loc[1] #arcmin
            tagSrc='%d_Gauss_alpha%.1f_%.1fJy_%d_%d'%(obsid,alpha,srcAmplitude,dRa,dDec)
            
            #check whether obs already exists
            if skipDone:
                try:
                    obsDone=getObservation(obsid=obsid,tag=tagSrc,poolName=poolSrc,instrument='SPIRE',poolLocation=poolPath)
                    if skipDone:
                        print '%d %s already constructed. Skipping'%(obsid,tagSrc)
                        continue
                    else:
                        print 'Reconstructing %d %s.'%(obsid,tagSrc)
                except:
                    pass
            print 'Constructing %d %s.'%(obsid,tagSrc)
            #
            #make blank observation
            print 'making blank observation context'
            obsOut=ObservationContext()
            obsOut.meta = obsIn.meta
            #
            #copy scans context
            scansIn=obsIn.level1.copy()
            #
            scansSrc=Level1Context()
            scansSrc.meta=scansIn.meta
            #
            nScan=0
            nScanTot=len(scansIn.refs)
            centre=[obsIn.meta['ra'].value +dRa/60./COS(obsIn.meta['dec'].value*math.pi/180.),obsIn.meta['dec'].value + dDec/60.]
            obsOut.meta['srcType']=StringParameter('Gaussian',description='Type of artificial source')
            obsOut.meta['srcRa']=DoubleParameter(centre[0],unit=herschel.share.unit.Angle.DEGREES,\
              description='RA of artificial source')
            obsOut.meta['srcDec']=DoubleParameter(centre[1],unit=herschel.share.unit.Angle.DEGREES,\
              description='Dec of artificial source')
            obsOut.meta['srcPeak']=DoubleParameter(srcAmplitude,unit=herschel.share.unit.FluxDensity.JANSKYS,\
              description='Peak flux density of artificial source')
            obsOut.meta['srcAlpha']=DoubleParameter(alpha,unit=herschel.share.unit.FluxDensity.JANSKYS,\
              description='Spectral index of artificial source')
            obsOut.meta['srcFwhm']=DoubleParameter(effBeamSrcFwhm[band],unit=herschel.share.unit.Angle.SECONDS_ARC,\
              description='FWHM of source')
            obsOut.meta['srcArea']=DoubleParameter(effBeamSrcArea[band],unit=herschel.share.unit.SolidAngle.SQUARE_SECONDS_ARC,\
              description='Solid angle of source')
            print 'adding %.1f Jy nu^%.1f source to %d scans of %d at (%.2f,%.2f) [%.1f,%.1f arcmin from centre]'%(srcAmplitude,alpha,nScanTot,obsid,centre[0],centre[1],dRa,dDec)
            #
            ###add source to obsid timeline
            #add source to timelines
            #first=True
            for pdtref in scansIn.refs:
                nScan+=1
                #print 'adding source to scan %d of %d'%(nScan,nScanTot)
                #
                pdt=pdtref.product
                bolometers=pdt.getChannelNames()
                for bolo in bolometers:
                    ## if first time then calculate where the centre is best
                    #if first:
                    #    # create centre
                    #    centre = [0.0,0.0]
                    #    centre[0] = (max(pdt.getRa("PSWE8")) - min(pdt.getRa("PSWE8"))) / 2.0 + min(pdt.getRa("PSWE8"))
                    #    centre[1] = (max(pdt.getDec("PSWE8")) - min(pdt.getDec("PSWE8"))) / 2.0 + min(pdt.getDec("PSWE8"))
                    #    first = False
                    #
                    # skip anything that is not a bolometer
                    if bolo[3:4] == "T" or bolo[3:5] == "DP" or bolo[3:4] == "R":
                        continue
                    #
                    if bolo.find('PSW')>=0:
                        band='PSW'
                    elif bolo.find('PMW')>=0:
                        band='PMW'
                    elif bolo.find('PLW')>=0:
                        band='PLW'
                    #
                    # locate all samples that are in the range required
                    boloRA = pdt.getRa(bolo)
                    boloDEC = pdt.getDec(bolo)
                    boloSignal = pdt.getSignal(bolo)
                    selection = boloRA.where(((boloRA - centre[0]) * math.cos(centre[1] * math.pi / 180.0))**2.0 + (boloDEC - centre[1])**2.0 < extent[band]**2.0)
                    objRA = boloRA[selection]
                    objDEC = boloDEC[selection]
                    objSignal = boloSignal[selection]
                    #
                    # add the source to the timeline
                    for i in range(0,len(objSignal)):
                        #calculate rad in degrees
                        radSqu = ((objRA[i] -centre[0])*math.cos(objDEC[i] * math.pi/180.0))**2.0 + (objDEC[i]-centre[1])**2.0
                        #add Gaussian to timeline
                        objSignal[i] = objSignal[i] + srcAmp[band] * math.exp(-radSqu / (2.0*(effBeamSrcSigma[band]/3600.)**2.0))
                    # add modified values back to full array
                    boloSignal[selection] = objSignal
                    #
                    # set modified array back to pdt
                    pdt.setSignal(bolo, boloSignal)
                    #
                scansSrc.addProduct(pdt)
            #
            #add scans to obsOut
            obsOut.setLevel1(scansSrc)
            #
            mapsOut=MapContext()
            mapsOut.meta = obsIn.level2.meta
            #make maps from scans
            for band in ['PSW','PMW','PLW']:
                print 'making psrc %s map'%band
                map=naiveScanMapper(input=scansSrc,array=band,minVel=5.0)
                #add metadata
                map.meta['srcType']=StringParameter('Gaussian',description='Type of artificial source')
                map.meta['srcRa']=DoubleParameter(centre[0],unit=herschel.share.unit.Angle.DEGREES,\
                  description='RA of artificial source')
                map.meta['srcDec']=DoubleParameter(centre[1],unit=herschel.share.unit.Angle.DEGREES,\
                  description='Dec of artificial source')
                map.meta['srcPeak']=DoubleParameter(srcAmplitude,unit=herschel.share.unit.FluxDensity.JANSKYS,\
                  description='Peak flux density of artificial source')
                map.meta['srcAlpha']=DoubleParameter(alpha,unit=herschel.share.unit.FluxDensity.JANSKYS,\
                  description='Spectral index of artificial source')
                map.meta['srcFwhm']=DoubleParameter(effBeamSrcFwhm[band],unit=herschel.share.unit.Angle.SECONDS_ARC,\
                  description='FWHM of source')
                obsOut.meta['srcArea']=DoubleParameter(effBeamSrcArea[band],unit=herschel.share.unit.SolidAngle.SQUARE_SECONDS_ARC,\
                  description='Solid angle of source')
                #convert to extended
                print 'making extd %s map'%band
                mapExtd=convertImageUnit(image=map,newUnit='MJy/sr',beamArea=areaPip[band])
                mapExtd=imageDivide(image1=mapExtd,scalar=k4P[band])
                mapExtd=imageMultiply(image1=mapExtd,scalar=k4E[band])
                mapExtd=imageAdd(image1=mapExtd,scalar=zPointOffset[band].value)
                mapExtd.meta['zPointOffset']=zPointOffset[band]
                #add to MapContext
                mapsOut.setProduct('psrc%ssrc'%band,map)
                mapsOut.setProduct('extd%ssrc'%band,mapExtd)
            #
            obsOut.setLevel2(mapsOut)
            spiaSaveObs(obsOut,nameTag=tagSrc,Pool=poolSrc,PoolPath=poolPath)


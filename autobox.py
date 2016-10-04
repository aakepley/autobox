from taskinit import * # This gets me the toolkit tasks.

def runTclean(paramList, sidelobeThreshold,noiseThreshold, peakThreshold,
              smoothFactor,cutThreshold, save=False, minBeamFrac=-1):
    '''
    run clean
    '''

    import os
    import shutil
    import numpy
    import copy
    import time;
    
    from refimagerhelper import PySynthesisImager
    from refimagerhelper import ImagerParameters, PerformanceMeasure

    paramList.printParameters()

    imager = PySynthesisImager(params=paramList)

    # init major cycle elements
    imager.initializeImagers()    
    imager.initializeNormalizers()
    imager.setWeighting()

    ## Init minor cycle elements
    imager.initializeDeconvolvers()
    imager.initializeIterationControl()            
    imager.makePSF()
    imager.makePB()
    
    ## get image names
    (image,psfImage,residImage,maskImage) = getImageNames(paramList)
    
    ## Make dirty image
    imager.runMajorCycle()

    # determine sidelobe level
    sidelobeLevel = fitPSF(psfImage)
    
    # calculate size of smoothing beam
    smoothKernel = calcSmooth(psfImage,smoothFactor)

    # determine clean threshold
    iterPars = paramList.getIterPars()
    cleanThreshold = float(iterPars['threshold'][0:-2])
    
    ## Do deconvolution and iterations
    while ( not imager.hasConverged() ):

        # determine peak of residual
        (residPeak, residRMS) = findResidualStats(residImage)
        
        # calculate thresholds
        thresholdNameList = ['Sidelobe','Peak','Noise']
        thresholdList = [sidelobeThreshold*sidelobeLevel*residPeak,
                         peakThreshold*residPeak,
                         max(noiseThreshold*residRMS,cleanThreshold)]

        ## compare various thresholds -- Do I want an absolute value here?
        maskThreshold = max(thresholdList)
        maskThresholdIdx = thresholdList.index(maskThreshold)
        
        casalog.post("Using " + thresholdNameList[maskThresholdIdx] + " threshold: " +  str(maskThreshold), origin='autobox')

        # Print out values for all thresholds    
        for (name,value) in zip(thresholdNameList,thresholdList):
            casalog.post( name + " threshold is " + str(value), origin='autobox')
     
        # create a new mask. If I need to do pruning, it should be done here.
        maskRoot = 'tmp'
        calcMask(residImage,psfImage,
                 maskThreshold,smoothKernel,cutThreshold,maskRoot=maskRoot,
                 minBeamFrac=minBeamFrac)
        
        # add masks together
        if imager.ncycle > 0:
            addMasks(maskImage+str(imager.ncycle-1),maskRoot+'_final_mask',maskRoot+'_final_mask_sum')
            casalog.post( 'adding mask '+ maskImage+str(imager.ncycle-1) + ' and ' + maskRoot+'_final_mask',origin='autobox')
            finalMaskImage = maskRoot+'_final_mask_sum'
        else: 
            finalMaskImage = maskRoot + '_final_mask'

        if os.path.exists(maskImage):
            shutil.rmtree(maskImage)

        casalog.post( 'copying ' + finalMaskImage + ' to '+ maskImage, origin='autobox')
        shutil.copytree(finalMaskImage,maskImage)
        casalog.post( 'copying ' +maskImage + ' to '+ maskImage+str(imager.ncycle),origin='autobox')
        shutil.copytree(maskImage,maskImage+str(imager.ncycle))
    
        # save intermediate masks and residuals for diagnostics        
        if save:
            shutil.copytree(maskRoot+'_mask',maskRoot+'_mask'+str(imager.ncycle))
            shutil.copytree(maskRoot+'_smooth_mask',maskRoot+'_smooth_mask'+str(imager.ncycle))
            shutil.copytree(maskRoot+'_final_mask',maskRoot+'_final_mask'+str(imager.ncycle))
            if imager.ncycle > 0:
                shutil.copytree(maskRoot+'_final_mask_sum',maskRoot+'_final_mask_sum'+str(imager.ncycle))
            shutil.copytree(residImage,residImage+str(imager.ncycle))



        imager.runMinorCycle() 
        imager.runMajorCycle()
        
    imager.restoreImages()
    imager.pbcorImages()
                    
    ## Close tools.
    imager.deleteTools() 


def getImageNames(paramList):
    '''
    Get image names
    '''    

    ## get image names
    if paramList.allimpars['0']['nterms'] < 2:
        image = paramList.alldecpars['0']['imagename']+'.image'
        psfImage = paramList.alldecpars['0']['imagename']+'.psf'
        residImage = paramList.alldecpars['0']['imagename']+'.residual'
    else:
        image = paramList.alldecpars['0']['imagename']+'.image.tt0'
        psfImage = paramList.alldecpars['0']['imagename']+'.psf.tt0'
        residImage = paramList.alldecpars['0']['imagename']+'.residual.tt0'

    #get name of current mask
    if paramList.alldecpars['0']['mask']:
        maskImage = paramList.alldecpars['0']['mask']
    else:
        maskImage = paramList.alldecpars['0']['imagename']+'.mask'

    return (image,psfImage,residImage,maskImage)

def fitPSF(psfname,boxpixels=20):
    '''
    fit psf and then subtract psf. return the minimum psf sidelobe
    value. Code is based on VLASS example code.
    '''

    # create residual image
    psfresid = psfname + '.resid'

    # fit and subtract gaussian from PSF
    psf = ia.newimagefromimage(infile=psfname,outfile=psfresid,overwrite=True)
    psfimstat = psf.statistics()
    blcx=psfimstat['maxpos'][0]-boxpixels
    trcx=psfimstat['maxpos'][0]+boxpixels
    blcy=psfimstat['maxpos'][1]-boxpixels
    trcy=psfimstat['maxpos'][1]+boxpixels
    blctrc=str(blcx)+','+str(blcy)+','+str(trcx)+','+str(trcy)
    clrec= psf.fitcomponents(box=blctrc)
    psf.modify(clrec['results'],subtract=True)
    
    # find the maximum sidelobe.
    psfresidimstat=psf.statistics()
    psfmin=max(abs(psfresidimstat['min'][0]),psfresidimstat['max'][0])
    psf.done()

    return psfmin

def calcSmooth(psfImage, smoothFactor):
    '''
    calculate the size of the kernel to smooth to.
    '''
    
    ia.open(psfImage)
    beam = ia.commonbeam()
    ia.done()

    beam['major']['value'] = beam['major']['value'] * smoothFactor
    beam['minor']['value'] = beam['minor']['value'] * smoothFactor

    return beam
        
def findResidualStats(residImage):
    '''
    calculate the peak of the residual
    '''

    ia.open(residImage)
    residStats = ia.statistics()
    ia.done()
    residPeak = residStats['max'][0] ## do I want to make this the absolute value of the max/min??
    residRMS = residStats['rms'][0]

    return (residPeak, residRMS)




def calcMask(residImage,psfImage, maskThreshold,smoothKernel,cutThreshold,
             minBeamFrac=-1,
             maskRoot='tmp'):
    '''
    calculate the mask based on the residual image and the mask threshold
    '''

    import math

    ia.open(residImage)
    tmpMaskName = maskRoot+'_mask'
    tmpMask = ia.imagecalc(tmpMaskName,'iif('+residImage+'>'+str(maskThreshold)+',1.0,0.0)',overwrite=True)
    ia.done()
    tmpMask.done()

    if minBeamFrac > 0:
        casalog.post("pruning regions smaller than " + str(minBeamFrac) + "times the beam size",origin='autobox')
        pruneRegions(psfImage,tmpMaskName,minBeamFrac)

    major = smoothKernel['major']['value']
    minor = smoothKernel['minor']['value']
    pa = smoothKernel['pa']['value']

    majorStr = str(major)+smoothKernel['major']['unit']
    minorStr = str(minor)+smoothKernel['minor']['unit']
    paStr = str(pa)+smoothKernel['pa']['unit']

    casalog.post("smoothing by " + majorStr + " by " + minorStr,origin='autobox')

    tmpSmoothMaskName = maskRoot+'_smooth_mask'
    ia.open(tmpMaskName)
    tmpSmoothMask = ia.convolve2d(outfile=tmpSmoothMaskName,axes=[0,1],type='gauss',
                                       major=majorStr,minor=minorStr,pa=paStr,
                                       overwrite=True)
    
    tmpSmoothMaskStats = tmpSmoothMask.statistics()
    tmpSmoothMaskPeak = tmpSmoothMaskStats['max'][0]
    
    tmpFinalMaskName = maskRoot+'_final_mask'
    tmpFinalMask = tmpSmoothMask.imagecalc(tmpFinalMaskName,
                                           'iif('+tmpSmoothMaskName+'>'+str(cutThreshold*tmpSmoothMaskPeak)+',1.0,0.0)',
                                           overwrite=True)
    
    tmpSmoothMask.done()
    tmpFinalMask.done()
                                        
    ia.done()
  
def addMasks(mask1,mask2,mask3):
    '''
    Add masks together to create final mask
    '''

    ia.open(mask1)
    ia.close()
    ia.open(mask2)
    ia.close()

    tmp = ia.imagecalc(mask3, 'iif('+mask1+'+'+mask2+'>=1.0,1.0,0.0)', overwrite=True)

    tmp.done()
    ia.done()

def pruneRegions(psfImage,maskImage,minBeamFrac):

    ''' 
    code to prune regions that are too small. Based on code in M100 casaguide
    developed by Crystal Brogan and Jen Donovan Meyer.
    '''

    import scipy.ndimage
    import analysisUtils as au

    # get beam area in pixels using analysis utils
    beamArea = au.pixelsPerBeam(psfImage) 
    npix = beamArea * minBeamFrac
    
    # open image and get mask
    ia.open(maskImage)
    mask = ia.getchunk()
    
    # divide the mask up into labeled regions
    labeled, j = scipy.ndimage.label(mask)
    myhistogram = scipy.ndimage.measurements.histogram(labeled,0,j+1,j+1)
    object_slices = scipy.ndimage.find_objects(labeled)

    # go through each region and set the mask to 0 if there aren't enough pixels
    nremoved = 0
    for i in range(j):
        if myhistogram[i+1] < npix:
            mask[object_slices[i]] = 0
            nremoved =+ 1

    casalog.post("removed " + str(nremoved) + " regions",origin='autobox')

    # put the mask back
    ia.putchunk(mask)
    ia.done()
    



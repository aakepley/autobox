from taskinit import * # This gets me the toolkit tasks.

def runTclean(paramList, 
              sidelobeThreshold,lowNoiseThreshold, noiseThreshold, minBeamFrac=-1,
              peakThreshold, 
              smoothFactor, cutThreshold,  save=False):
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

    import analyzemsimage

    # set some defaults 
    if not lowNoiseThreshold:
        lowNoiseThreshold = noiseThreshold

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
    (image,psfImage,residImage,maskImage,pbImage) = getImageNames(paramList)
    
    ## Make dirty image
    imager.runMajorCycle()

    # determine sidelobe level
    sidelobeLevel = fitPSF(psfImage)
    
    # calculate size of smoothing beam
    smoothKernel = calcSmooth(psfImage,smoothFactor)

    # determine clean threshold
    iterPars = paramList.getIterPars()
    cleanThreshold = float(iterPars['threshold'][0:-2])
    normPars = paramList.getNormPars()
    pbLevel = normPars['0']['pblimit']

    # initially do not grow masks.
    growmask = False

    ## Do deconvolution and iterations
    while ( not imager.hasConverged() ):

        # determine peak and RMS of residual image
        (residPeak, residRMS) = findResidualStats(residImage)

        # determine noise in image in annulus.
        imageStats = analyzemsimage.imstatAnnulus(residImage,pbimage=pbImage,innerLevel=0.3,outerLevel=0.2,verbose=True)

        if len(imageStats['medabsdevmed']) > 0:
            imageRMS = imageStats['medabsdevmed'][0] * 1.4826
        else:
            imageRMS = 0.0

        # calculate thresholds -- I've switched to using imageRMS
        # here, but could go back to residual rms.
        thresholdNameList = ['Sidelobe','Peak','Noise']
        thresholdList = [sidelobeThreshold*sidelobeLevel*residPeak,
                         peakThreshold*residPeak,
                         noiseThreshold*imageRMS] 

        ## compare various thresholds -- Do I want an absolute value here?
        maskThreshold = max(thresholdList)
        maskThresholdIdx = thresholdList.index(maskThreshold)
        
        ## set the growth flag if residPeak is below the
        ## noiseThreshold*imageRMS. Note not setting per minor cycle, setting to OFF, then ON.
        if ((residPeak < noiseThreshold*imageRMS) and (thresholdNameList[maskThresholdIdx] == 'Noise')):
            growmask = True
        
        lowThreshold = lowNoiseThreshold * imageRMS

        casalog.post("Using " + thresholdNameList[maskThresholdIdx] + " threshold: " +  str(maskThreshold), origin='autobox')

        # Print out values for all thresholds    
        for (name,value) in zip(thresholdNameList,thresholdList):
            casalog.post( name + " threshold is " + str(value), origin='autobox')

        # let user know if growing mask
        casalog.post("Grow Mask: " + str(growmask),origin='autobox')

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

        # Grow masks here? I think I want the masks from previous cycles to be present.
        if growmask:
            growMask(residImage,finalMaskImage,lowThreshold)
            finalMaskImage = 'tmp_growmask'

        # moving masks around to make the available for the final clean cycle
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
        pbImage = paramList.alldecpars['0']['imagename']+'.pb'
    else:
        image = paramList.alldecpars['0']['imagename']+'.image.tt0'
        psfImage = paramList.alldecpars['0']['imagename']+'.psf.tt0'
        residImage = paramList.alldecpars['0']['imagename']+'.residual.tt0'
        pbImage = paramList.alldecpars['0']['imagename']+'.pb.tt0'

    #get name of current mask
    if paramList.alldecpars['0']['mask']:
        maskImage = paramList.alldecpars['0']['mask']
    else:
        maskImage = paramList.alldecpars['0']['imagename']+'.mask'

    return (image,psfImage,residImage,maskImage,pbImage)

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
    residStats = ia.statistics(robust=True)
    ia.done()
    residPeak = residStats['max'][0] ## do I want to make this the absolute value of the max/min??
    #residRMS = residStats['rms'][0]
    
    residRMS = residStats['medabsdevmed'][0] * 1.4826

    casalog.post('Using medabsdev of ' + str(residStats['medabsdevmed'][0]) +  ', which corresponds to sigma of ' + str(residRMS), origin='autobox')


    return (residPeak, residRMS)



def calcMask(residImage,psfImage, maskThreshold,smoothKernel,cutThreshold,
             minBeamFrac=-1,
             maskRoot='tmp'):
    '''

    calculate the mask based on the residual image and the mask threshold. 
    
    This should probably be refactored to make it easier to use. Code
    is pretty spaghetti line now.

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

    # smoothing the mask
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
            nremoved = nremoved + 1

    casalog.post("removed " + str(nremoved) + " regions",origin='autobox')

    # put the mask back
    ia.putchunk(mask)
    ia.done()
    
def growMask(residImage,maskImage,lowThreshold,iterations=10):

    ''' 
    grow mask through binary dilation
    '''

    import matplotlib.pyplot as plt
    from scipy.ndimage import binary_dilation, generate_binary_structure
    import numpy as np

    # getting the residual image and the existing mask
    ia.open(residImage)
    residArray = ia.getchunk(dropdeg=True)
    ia.close()
    ia.open(maskImage)
    highMask = ia.getchunk(dropdeg=True)
    ia.close()
    ia.done()
    
    # setting the lower contour mask
    lowMask = residArray > lowThreshold

    # dilating the binary mask into the low contour
    struct = generate_binary_structure(2,1).astype(highMask.dtype)
    growMask = binary_dilation(highMask,structure=struct,iterations=iterations,mask=lowMask).astype(highMask.dtype)

    # saving the new mask
    tmp = ia.newimagefromimage(infile=maskImage,outfile='tmp_growmask',overwrite=True)
    tmp.putchunk(growMask)
    tmp.done()
    ia.done()

        

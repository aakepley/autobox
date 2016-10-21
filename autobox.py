from taskinit import * # This gets me CASA toolkit tasks.

def runTclean(paramList, 
              sidelobeThreshold,lowNoiseThreshold, noiseThreshold,  
              smoothFactor, cutThreshold,
              minBeamFrac=-1):
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
    (image,psfImage,residImage,maskImage,pbImage) = getImageNames(paramList)
    
    ## Make dirty image
    imager.runMajorCycle()

    # determine sidelobe level
    sidelobeLevel = fitPSF(psfImage)

    casalog.post("Fractional Sidelobe Level: " + str(sidelobeLevel),origin='autobox')
    
    # calculate size of smoothing beam
    smoothKernel = calcSmooth(psfImage,smoothFactor)

    # Do deconvolution 
    while ( not imager.hasConverged()):
        
        # determine peak and RMS of residual image
        (residPeak, residRMS) = findResidualStats(residImage,pbImage,annulus=True)
    
        casalog.post("Peak Residual: " + str(residPeak),origin='autobox')
        casalog.post("Residual RMS: " + str(residRMS),origin='autobox')

        # calculate thresholds
        sidelobeThresholdValue = sidelobeThreshold*sidelobeLevel*residPeak
        noiseThresholdValue = noiseThreshold*residRMS

        thresholdNameList = ['sidelobe','noise']
        thresholdList = [sidelobeThresholdValue, noiseThresholdValue]

        maskThreshold = max(thresholdList)
        maskThresholdIdx = thresholdList.index(maskThreshold)

        # report what threshold we're using.
        casalog.post("Using " + thresholdNameList[maskThresholdIdx] + " threshold: " +  str(maskThreshold), origin='autobox')

        # Print out values for all thresholds    
        for (name,value) in zip(thresholdNameList,thresholdList):
            casalog.post( name + " threshold is " + str(value), origin='autobox')

        casalog.post("Creating a new mask",origin='autobox')
        
        # Create a simple threshold mask
        outMask = 'tmp_mask_thresh'+str(imager.ncycle)
        calcThresholdMask(residImage,maskThreshold,outMask)
        
        # If requested, prune regions that are smaller than the beam
        if minBeamFrac > 0:
            casalog.post("pruning regions smaller than " + str(minBeamFrac) + "times the beam size",origin='autobox')
            inMask=outMask
            outMask = 'tmp_mask_prune'+str(imager.ncycle)
            pruneRegions(psfImage,inMask,minBeamFrac,outMask)

        # Smooth mask
        inMask=outMask
        outMask = 'tmp_mask_smooth'+str(imager.ncycle)
        smoothMask(inMask,smoothKernel,outMask)
        
        # Convert smoothed mask to 1's and 0's
        inMask = outMask
        outMask =  'tmp_mask_cut'+str(imager.ncycle)
        cutMask(inMask,cutThreshold,outMask)
            
        # Add masks together if this isn't the first cycle
        if (imager.ncycle > 0):
            inMask = outMask
            outMask = 'tmp_mask_add'+str(imager.ncycle)
            addMasks(maskImage+str(imager.ncycle-1),inMask,outMask)
            casalog.post( 'adding mask '+ maskImage+str(imager.ncycle-1) + ' and ' + inMask,origin='autobox')

        
        #import pdb
        #pdb.set_trace()
            
        # lowThreshold = max(lowNoiseThreshold * residRMS, sidelobeThresholdValue) 
        
        # # creating the constraint mask
        # outConstraintMask = 'tmp_mask_constraint' + str(imager.ncycle)
        # calcThresholdMask(residImage,lowThreshold,outConstraintMask)
        
        # # run a binary dilation on the mask 
        # inMask = maskImage+str(imager.ncycle-1)
        # outMask = 'tmp_mask_growmask'+str(imager.ncycle)
        # growMask(inMask,outConstraintMask,outMask,iterations=100)
        
        # # Smooth the constraint mask
        # inMask=outMask
        # outMask = 'tmp_mask_grow_smooth'+str(imager.ncycle)
        # smoothMask(inMask,smoothKernel,outMask)
        
        # # Convert the smoothed constraint mask to 1's and 0's
        # inMask = outMask
        # outMask =  'tmp_mask_cut'+str(imager.ncycle)
        # cutMask(inMask,cutThreshold,outMask)

        
        # moving masks around to make the right mask available for 
        # the minor/major cycle
        if os.path.exists(maskImage):
            shutil.rmtree(maskImage)
        casalog.post( 'copying ' + outMask + ' to '+ maskImage, origin='autobox')
        shutil.copytree(outMask,maskImage)
        

        # saving the masks and residuals
        casalog.post( 'copying ' +maskImage + ' to '+ maskImage+str(imager.ncycle),origin='autobox')
        shutil.copytree(maskImage,maskImage+str(imager.ncycle))
        shutil.copytree(residImage,residImage+str(imager.ncycle))

        # run a major minor cycle part
        imager.runMinorCycle() 
        imager.runMajorCycle()
    
    # clean up
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

    # get name of current mask
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

    # select only the central channel if it's a cube
    # code is from Remy and needs to be tested.
    #  I might want to do a psf level for each plane of the cube?
    ia.open(psfname)
    mycs=ia.coordsys()
    shape=ia.shape()
    spax=mycs.findaxisbyname("Spectral")
    if shape[spax]>1:        
        psf0=ia.subimage(outfile=psfname+"0",region=rg.frombcs(mycs.torecord(),shape,box="0,0,%i,%i"%(shape[0]-1,shape[1]-1),chans="%i"%int(shape[spax]/2)),overwrite=True)
        mycs.done()
        mycs=psf0.coordsys()
        psf0.done()
        psfname=psfname+"0"
    ia.done()

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
    
    ## add PSF core stuff here??

    # find the maximum sidelobe.
    psfresidimstat=psf.statistics() # Do this in same box as fix?
    psfmin=max(abs(psfresidimstat['min'][0]),abs(psfresidimstat['max'][0]))
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
        
def findResidualStats(residImage, pbImage, annulus=True):
    '''
    calculate the peak of the residual
    '''

    import analyzemsimage

    MADtoRMS =  1.4826 # conversion factor between MAD and RMS. Ref: wikipedia

    ia.open(residImage)
    residStats = ia.statistics(robust=True)
    ia.done()
    residPeak = residStats['max'][0] ## do I want to make this the absolute value of the max/min??
    
    # calculate RMS in annulus
    if annulus:

        imageStats = analyzemsimage.imstatAnnulus(residImage,pbimage=pbImage,innerLevel=0.3,outerLevel=0.2,verbose=True)

        # if nothing in region set RMS to 0.0
        if len(imageStats['medabsdevmed']) > 0:
            residRMS = imageStats['medabsdevmed'][0] * MADtoRMS
        else:
            residRMS = 0.0

    # just calculate the RMS in the whole image
    else:

        residRMS = residStats['medabsdevmed'][0] * MADtoRMS

    casalog.post('Noise in residual image is ' + str(residRMS), origin='autobox')

    return (residPeak, residRMS)
 
def calcThresholdMask(residImage,maskThreshold,outMask):
    '''
    Calculate the mask based on the residual image and mask threshold
    '''

    ia.open(residImage)
    tmpMask = ia.imagecalc(outMask,'iif('+residImage+'>'+str(maskThreshold)+',1.0,0.0)',overwrite=True)
    tmpMask.done()
    ia.done()
   
def smoothMask(maskImage, smoothKernel, outMask):
    '''
    create a smoothed mask
    '''

    # getting the size of the kernel to smooth by
    major = smoothKernel['major']['value']
    minor = smoothKernel['minor']['value']
    pa = smoothKernel['pa']['value']

    # writing out this info to the log
    majorStr = str(major)+smoothKernel['major']['unit']
    minorStr = str(minor)+smoothKernel['minor']['unit']
    paStr = str(pa)+smoothKernel['pa']['unit']

    casalog.post("smoothing by " + majorStr + " by " + minorStr,origin='autobox')

    # convolve the mask by the kernel
    ia.open(maskImage)
    tmp = ia.convolve2d(outfile=outMask,axes=[0,1],type='gauss',
                                       major=majorStr,minor=minorStr,pa=paStr,
                                       overwrite=True)
    tmp.done()
    ia.done()

   
def cutMask(maskImage,cutThreshold, outMask):
    '''
    convert smoothed mask to 1/0
    '''

    ia.open(maskImage)
    maskStats = ia.statistics()
    maskPeak = maskStats['max'][0]

    tmp = ia.imagecalc(outMask,
                       'iif('+maskImage+'>'+str(cutThreshold*maskPeak)+',1.0,0.0)',
                                           overwrite=True)

    ia.done()
    tmp.done()

    
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

def pruneRegions(psfImage,maskImage,minBeamFrac,outMask):

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
    ia.done()

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
    tmp = ia.newimagefromimage(infile=maskImage,outfile=outMask,overwrite=True)
    tmp.putchunk(mask)
    tmp.done()
    ia.done()
    
def growMask(maskImage, constraintMask, outMask, iterations=10):

    ''' 
    grow mask through binary dilation
    '''

    import matplotlib.pyplot as plt
    from scipy.ndimage import binary_dilation, generate_binary_structure
    import numpy as np

    # getting the residual image and the existing mask   
    ia.open(maskImage)
    highMask = ia.getchunk(dropdeg=True)
    ia.close()
    
    ia.open(constraintMask)
    constraintArray = ia.getchunk(dropdeg=True)
    ia.close()
    ia.done()

    # dilating the binary mask into the low contour
    struct = generate_binary_structure(2,1).astype(highMask.dtype)
    growMask = binary_dilation(highMask,structure=struct,iterations=iterations,mask=constraintArray).astype(highMask.dtype)

    # saving the new mask
    tmp = ia.newimagefromimage(infile=maskImage,outfile=outMask,overwrite=True)
    tmp.putchunk(growMask)
    tmp.done()
    ia.done()

        

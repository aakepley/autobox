from taskinit import * # This gets me CASA toolkit tasks.

if casalog.version() > '5.0.0': # fix needed for casa 5.0.0
   ia = iatool() 

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

#Original    
#    from refimagerhelper import PySynthesisImager
#    from refimagerhelper import ImagerParameters, PerformanceMeasure
#

#From Urvashi
    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from imagerhelpers.input_parameters import ImagerParameters
#

    import analyzemsimage
    import pdb

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

    # determine sidelobe level
    sidelobeLevel = fitPSF(psfImage)

    casalog.post("Fractional Sidelobe Level: " + str(sidelobeLevel),origin='autobox')
    
    # calculate size of smoothing beam
    smoothKernel = calcSmooth(psfImage,smoothFactor)
    
    ## Make dirty image
    imager.runMajorCycle()
    
    # make initial threshold mask
    (maskThreshold, lowMaskThreshold) = calcThresholds(residImage,pbImage, sidelobeThreshold, sidelobeLevel,noiseThreshold,lowNoiseThreshold)
    thresholdMask = createThresholdMask(residImage,psfImage,maskThreshold,minBeamFrac,smoothKernel,cutThreshold,ncycle=imager.ncycle)
    outMask = thresholdMask

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

    # Do deconvolution 
    while ( not imager.hasConverged()):
          
        # run a major minor cycle part
        imager.runMinorCycle() 
        imager.runMajorCycle()
    
        # make threshold mask
        (maskThreshold, lowMaskThreshold) = calcThresholds(residImage,pbImage, sidelobeThreshold, sidelobeLevel,noiseThreshold,lowNoiseThreshold)

        thresholdMask = createThresholdMask(residImage,psfImage,maskThreshold,minBeamFrac,smoothKernel,cutThreshold,ncycle=imager.ncycle)

        # figure out if anything was masked in the threshold mask.
        ia.open(thresholdMask)
        maskStats = ia.statistics(axes=[0,1])
        ia.done()
        #doGrow = maskStats['max'] < 1
        doGrow = numpy.ones(len(maskStats['max'])).astype('bool')

        outConstraintMask = 'tmp_mask_constraint'+str(imager.ncycle)
        calcThresholdMask(residImage,lowMaskThreshold,outConstraintMask)
      
        # run a binary dilation on the mask 
        prevMask = maskImage + str(imager.ncycle - 1)
        outMask = 'tmp_mask_grow'+str(imager.ncycle)
        growMask(prevMask,outConstraintMask,outMask,doGrow,iterations=100)
                
        # multiply the binary dilated mask by the constraint
        # mask. Returns lowNoiseThreshold mask that is connected to
        # previous regions. Need this to avoid smoothing things twice.
        inMask = outMask
        outMask = 'tmp_mask_grow_multiply'+str(imager.ncycle)
        multiplyMasks(outConstraintMask,inMask,outMask)

        # prune regions smaller than beam
        if minBeamFrac > 0:
            casalog.post("pruning regions smaller than " + str(minBeamFrac) + "times the beam size",origin='autobox')
            inMask = outMask
            outMask = 'tmp_mask_grow_prune'+str(imager.ncycle)
            pruneRegions(psfImage,inMask,minBeamFrac,outMask)

        # Smooth the resulting mask
        inMask = outMask
        outMask = 'tmp_mask_grow_smooth'+str(imager.ncycle)
        smoothMask(inMask,smoothKernel,outMask)
        
        # Convert smoothed mask to 1's and 0's
        inMask = outMask
        outMask =  'tmp_mask_grow_cut'+str(imager.ncycle)
        cutMask(inMask,cutThreshold,outMask)
        
        # add masks
        previousMask = maskImage + str(imager.ncycle - 1)    
        inMask = outMask
        outMask = 'tmp_mask_prev_add'+str(imager.ncycle)
        addMasks(previousMask,inMask,outMask)
        casalog.post( 'adding mask ' + previousMask + ' and ' + inMask,origin='autobox')

        inMask = outMask
        outMask = 'tmp_mask_threshold_add'+str(imager.ncycle)
        addMasks(thresholdMask,inMask,outMask)
        casalog.post( 'adding mask ' + thresholdMask + ' and ' + inMask,origin='autobox')
            
        # copy mask into final mask
        copyMask(outMask,maskImage)

        # saving the masks and residuals
        casalog.post( 'copying ' +maskImage + ' to '+ maskImage+str(imager.ncycle),origin='autobox')
        shutil.copytree(maskImage,maskImage+str(imager.ncycle))
        shutil.copytree(residImage,residImage+str(imager.ncycle))

    # run a final major/minor cycle
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
     
def calcThresholds(residImage, pbImage, sidelobeThreshold, sidelobeLevel,noiseThreshold,lowNoiseThreshold):
    '''
    calculate the various thresholds for each mask
    '''
   
    import numpy

    # determine peak and RMS of residual image
    (residPeak, residRMS) = findResidualStats(residImage,pbImage,annulus=False)
    
    casalog.post("Peak Residual: " + str(residPeak),origin='autobox')
    casalog.post("Residual RMS: " + str(residRMS),origin='autobox')
    
    # calculate thresholds
    sidelobeThresholdValue = sidelobeThreshold*sidelobeLevel*residPeak
    noiseThresholdValue = noiseThreshold*residRMS
    
    maskThreshold = numpy.maximum(sidelobeThresholdValue,noiseThresholdValue)    
    
    lowNoiseThresholdValue = lowNoiseThreshold*residRMS
    lowMaskThreshold = numpy.maximum(sidelobeThresholdValue,lowNoiseThresholdValue)
    
    # Report what threshold we are using
    for (i,threshold) in zip(numpy.arange(len(maskThreshold)),maskThreshold):
        casalog.post("Channel " + str(i) + " threshold: " + str(threshold),origin='autobox')

    return (maskThreshold, lowMaskThreshold)


def findResidualStats(residImage, pbImage, annulus=True):
    '''
    calculate the peak of the residual
    '''

    import analyzemsimage

    MADtoRMS =  1.4826 # conversion factor between MAD and RMS. Ref: wikipedia

    ia.open(residImage)
    residStats = ia.statistics(robust=True,axes=[0,1])
    ia.done()
    residPeak = residStats['max'] ## do I want to make this the absolute value of the max/min??
    
    # calculate RMS in annulus
    if annulus:

        imageStats = analyzemsimage.imstatAnnulus(residImage,pbimage=pbImage,innerLevel=0.3,outerLevel=0.2,verbose=True,perChannel=True)

        residRMS = imageStats['medabsdevmed'] * MADtoRMS

    # just calculate the RMS in the whole image
    else:

        residRMS = residStats['medabsdevmed'] * MADtoRMS

    casalog.post('Noise in residual image is ' + str(residRMS), origin='autobox')

    return (residPeak, residRMS)
 

def createThresholdMask(residImage,psfImage,maskThreshold,minBeamFrac,smoothKernel,cutThreshold,ncycle=1):
    '''
    Create a threshold mask that is pruned, smoothed, and cut.
    '''

    # Create a simple threshold mask
    outMask = 'tmp_mask_thresh'+str(ncycle)
    calcThresholdMask(residImage,maskThreshold,outMask)
    
    # If requested, prune regions that are smaller than the beam
    if minBeamFrac > 0:
        casalog.post("pruning regions smaller than " + str(minBeamFrac) + "times the beam size",origin='autobox')
        inMask=outMask
        outMask = 'tmp_mask_prune'+str(ncycle)
        pruneRegions(psfImage,inMask,minBeamFrac,outMask)
       
    # Smooth mask
    inMask=outMask
    outMask = 'tmp_mask_smooth'+str(ncycle)
    smoothMask(inMask,smoothKernel,outMask)
        
    # Convert smoothed mask to 1's and 0's
    inMask = outMask
    outMask =  'tmp_mask_cut'+str(ncycle)
    cutMask(inMask,cutThreshold,outMask)

    thresholdMask = outMask

    return thresholdMask
    

def calcThresholdMask(residImage,maskThreshold,outMask):
    '''
    Calculate the mask based on the residual image and mask threshold
    '''

    import numpy as np

    ia.open(residImage)
    resid = ia.getchunk(dropdeg=True)
    ia.done()

    # add another axis if only RA/Dec
    if resid.ndim == 2:
        resid = resid[:,:,np.newaxis]
    
    # create mask
    if resid.ndim == 3:
        mask = np.zeros(resid.shape)
        nchan = (resid.shape)[2]        
        for i in np.arange(nchan):
            mask[:,:,i] = resid[:,:,i] > maskThreshold[i]

        # save image
        tmp = ia.newimagefromimage(infile=residImage,outfile=outMask,overwrite=True)
        
        # add stokes axis
        mask=np.expand_dims(mask,axis=2)
        tmp.putchunk(mask)
        tmp.done()
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

    import numpy as np

    ia.open(maskImage)
    mask = ia.getchunk(dropdeg=True)
    maskStats = ia.statistics(axes=[0,1])
    maskPeak = maskStats['max']
    ia.done()

    # add another axis if only RA/Dec
    if mask.ndim == 2:
        mask = mask[:,:,np.newaxis]

    if mask.ndim == 3:
        newmask = np.zeros(mask.shape)
        nchan = (newmask.shape)[2]
        for i in np.arange(nchan):
            newmask[:,:,i] = mask[:,:,i] > cutThreshold*maskPeak[i]
        
        # save image
        tmp = ia.newimagefromimage(infile=maskImage,outfile=outMask,overwrite=True)
        newmask = np.expand_dims(newmask,axis=2)
        tmp.putchunk(newmask)
        tmp.done()
        ia.done()

    else: 
        casalog.post("Don't know now to deal with an image with this number of dimension",origin='autobox')

    
def addMasks(mask1,mask2,mask3):
    '''

    Add masks together to create final mask. Once the mask is created,
    you need to move the mask in and out of the mask file. Copying new
    images over to the mask file will cause them not to be used.

    '''

    ia.open(mask1)
    ia.close()
    ia.open(mask2)
    ia.close()

    tmp = ia.imagecalc(mask3, 'iif('+mask1+'+'+mask2+'>=1.0,1.0,0.0)', overwrite=True)

    tmp.done()
    ia.done()    

def subtractMasks(mask1,mask2,mask3):
    '''

    Subtract masks together to create final mask. Once the mask is created,
    you need to move the mask in and out of the mask file. Copying new
    images over to the mask file will cause them not to be used.

    Order matters here!!! Mask1 should be the bigger mask.

    '''

    ia.open(mask1)
    ia.close()
    ia.open(mask2)
    ia.close()

    tmp = ia.imagecalc(mask3, 'iif('+mask1+'-'+mask2+'>=1.0,1.0,0.0)', overwrite=True)

    tmp.done()
    ia.done()    


def multiplyMasks(mask1,mask2,mask3):
    '''

    Multiply masks together

    '''

    ia.open(mask1)
    ia.close()
    ia.open(mask2)
    ia.close()

    tmp = ia.imagecalc(mask3, 'iif('+mask1+'*'+mask2+'>=1.0,1.0,0.0)', overwrite=True)

    tmp.done()
    ia.done()    

def pruneRegions(psfImage,maskImage,minBeamFrac,outMask):

    ''' 
    code to prune regions that are too small. Based on code in M100 casaguide
    developed by Crystal Brogan and Jen Donovan Meyer.
    '''

    import scipy.ndimage
    import analysisUtils as au
    import numpy as np

    # get beam area in pixels using analysis utils
    beamArea = au.pixelsPerBeam(psfImage) 
    npix = beamArea * minBeamFrac

    casalog.post("beam area " + str(beamArea) + " pixels",origin='autobox')
    casalog.post("prune area less then " + str(npix) + " pixels",origin='autobox')
    
    # open image and get mask
    ia.open(maskImage)
    mask = ia.getchunk(dropdeg=True)
    ia.done()

    if mask.ndim == 2:
        mask = mask[:,:,np.newaxis]
    
    nchan = (mask.shape)[2]
    
    for k in np.arange(nchan):
        maskPlane = mask[:,:,k]
        # divide the mask up into labeled regions
        labeled, j = scipy.ndimage.label(maskPlane)
        myhistogram = scipy.ndimage.measurements.histogram(labeled,0,j+1,j+1)
        object_slices = scipy.ndimage.find_objects(labeled)

        # go through each region and set the mask to 0 if there aren't enough pixels
        nremoved = 0
        for i in range(j):
            if myhistogram[i+1] < npix:
                maskPlane[object_slices[i]] = 0
                nremoved = nremoved + 1

        mask[:,:,k] = maskPlane

    casalog.post("removed " + str(nremoved) + " regions",origin='autobox')

    # put the mask back
    tmp = ia.newimagefromimage(infile=maskImage,outfile=outMask,overwrite=True)
    mask = np.expand_dims(mask,axis=2)
    tmp.putchunk(mask)
    tmp.done()
    ia.done()
    
def growMask(maskImage, constraintMask, outMask, doGrow, iterations=10):

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
    
    # open the constraintMask
    ia.open(constraintMask)
    constraintArray = ia.getchunk(dropdeg=True)
    ia.close()
    ia.done()

    if highMask.ndim == 2:
        highMask = highMask[:,:,np.newaxis]
        constraintArray = constraintArray[:,:,np.newaxis]

    nchan = (highMask.shape)[2]
    newmask = np.zeros(highMask.shape)

    # dilating the binary mask into the low contour.    
    for k in np.arange(nchan):
        if doGrow[k]:
            struct = generate_binary_structure(2,1).astype(highMask.dtype) 
            newmask[:,:,k] = binary_dilation(highMask[:,:,k],structure=struct,iterations=iterations,mask=constraintArray[:,:,k]).astype(highMask.dtype)

    # saving the new mask
    tmp = ia.newimagefromimage(infile=maskImage,outfile=outMask,overwrite=True)
    # assume that if 3d, we have a degenerate stokes before the spex
    if newmask.ndim > 2:
        newmask=np.expand_dims(newmask,axis=2)
    tmp.putchunk(newmask)
    tmp.done()
    ia.done()   

def copyMask(inMask,outMask):
    '''
    copy mask into the mask file. Need to use getchunk and putchunk to 
    avoid deleting masks
    '''

    ia.open(inMask)
    mask = ia.getchunk()
    ia.close()

    ia.open(outMask)
    ia.putchunk(mask.astype('float'))
    ia.close()
    ia.done()
   

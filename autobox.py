from taskinit import * # This gets me the toolkit tasks.

def runTclean(paramList, sidelobeThreshold,floorThreshold, peakThreshold,
              smoothFactor,cutThreshold):
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
    
    # determine sidelobe level
    sidelobeLevel = fitPSF(psfImage)
    
    ## Make dirty image
    imager.runMajorCycle()
    
    # calculate size of smoothing beam
    smoothKernel = calcSmooth(psfImage,smoothFactor)
    
    ## Do deconvolution and iterations
    while ( not imager.hasConverged() ):

        # determine peak of residual
        residPeak = findResidualPeak(residImage)

        # calculate thresholds
        thresholdNameList = ['Sidelobe','Peak','Floor']
        thresholdList = [sidelobeThreshold*sidelobeLevel*residPeak,
                         peakThreshold*residPeak,
                         floorThreshold]

        ## compare various thresholds -- Do I want an absolute value here?
        maskThreshold = max(thresholdList)
        maskThresholdIdx = thresholdList.index(maskThreshold)

        print "Using " + thresholdNameList[maskThresholdIdx] + " threshold: ", maskThreshold

        # Print out values for all thresholds    
        for (name,value) in zip(thresholdNameList,thresholdList):
            print name, " threshold is ", value
     
        # create a new mask
        calcMask(residImage,maskThreshold,smoothKernel,cutThreshold)

        # move things around so that the mask is read
        if imager.ncycle > 0:
            shutil.copytree(maskImage,maskImage+str(imager.ncycle))

        if os.path.exists(maskImage):
            shutil.rmtree(maskImage)

        shutil.copytree('tmp_final_mask',maskImage)
    
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
        
def findResidualPeak(residImage):
    '''
    calculate the peak of the residual
    '''

    ia.open(residImage)
    residStats = ia.statistics()
    ia.done()
    residPeak = residStats['max'][0] ## do I want to make this the absolute value of the max/min??
    return residPeak

def calcMask(residImage, maskThreshold,smoothKernel,cutThreshold):
    '''
    calculate the mask based on the residual image and the mask threshold
    '''

    ia.open(residImage)
    tmpMask = ia.imagecalc('tmp_mask','iif('+residImage+'>'+str(maskThreshold)+',1.0,0.0)',overwrite=True)

    # this is to avoid a bug with ia.convolve2d. bug reported in CAS-8928. turns out the
    # default is something assigned to the value which is annoying. Need to assign to null 
    # value to get to work.
    #majorStr = str(smoothKernel['major']['value'])+smoothKernel['major']['unit']
    #minorStr = str(smoothKernel['minor']['value'])+smoothKernel['minor']['unit']
    #paStr = str(smoothKernel['pa']['value'])+smoothKernel['pa']['unit']

    tmpSmoothMask = tmpMask.convolve2d(outfile='tmp_smooth_mask',axes=[0,1],type='gauss',
                                     major="",minor="",pa="",beam=smoothKernel,overwrite=True)

    tmpMask.done()
    
    tmpSmoothMaskStats = tmpSmoothMask.statistics()
    tmpSmoothMaskPeak = tmpSmoothMaskStats['max'][0]
    
    tmpFinalMask = tmpSmoothMask.imagecalc('tmp_final_mask',
                                           'iif(tmp_smooth_mask'+'>'+str(cutThreshold*tmpSmoothMaskPeak)+',1.0,0.0)',
                                           overwrite=True)
    
    tmpSmoothMask.done()
    tmpFinalMask.done()
                                        
    ia.done()
  

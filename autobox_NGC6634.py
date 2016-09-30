# Purpose: run autoboxing code with a bunch of different parameters to
# see what works best.

## Data taken from September 7 run.

##/lustre/naasc/sciops/comm/amcnicho/pipeline/root/2015.1.00150.S_2016_09_09T19_56_47.341/SOUS_uid___A001_X2f6_X263/GOUS_uid___A001_X2f6_X264/MOUS_uid___A001_X2f6_X265/working


# tclean(phasecenter='ICRS 17:20:53.3500 -035.47.01.500', scan=['9, 11, 14, 16, 19'], calcres=False, spw=['22:301.223243871~301.22666184GHz;301.369728246~301.375587621GHz;301.426368871~301.431739965GHz;301.503029027~301.508400121GHz;301.667091527~301.67588059GHz;301.750587621~301.755958715GHz;301.799415746~301.802833715GHz;302.421486058~302.426857152GHz;302.648048558~302.651466527GHz;302.955665746~302.959083715GHz,24:302.012543716~302.015961685GHz;302.218598403~302.223969497GHz;302.340668716~302.344086685GHz;302.356293716~302.36557106GHz;302.413422622~302.42416481GHz;302.436860122~302.446137466GHz;302.506684341~302.513520278GHz;303.576020278~303.585785903GHz,26:303.75432342~303.760671076GHz;303.847585139~303.858815607GHz;303.893483576~303.89885467GHz;303.989186701~303.99651092GHz;304.032643732~304.038503107GHz;304.136159357~304.142507014GHz;304.175710139~304.178639826GHz'], vis=['uid___A002_Xaef195_X458c_target.ms'], imagename='uid___A001_X2f6_X265.s29_0.NGC_6334I_sci.spw22_24_26.cont.I.iter1', threshold='0.11054552714Jy', imsize=[240, 240], pbcor=True, npixels=0, calcpsf=False, cell=['0.13arcsec'], outframe='LSRK', gridder='standard', stokes='I', datacolumn='data', savemodel='none', restoration=True, intent='OBSERVE_TARGET#ON_SOURCE', robust=0.5, usemask='user', parallel=False, restart=True, niter=900000, nchan=-1, deconvolver='hogbom', weighting='briggs', mask='/lustre/naasc/sciops/comm/amcnicho/pipeline/root/2015.1.00150.S_2016_09_09T19_56_47.341/SOUS_uid___A001_X2f6_X263/GOUS_uid___A001_X2f6_X264/MOUS_uid___A001_X2f6_X265/working/uid___A001_X2f6_X265.s29_0.NGC_6334I_sci.spw22_24_26.cont.I.iter1.cleanmask', pblimit=0.21925532817840576, restoringbeam='common', specmode='mfs', chanchunks=-1, interactive=0)

import autobox

from refimagerhelper import PySynthesisImager
from refimagerhelper import ImagerParameters, PerformanceMeasure

#####################################################
#### Autoboxing parameters
#####################################################

peakThresholdList = [0.7] # N times peak residual.

# smoothing
smoothFactorList = [3.0] # smooth by this beam size
cutThreshold = 0.5 # cut from peak

# threshold
sidelobeThreshold = 3 # N times sidelobe level
floorThreshold =  0.11054552714*(2.5/4.0) #Jy
#floorThreshold =  0.11054552714 #Jy

# name root
nameRoot = 'NGC6634'

# visibilities
myvislist = ['uid___A002_Xaef195_X458c.ms']


#####################################################
#### Construct ImagerParameters object
#####################################################

for peakThreshold in peakThresholdList:
    for smoothFactor in smoothFactorList:

        currentImageName = nameRoot + '_pT'+str(peakThreshold) + '_smoT' + str(smoothFactor) + '_cutT' + str(cutThreshold)

        print '***************creating ' + currentImageName + '***************'

        imager = None
        paramList = None

        # Put all parameters into dictionaries and check them. 
        paramList = ImagerParameters(
            msname =myvislist,
            field='3',
            spw=['22:301.223243871~301.22666184GHz;301.369728246~301.375587621GHz;301.426368871~301.431739965GHz;301.503029027~301.508400121GHz;301.667091527~301.67588059GHz;301.750587621~301.755958715GHz;301.799415746~301.802833715GHz;302.421486058~302.426857152GHz;302.648048558~302.651466527GHz;302.955665746~302.959083715GHz,24:302.012543716~302.015961685GHz;302.218598403~302.223969497GHz;302.340668716~302.344086685GHz;302.356293716~302.36557106GHz;302.413422622~302.42416481GHz;302.436860122~302.446137466GHz;302.506684341~302.513520278GHz;303.576020278~303.585785903GHz,26:303.75432342~303.760671076GHz;303.847585139~303.858815607GHz;303.893483576~303.89885467GHz;303.989186701~303.99651092GHz;304.032643732~304.038503107GHz;304.136159357~304.142507014GHz;304.175710139~304.178639826GHz'],
            ### Image....
            imagename=currentImageName,
            #### Direction Image Coords
            imsize=[240,240], 
            cell=['0.13arcsec'], 
            phasecenter='ICRS 17:20:53.3500 -035.47.01.500',
            
            ### Spectral Image Coords
            specmode='mfs',
            outframe='LSRK',
    
            gridder='standard',
            
            ### Gridding....
            pblimit=0.21925532817840576,
            
            restart=True,
            
            weighting='briggs',
            robust=0.5,
    
            ### Deconvolution
            niter=600000,
            threshold=str(floorThreshold)+'Jy',
            interactive=False,
    
            deconvolver='hogbom',
            #nterms=2,
            restoringbeam='common',
        
            ### new mask params
            usemask='user'
            #mask=mask,
            #pbmask=pbmask,
            #maskthreshold=maskthreshold,
            #maskresolution=maskresolution,
            #nmask=nmask,
            #autoadjust=autoadjust,
        )




        autobox.runTclean(paramList,sidelobeThreshold=sidelobeThreshold,
                          floorThreshold=floorThreshold,
                          peakThreshold=peakThreshold,
                          smoothFactor=smoothFactor,
                          cutThreshold=cutThreshold,
                          circle=False,save=True)
    


##################################################

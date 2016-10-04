# Purpose: run autoboxing code with a bunch of different parameters to
# see what works best.

## Data taken from September 7 run: 
#/lustre/naasc/sciops/comm/amcnicho/pipeline/root/2013.1.00014.S_2016_09_09T14_50_31.990/SOUS_uid___A001_X12b_X22f/GOUS_uid___A001_X12b_X230/MOUS_uid___A001_X12b_X231/working/pipeline-20160909T150754/html/

# tclean(phasecenter='ICRS 02:42:40.7081 -000.00.47.928', scan=['11, 13, 18, 20', '11, 13, 18, 20'], calcres=False, spw=['17:352.25883406~353.014693435GHz;353.473189529~354.126509842GHz,19:350.59867781~352.20805281GHz,21:340.57230281~341.55667781GHz;341.75980281~342.13480281GHz,23:339.22855281~340.32230281GHz', '17:352.245967119~353.002314775GHz;353.460322588~354.114131181GHz,19:350.57067415~352.19567415GHz,21:340.57554915~341.54429915GHz;341.74742415~342.12242415GHz,23:339.23179915~340.30992415GHz'], vis=['uid___A002_Xa5df2c_X807f_target.ms', 'uid___A002_Xaab3b0_X118d_target.ms'], imagename='uid___A001_X12b_X231.s29_0.NGC_1068_sci.spw17_19_21_23.cont.I.iter1', threshold='0.000372165076315Jy', imsize=[1568, 1568], pbcor=True, npixels=0, calcpsf=False, cell=['0.018arcsec'], outframe='LSRK', gridder='standard', stokes='I', datacolumn='data', savemodel='none', restoration=True, intent='OBSERVE_TARGET#ON_SOURCE', robust=0.5, usemask='user', parallel=False, restart=True, niter=20000000, nchan=-1, deconvolver='hogbom', weighting='briggs', mask='/lustre/naasc/sciops/comm/amcnicho/pipeline/root/2013.1.00014.S_2016_09_09T14_50_31.990/SOUS_uid___A001_X12b_X22f/GOUS_uid___A001_X12b_X230/MOUS_uid___A001_X12b_X231/working/uid___A001_X12b_X231.s29_0.NGC_1068_sci.spw17_19_21_23.cont.I.iter1.cleanmask', pblimit=0.20028750598430634, restoringbeam='common', specmode='mfs', chanchunks=-1, interactive=0)

import autobox

from refimagerhelper import PySynthesisImager
from refimagerhelper import ImagerParameters, PerformanceMeasure

#####################################################
#### Autoboxing parameters
#####################################################

peakThreshold = 0.5 # N times peak residual.

# smoothing 
smoothFactor = 1.0 # smooth by this beam size
cutThreshold = 0.01 # cut from peak

# threshold
sidelobeThreshold = 3 # N times sidelobe level
floorThreshold = 0.000372165076315 *(1.5/4.0) #Jy
noiseThreshold = 5.0 # factor to multiply RMS by for floor threshold.

# pruning
minBeamFrac=0.5

# name root
nameRoot = 'NGC1068'

# visibilities
myvislist = ['uid___A002_Xa5df2c_X807f.ms', 'uid___A002_Xaab3b0_X118d.ms']


#####################################################
#### Construct ImagerParameters object
#####################################################


currentImageName = nameRoot + '_pT'+str(peakThreshold) + '_smoT' + str(smoothFactor) + '_cutT' + str(cutThreshold)

print '***************creating ' + currentImageName + '***************'

imager = None
paramList = None

# Put all parameters into dictionaries and check them. 
paramList = ImagerParameters(
    msname =myvislist,
    field='4',
    spw=['17:352.25883406~353.014693435GHz;353.473189529~354.126509842GHz,19:350.59867781~352.20805281GHz,21:340.57230281~341.55667781GHz;341.75980281~342.13480281GHz,23:339.22855281~340.32230281GHz', '17:352.245967119~353.002314775GHz;353.460322588~354.114131181GHz,19:350.57067415~352.19567415GHz,21:340.57554915~341.54429915GHz;341.74742415~342.12242415GHz,23:339.23179915~340.30992415GHz'],
    ### Image....
    imagename=currentImageName,
    #### Direction Image Coords
    imsize=[1568,1568], 
    cell=['0.018arcsec'], 
    phasecenter='ICRS 02:42:40.7081 -000.00.47.928',
    
    ### Spectral Image Coords
    specmode='mfs',
    outframe='LSRK',
    
    gridder='standard',
    
    ### Gridding....
    pblimit=0.20028750598430634,
    
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
                  noiseThreshold=noiseThreshold,
                  peakThreshold=peakThreshold,
                  smoothFactor=smoothFactor,
                  cutThreshold=cutThreshold,
                  save=True,minBeamFrac=minBeamFrac)
    


##################################################

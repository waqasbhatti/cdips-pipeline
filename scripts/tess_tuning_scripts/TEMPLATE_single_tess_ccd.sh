#!/usr/bin/env bash

##########################################
#
# USAGE: execute from within /pipe-trex/scripts/tess_tuning_scripts/
#
# PURPOSE:
#     Tune the pipeline to get "best" lightcurves (by various metrics).
#
##########################################

############################
# data-specific parameters #
############################
camnum=2
ccdnum=2
projectid=1001
sector='s0001'

#################################
# reduction-specific parameters #
#################################
tuneparameters=true
nworkers=32
aperturelist="1:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0"
epdsmooth=11
epdsigclip=10000
photdisjointradius=2
anetfluxthreshold=50000
anettweak=6
anetradius=30
initccdextent="0:2048,0:2048"
kernelspec="i/2;d=3/2"
catalog_faintrmag=16
fiphotfluxthreshold=500
photreffluxthreshold=500
extractsources=0
binlightcurves=0
translateimages=1
reversesubtract=1
skipepd=1

##########################################
##########################################
##########################################

#########################
# interpret the options #
#########################
# spacecraft configuration map used to produce CAL images.
if [ "$sector" = 's0001' ] ; then
  scid=0120
elif [ "$sector" = 's0002' ] ; then
  scid=0121
elif [ "$sector" = 's0003' ] ; then
  scid=0123
elif [ "$sector" = 's0004' ] ; then
  scid=0124
else
  echo 'error: need to hard-code in scid for sector s000X'
  exit 42
fi

if [ "$binlightcurves" = 1 ] ; then
  binlcoption=binlightcurves
else
  binlcoption=no-binlightcurves
fi

if [ "$translateimages" = 1 ] ; then
  translateoption=translateimages
else
  translateoption=no-translateimages
fi

if [ "$reversesubtract" = 1 ] ; then
  rsuboption=reversesubtract
else
  rsuboption=no-reversesubtract
fi

if [ "$skipepd" = 1 ] ; then
  skipepdoption=skipepd
else
  skipepdoption=no-skipepd
fi

###############################################################################
# define paths. trimmed, single-extension fits images get worked on in fitsdir.
# they match fitsdir+fitsglob. their lightcurvs are written to lcdir.
###############################################################################
if [ "$tuneparameters" = true ] ; then
  tunefullstr='TUNE'
else
  tunefullstr='FULL'
fi

LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/FFI/RED_IMGSUB/"${tunefullstr}
sectordir=$LOCAL_IMGBASE"/"${sector}"/"
fitsdir=$sectordir"RED_"${camnum}"-"${ccdnum}"-"${projectid}"_ISP/"
LOCAL_GLOBPATTERN='tess?????????????-'${sector}'-'${camnum}'-'${ccdnum}'-'${scid}'_cal_img_bkgdsub.fits'
fitsglob=$LOCAL_GLOBPATTERN
lcbase="/nfs/phtess1/ar1/TESS/FFI/LC/"${tunefullstr}
lcsector=$lcbase"/"${sector}"/"
lcdir=${lcsector}"ISP_"${camnum}"-"${ccdnum}"-"${projectid}"/"

####################
# make directories #
####################
if [ ! -d "$LOCAL_IMGBASE" ]; then
  mkdir $LOCAL_IMGBASE
fi
if [ ! -d "$sectordir" ]; then
  mkdir $sectordir
fi
if [ ! -d "$fitsdir" ]; then
  mkdir $fitsdir
fi
if [ ! -d "$lcbase" ]; then
  mkdir $lcbase
fi
if [ ! -d "$lcsector" ]; then
  mkdir $lcsector
fi
if [ ! -d "$lcdir" ]; then
  mkdir $lcdir
fi

################################
# turn images into lightcurves #
################################
logname=${sector}'-cam'${camnum}'-ccd'${ccdnum}'-projid'${projectid}'.log'

cd ../

python -u TESS_reduction.py \
  --projectid $projectid \
  --fitsdir $fitsdir --fitsglob $fitsglob --outdir $fitsdir --field $sector\
  --nworkers $nworkers --aperturelist $aperturelist --lcdirectory $lcdir \
  --convert_to_fitsh_compatible --epdsmooth $epdsmooth \
  --kernelspec $kernelspec \
  --epdsigclip $epdsigclip --photdisjointradius $photdisjointradius \
  --tuneparameters $tuneparameters --anetfluxthreshold $anetfluxthreshold \
  --anettweak $anettweak --initccdextent $initccdextent \
  --anetradius $anetradius --catalog_faintrmag $catalog_faintrmag \
  --fiphotfluxthreshold $fiphotfluxthreshold \
  --photreffluxthreshold $photreffluxthreshold \
  --extractsources $extractsources --$binlcoption \
  --camnum $camnum --ccdnum $ccdnum --$translateoption \
  --$rsuboption --$skipepdoption \
  &> logs/$logname &

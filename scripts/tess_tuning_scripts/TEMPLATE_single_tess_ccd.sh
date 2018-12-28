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
nworkers=20
aperturelist="1.45:7.0:6.0,2.2:7.0:6.0,2.95:7.0:6.0"
epdsmooth=11
epdsigclip=10000
photdisjointradius=2
anetfluxthreshold=50000
anettweak=6
anetradius=30
initccdextent="0:2048,0:2048"
kernelspec="b/4;i/4;d=5/2"
catalog_faintrmag=13
fiphotfluxthreshold=3000
photreffluxthreshold=3000
extractsources=0
binlightcurves=0

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
else
  echo 'error: need to hard-code in scid for sector s000X'
  exit 42
fi

if [ "$binlightcurves" = 1 ] ; then
  binlcoption=binlightcurves
else
  binlcoption=no-binlightcurves
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
LOCAL_GLOBPATTERN='tess?????????????-'${sector}'-'${camnum}'-'${ccdnum}'-'${scid}'_cal_img.fits'
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
  --camnum $camnum --ccdnum $ccdnum \
  &> logs/$logname &

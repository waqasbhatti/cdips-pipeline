#!/usr/bin/env bash

##########################################
#
# USAGE:
#   ./reduce_ete6_field.sh (log piping is automatic)
#
# PURPOSE:
#   make lightcurves from images
#
##########################################

# data-specific parameters
camnum=1
ccdnum=2
projectid=43
orbit='orbit-10'

# reduction-specific parameters
tuneparameters=true
nworkers=20
aperturelist="1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0"
epdsmooth=11    # 11*30min = 5.5 hour median smooth in EPD pre-processing.
epdsigclip=10
photdisjointradius=2
anetfluxthreshold=50000
anettweak=6
anetradius=30
initccdextent="0:2048,0:2048"
kernelspec="b/4;i/4;d=4/4"
catalog_faintrmag=13
fiphotfluxthreshold=3000
photreffluxthreshold=3000
## catalog_faintrmag=16
## fiphotfluxthreshold=300
## photreffluxthreshold=300
extractsources=0
binlightcurves=0

#########################
# interpret the options #
#########################
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

LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/SIMFFI/RED_IMGSUB/"${tunefullstr}
orbitdir=$LOCAL_IMGBASE"/"${orbit}"/"
fitsdir=$orbitdir"RED_"${camnum}"-"${ccdnum}"_ISP/"
LOCAL_GLOBPATTERN='tess?????????????-'${camnum}'-'${ccdnum}'-0016_cal_img.fits'
fitsglob=$LOCAL_GLOBPATTERN
lcbase="/nfs/phtess1/ar1/TESS/SIMFFI/LC/"${tunefullstr}
lcorbit=$lcbase"/"${orbit}"/"
lcdir=${lcorbit}"ISP_"${camnum}"-"${ccdnum}"/"

####################
# make directories #
####################
if [ ! -d "$LOCAL_IMGBASE" ]; then
  mkdir $LOCAL_IMGBASE
fi
if [ ! -d "$orbitdir" ]; then
  mkdir $orbitdir
fi
if [ ! -d "$fitsdir" ]; then
  mkdir $fitsdir
fi
if [ ! -d "$lcbase" ]; then
  mkdir $lcbase
fi
if [ ! -d "$lcorbit" ]; then
  mkdir $lcorbit
fi
if [ ! -d "$lcdir" ]; then
  mkdir $lcdir
fi

################################
# turn images into lightcurves #
################################
logname=${orbit}'-cam'${camnum}'-ccd'${ccdnum}'.log'

python TESS_ETE6_reduction.py \
  --projectid $projectid \
  --fitsdir $fitsdir --fitsglob $fitsglob --outdir $fitsdir --field $orbit\
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
  &> logs/$logname &

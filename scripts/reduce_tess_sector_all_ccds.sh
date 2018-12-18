#!/usr/bin/env bash

##########################################
#
# USAGE:
#   ./reduce_all_ete6_ccds.sh > logs/times.txt &
#
# PURPOSE:
#   FULL reduction.
#
##########################################

# comment-lines for parity w/ reduce_single_ete6_ccd.sh
#
# data-specific parameters
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
catalog_faintrmag=16
fiphotfluxthreshold=300
photreffluxthreshold=300
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
# define paths, make directories.

if [ "$tuneparameters" = true ] ; then
  tunefullstr='TUNE'
else
  tunefullstr='FULL'
fi

for camnum in {1..4}; do
  for ccdnum in {1..4}; do

    # define paths
    LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/SIMFFI/RED_IMGSUB/"${tunefullstr}
    orbitdir=$LOCAL_IMGBASE"/"${orbit}"/"
    fitsdir=$orbitdir"RED_"${camnum}"-"${ccdnum}"_ISP/"
    LOCAL_GLOBPATTERN="tess?????????????-"${camnum}"-"${ccdnum}
    LOCAL_GLOBPATTERN+="-0016_cal_img.fits"
    fitsglob=$LOCAL_GLOBPATTERN
    lcbase="/nfs/phtess1/ar1/TESS/SIMFFI/LC/"${tunefullstr}
    lcorbit=$lcbase"/"${orbit}"/"
    lcdir=${lcorbit}"ISP_"${camnum}"-"${ccdnum}"/"

    # make directories
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

    logname=${orbit}'-cam'${camnum}'-ccd'${ccdnum}'.log'

    echo `date`
    echo "launching reduction for "${fitsdir}
    echo

    python TESS_ETE6_reduction.py \
      --projectid $projectid \
      --fitsdir $fitsdir --fitsglob $fitsglob --outdir $fitsdir --field $orbit\
      --nworkers $nworkers --aperturelist $aperturelist --lcdirectory $lcdir \
      --convert_to_fitsh_compatible --epdsmooth $epdsmooth \
      --epdsigclip $epdsigclip --photdisjointradius $photdisjointradius \
      --tuneparameters $tuneparameters --kernelspec $kernelspec \
      --anetfluxthreshold $anetfluxthreshold --anettweak $anettweak \
      --initccdextent $initccdextent --anetradius $anetradius \
      --catalog_faintrmag $catalog_faintrmag \
      --fiphotfluxthreshold $fiphotfluxthreshold \
      --photreffluxthreshold $photreffluxthreshold \
      --extractsources $extractsources --$binlcoption \
      &> logs/$logname

  done
done

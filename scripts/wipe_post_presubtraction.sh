#!/usr/bin/env bash

# USE WITH CAUTION.
#
# The pre-subtraction steps (source extraction, astrometry, etc) look fine.
# However the post-subtraction was crap, and you need to remove the bad files.

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

###############################################################################
# define paths, make directories.

if [ "$tuneparameters" = true ] ; then
  tunefullstr='TUNE'
else
  tunefullstr='FULL'
fi

for camnum in {1..4}; do
  for ccdnum in {1..4}; do

    # get paths
    LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/SIMFFI/RED_IMGSUB/"${tunefullstr}
    orbitdir=$LOCAL_IMGBASE"/"${orbit}"/"
    fitsdir=$orbitdir"RED_"${camnum}"-"${ccdnum}"_ISP/"
    LOCAL_GLOBPATTERN="tess?????????????-"${camnum}"-"${ccdnum}
    LOCAL_GLOBPATTERN+="-0016_cal_img.fits"
    fitsglob=$LOCAL_GLOBPATTERN
    lcbase="/nfs/phtess1/ar1/TESS/SIMFFI/LC/"${tunefullstr}
    lcorbit=$lcbase"/"${orbit}"/"
    lcdir=${lcorbit}"ISP_"${camnum}"-"${ccdnum}"/"

    echo "removing post-presubtraction files from "${fitsdir}

    # from imagesubphot
    rm ${fitsdir}*XTRNS*jpg
    rm ${fitsdir}*xtrns.fits
    rm ${fitsdir}*.itrans
    rm ${fitsdir}*.xysdk
    rm ${fitsdir}*ASTOMREF*jpg
    rm ${fitsdir}rsub-*
    rm ${fitsdir}JPEG-SUBTRACTEDCONV*jpg

    echo "removing lightcurves from "${lcdir}

    # from lightcurves. takes a while (millions of files) -> run in "parallel"
    echo ${lcdir}* | xargs rm -rf &

  done
done

# remove all reference files (otherwise bad old cached files get collected)
echo "removing all reference files"
rm /nfs/phtess1/ar1/TESS/SIMFFI/BASE/reference-frames/*
rm -rf /nfs/phtess1/ar1/TESS/SIMFFI/BASE/frameinfo-cache/*

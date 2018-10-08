#!/usr/bin/env bash

##########################################
#
# USAGE:
#   ./reduce_ete6_field.sh &> logs/log1.txt &
#
# PURPOSE:
#   make lightcurves from images
#
##########################################

# data-specific parameters
camnum=1
ccdnum=3
projectid=42
field='sector1'

# reduction-specific parameters
nworkers=20
aperturelist="1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0"
epdsmooth=11    # 5.5 hour median smoothing in EPD pre-processing.
epdsigclip=10
photdisjointradius=2

# paths
LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/SIMFFI/" # tess ete6 data is here
fitsdir=$LOCAL_IMGBASE"RED_"${camnum}"-"${ccdnum}"_ISP/"
LOCAL_GLOBPATTERN='tess?????????????-'${camnum}'-'${ccdnum}'-0016_cal_img.fits'
fitsglob=$LOCAL_GLOBPATTERN
lcdir=$LOCAL_IMGBASE"LC/ISP_"${camnum}"-"${ccdnum}"/"

# make some lightcurves
python TESS_ETE6_reduction.py \
  --projectid $projectid \
  --fitsdir $fitsdir --fitsglob $fitsglob --outdir $fitsdir --field $field\
  --nworkers $nworkers --aperturelist $aperturelist --lcdirectory $lcdir \
  --convert_to_fitsh_compatible --epdsmooth $epdsmooth \
  --epdsigclip $epdsigclip --photdisjointradius $photdisjointradius

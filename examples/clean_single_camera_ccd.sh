#!/usr/bin/env bash

# PURPOSE: USE WITH CAUTION.
#
# remove post-presubtraction files:
#   (XTRNS, .itrans, .xysdk, astromref, rsub*, subtractedconv)
# remove select lightcurves
# remove select reference images
# clear out astromrefs and photrefs table.
# clear out calibratedframes table.
# clear frameinfo cache from past day.

camnum=2
ccdnum=4

##########

psql -U hpx -h xphtess1 hpx -c \
  "delete from astromrefs where camera = "${camnum}" and ccd = "${ccdnum}";"

psql -U hpx -h xphtess1 hpx -c \
  "delete from photrefs where camera = "${camnum}" and ccd = "${ccdnum}";"

##########

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

# define paths to remove

if [ "$tuneparameters" = true ] ; then
  tunefullstr='TUNE'
else
  tunefullstr='FULL'
fi

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

# remove specific reference files (otherwise bad old cached files get collected)
echo "removing all reference files"

rm /nfs/phtess1/ar1/TESS/SIMFFI/BASE/reference-frames/*proj${projectid}-camera${camnum}-ccd${ccdnum}*

# for the frameinfo-cache, the name is tricky. But if you're running this, it's
# probably safe to remove everything from the past day.
echo "removing frameinfo cache from past day"

rm -rf `find /nfs/phtess1/ar1/TESS/SIMFFI/BASE/frameinfo-cache/* -mtime -1 -print`

# clean calibratedframes. regex reference:
# https://www.postgresql.org/docs/9.5/static/functions-matching.html

psql -U hpx -h xphtess1 hpx -c \
  "DELEte from calibratedframes where fits like '"${fitsdir}"%';"

#!/usr/bin/env bash

# PURPOSE: USE WITH CAUTION.
#
# remove post-presubtraction files:
#   (XTRNS, .itrans, .xysdk, astromref, rsub*, subtractedconv)
# remove select lightcurves
# remove select reference images
# clear out astromrefs and photrefs table.
# clear out calibratedframes table.
# clear frameinfo cache from past week. 
camnum=2
ccdnum=1
projectid=1320
sector='s0003'
tuneparameters=false

##########

psql -U hpx -h xphtess1 hpx -c \
  "delete from astromrefs where projectid = "${projectid}" and camera = "${camnum}" and ccd = "${ccdnum}";"

psql -U hpx -h xphtess1 hpx -c \
  "delete from photrefs where projectid = "${projectid}" and camera = "${camnum}" and ccd = "${ccdnum}";"

##########

# define paths to remove

if [ "$tuneparameters" = true ] ; then
  tunefullstr='TUNE'
else
  tunefullstr='FULL'
fi

# get paths
LOCAL_IMGBASE="/nfs/phtess1/ar1/TESS/FFI/RED_IMGSUB/"${tunefullstr}
sectordir=$LOCAL_IMGBASE"/"${sector}"/"
fitsdir=$sectordir"RED_"${camnum}"-"${ccdnum}"-"${projectid}"_ISP/"
LOCAL_GLOBPATTERN='tess?????????????-'${sector}'-'${camnum}'-'${ccdnum}'-'${scid}'_cal_img_bkgdsub.fits'
fitsglob=$LOCAL_GLOBPATTERN
lcbase="/nfs/phtess1/ar1/TESS/FFI/LC/"${tunefullstr}
lcsector=$lcbase"/"${sector}"/"
lcdir=${lcsector}"ISP_"${camnum}"-"${ccdnum}"-"${projectid}"/"
lcdirold=${lcsector}"ISP_"${camnum}"-"${ccdnum}"-"${projectid}"_old/"

echo "removing post-presubtraction files from "${fitsdir}

# from imagesubphot
rm ${fitsdir}*XTRNS*jpg
rm ${fitsdir}*xtrns.fits
rm ${fitsdir}*.itrans
rm ${fitsdir}*.xysdk
rm ${fitsdir}*ASTOMREF*jpg
rm ${fitsdir}rsub-*
rm ${fitsdir}JPEG-SUBTRACTEDCONV*jpg

# remove lightcurves and stats. since this takes a while... first rename, then
# run!
echo "moving "${lcdir}" to "${lcdirold}
mv ${lcdir} ${lcdirold}
echo "removing lightcurves from "${lcdirold}
rm -rf ${lcdirold}

#FIXME: old
# remove specific reference files (otherwise bad old cached files get collected)
echo "removing all reference files"

rm /nfs/phtess1/ar1/TESS/FFI/BASE/reference-frames/*proj${projectid}-camera${camnum}-ccd${ccdnum}*
# 
# for the frameinfo-cache, the name is tricky. But if you're running this, it's
# probably safe to remove everything from the past week.
echo "removing frameinfo cache from past week"

rm -rf `find /nfs/phtess1/ar1/TESS/FFI/BASE/frameinfo-cache/* -mtime -7 -print`

# clean calibratedframes. regex reference:
# https://www.postgresql.org/docs/9.5/static/functions-matching.html

psql -U hpx -h xphtess1 hpx -c \
  "DELETE from calibratedframes where fits like '"${fitsdir}"%';"

# # clean cache
# echo "removing all reference files"
# rm /nfs/phtess1/ar1/TESS/FFI/BASE/reference-frames/*cam${camnum}-ccd${ccdnum}*
# rm /nfs/phtess1/ar1/TESS/FFI/BASE/reference-frames/*camera${camnum}-ccd${ccdnum}*
# rm -rf /nfs/phtess1/ar1/TESS/FFI/BASE/frameinfo-cache/*

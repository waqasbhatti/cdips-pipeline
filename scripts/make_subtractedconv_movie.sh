#!/usr/bin/env bash

# PURPOSE: once SUBTRACTEDCONV jpeg images exist, collect them into a movie.

# 13 days in ~1 minute makes ~12 frames per second

# individual ccd & camera
orbitnum=10
ccdnum=4
camnum=1

fdir='/nfs/phtess1/ar1/TESS/SIMFFI/RED_IMGSUB/TUNE/orbit-'${orbitnum}'/'${camnum}'-'${ccdnum}'_ISP/'
outdir='/nfs/phtess1/ar1/TESS/SIMFFI/MOVIES/'

# single ccd, single camera
ffmpeg -framerate 24 \
       -pattern_type \
       glob -i ${fdir}'JPEG-SUBTRACTEDCONV-rsub*-*-'${camnum}'-'${ccdnum}'-0016_*.jpg' \
       -c:v libx264 \
       -preset fast \
       ${outdir}ete6_orbit${orbitnum}_cam${camnum}-ccd${ccdnum}.mp4

## # a movie for each camera, merged ccds
## for camnum in {1..4}; do
##     ffmpeg -framerate 24 \
##            -pattern_type \
##            glob -i '../data/orbit-10_cal_ffi_jpg_logscale/tess*-*-'${camnum}'-*_fullcam.cal.jpg' \
##            -c:v libx264 \
##            -preset fast \
##            ../movies/orbit-10_cam${camnum}_fullcam_logscale_cal_ffi.mp4
## done

## # a movie for each camera and ccd
## for camnum in {1..4}; do
##     for ccdnum in {1..4}; do
##         ffmpeg -framerate 24 \
##                -pattern_type \
##                glob -i '../data/orbit-10_cal_ffi_jpg/tess*-*-'${camnum}'-*_ccd'${ccdnum}'*.jpg' \
##                -c:v libx264 \
##                -preset fast \
##                ../movies/orbit-10_cam${camnum}-ccd${ccdnum}_cal_ffi.mp4
##     done
## done

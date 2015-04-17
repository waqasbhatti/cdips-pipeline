#!/usr/bin/env python

'''
imagesubphot.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - March 2015

This contains functions to do image subtraction photometry.

'''

#############
## IMPORTS ##
#############

import os
import os.path
import glob
import multiprocessing as mp
import subprocess
import shlex
from datetime import datetime
import re
import json
import shutil
import random
import cPickle as pickle

import numpy as np

from scipy.spatial import cKDTree as kdtree
from scipy.signal import medfilt
from scipy.linalg import lstsq
from scipy.stats import sigmaclip as stats_sigmaclip
from scipy.optimize import curve_fit

import scipy.stats
import numpy.random as nprand

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import pyfits

import imageutils
from imageutils import get_header_keyword


# set this to show extra info
DEBUG = False

# CCD minimum and maximum X,Y pixel coordinates
# used to strip things outside FOV from output of make_frame_sourcelist
CCDEXTENT = {'x':[0.0,2048.0],
             'y':[0.0,2048.0]}

# zeropoint mags for the HATPI lenses given exp time of 30 seconds
# from Chelsea's src directory on phs3: run_phot_astrom.py (2014-12-15)
# FIXME: check where these came from and fix if out of date, especially if
# cameras moved around
ZEROPOINTS = {5:17.11,
              6:17.11,
              7:17.11,
              8:16.63}

# used to get the station ID, frame number, and CCD number from a FITS filename
FRAMEREGEX = re.compile(r'(\d{1})\-(\d{6}\w{0,1})_(\d{1})')

# command string to do a 2massread for a specified FOV
TWOMASSREADCMD = ("2massread -r {ra:f} -d {dec:f} -s {boxlen:f} "
                  "--cat {catalogpath} -mr {brightrmag:f} -Mr {faintrmag:f} "
                  "--xieta-coords {ra:f} {dec:f} -o {outfile}")

# command string to do a 2massread for a specified FOV
UCAC4READCMD = ("ucac4read -r {ra:f} -d {dec:f} -s {boxlen:f} "
                  "--cat {catalogpath} -mr {brightrmag:f} -Mr {faintrmag:f} "
                  "-o {outfile}")

# locations of catalogs
CATALOGS = {
    '2MASS':{'cmd':TWOMASSREADCMD,
             'path':'/nfs/lcohpsrv1/ar0/P/HP0/CAT/2MASS/2MASS_JH_AP/data'},
    'UCAC4':{'cmd':UCAC4READCMD,
             'path':'/nfs/lcohpsrv1/ar0/P/HP0/CAT/UCAC4'}
    }

# command string to run fistar
# parameters:
# {frame}
# {extractedlist}
# {ccdgain}
# {zeropoint}
# {exptime}
# {fluxthreshold}
# {ccdsection}
FISTARCMD = ("{fistarexec} -i {frame} -o {extractedlist} "
             "--model elliptic --iterations symmetric=4,general=2 "
             "--algorithm uplink --format id,x,y,bg,amp,s,d,k,flux,s/n "
             "-g {ccdgain} --mag-flux {zeropoint},{exptime} --sort flux "
             "--flux-threshold {fluxthreshold} --section {ccdsection} "
             "--comment")

# command string to run transformation from RA/Dec in 2MASS catalogs to X/Y in
# FITs image to eventually run photometry at those locations using the
# transformations noted in each frame's .wcs file. we do the transform, then
# remove any objects outside the [0:2048, 0:2048] box for the CCD. we then use
# the resulting source list as input to fiphot. the parameters are:

# {transformer}:       'anrd2xy' or similar, coord -> pix converter executable
# {framewcsfile}:      WCS transformation file associated with frame
# {catalogsourcelist}: 2MASS catalog file associated with camera FOV
# {outfile}:           temporary output file
TRANSFORMCMD = ("{transformer} -w {framewcsfile} "
                "-o {outputfile} "
                "-c 2,3 "
                "{catalogsourcelist} ")

# fiphot command string to run fiphot. requires a sourcelist obtained from
# running TRANSFORMCMD and removing objects outside the CCD. the parameters are:
# FIXME: WTF is single-background and why is it 3?

# {fits}:         name of input FITS frame
# {sourcelist}:   name of the source list file
# {zeropoint}:    zeropoint magnitude from ZEROPOINTS above
# {xycols}:       comma-separated 1-indexed column numbers of x and y coords
# {ccdgain}:      gain of the CCD
# {ccdexptime}:   exposure time of the CCD in seconds
# {aperturelist}: aperture list in the following format (all units are pixels):
#                 aper1rad:sky1inner:sky1rad,...,aperNrad,skyNinner,skyNrad
# {fitsbase}:     the FITS base filename without any extensions
#                 e.g. for 1-377741e_5.fits, this is 1-377741e_5
#                 (used for later collection of fiphot files into an LC)
# {outfile}:      name of the output .phot file (binary format)
# the imagesub fiphot aperture string is:
# --apertures '2.5:7.0:6.0,2.95:7.0:6.0,3.35:7.0:6.0,3.95:7.0:6.0'
# but we should use 3 apertures, so this will become:
# --apertures '2.95:7.0:6.0,3.35:7.0:6.0,3.95:7.0:6.0'

FIPHOTCMD = ("fiphot --input {fits} --input-list {sourcelist} "
             "--col-id 1 --col-xy {xycols} --gain {ccdgain:f} "
             "--mag-flux {zeropoint:f},{ccdexptime:f} "
             "--apertures {aperturelist} "
             "--sky-fit 'mode,sigma=3,iterations=2' --disjoint-radius 2 "
             "--serial {fitsbase} "
             "--format 'ISXY,BbMms' --nan-string 'NaN' "
             "--aperture-mask-ignore 'saturated' --comment '--comment' "
             "--single-background 3 {binaryout} --output {outfile} -k")


## list of tasks

# 1. run fistar on everything, run astrometry on everything

# 2. find a good astrometry reference (sharp stars, low zenith distance, new
    # moon)

# 3. run grmatch to get shifts between each frame




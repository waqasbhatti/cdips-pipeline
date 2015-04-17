#!/usr/bin/env python

'''astrometry.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - March 2015

This contains functions to perform astrometry on reduced HATPI frames using
astrometry.net.

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


########################
## USEFUL DEFINITIONS ##
########################

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



###############################
## ANET ASTROMETRY FUNCTIONS ##
###############################

def anet_solve_frame(srclist,
                     wcsout,
                     ra,
                     dec,
                     width=13,
                     tweak=6,
                     radius=4,
                     cols=(2,3)):
    '''
    This uses anet to solve a frame by using its extracted sources and returns a
    .wcs file containing the astrometric transformation between frame x,y and
    RA/DEC.

    Example anet command:

    anet --ra 60. --dec -22.5 --radius 4 --width 13 --tweak 6 --cols 2,3 1-383272f_5.fistar --wcs 1-383272f_5.wcs

    assuming an ~/.anetrc file with the following contents:

    xsize = 2048
    ysize = 2048
    tweak = 3
    xcol = 2
    ycol = 3
    verify = 1
    log = 1
    indexpath = /P/HP0/CAT/ANET_INDEX/ucac4_2014/

    otherwise, we'll need to provide these as kwargs to the anet executable.

    The input sourcelist can come from fistar, with a fluxthreshold set 10000 to
    just get the bright stars. This makes anet faster.

    '''

    ANETCMDSTR = ("anet -r {ra} -d {dec} -w {width} "
                  "--tweak {tweak} --radius {radius} "
                  "--cols {colx},{coly} --wcs {outwcsfile} {sourcelist}")


    anetcmd = ANETCMDSTR.format(ra=ra,
                                dec=dec,
                                width=width,
                                tweak=tweak,
                                radius=radius,
                                colx=cols[0],
                                coly=cols[1],
                                outwcsfile=wcsout,
                                sourcelist=srclist)

    if DEBUG:
        print(anetcmd)

    # execute the anet shell command
    anetproc = subprocess.Popen(shlex.split(anetcmd),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

    # get results
    anet_stdout, anet_stderr = anetproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if anetproc.returncode == 0:

        print('%sZ: anet WCS %s generated for frame sourcelist %s' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(wcsout), os.path.abspath(srclist)))

        return os.path.abspath(wcsout)

    else:

        print('%sZ: anet WCS %s failed for frame sourcelist %s! Error was: %s' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(wcsout),
               os.path.abspath(srclist),
               anet_stderr))

        return None


def parallel_anet_worker(task):
    '''
    This expands the task arg into the args and kwargs necessary for
    extract_frame_sources.

    '''

    return (task[0], anet_solve_frame(task[0],task[1], task[2], task[3],
                                      **task[4]))



def parallel_anet(srclistdir,
                  outdir,
                  ra, dec,
                  nworkers=16,
                  maxtasksperworker=1000,
                  width=13,
                  tweak=6,
                  radius=4,
                  cols=(2,3)):
    '''
    This does parallel anet astrometry for all frames in srclistdir and
    generates their wcs files.

    '''

    # get a list of all fits files in the directory
    fistarlist = glob.glob(os.path.join(srclistdir,'*_?.fistar'))

    print('%sZ: found %s fistar files in %s, starting astrometry...' %
          (datetime.utcnow().isoformat(),
           len(fistarlist), srclistdir))

    if outdir and not os.path.exists(outdir):

        print('%sZ: making new output directory %s' %
              (datetime.utcnow().isoformat(),
               outdir))
        os.mkdir(outdir)

    pool = mp.Pool(nworkers, maxtasksperchild=maxtasksperworker)

    tasks = [
        [x, os.path.join(outdir, os.path.basename(x.replace('.fistar','.wcs'))),
                         ra, dec, {'width':width,
                                   'tweak':tweak,
                                   'radius':radius,
                                   'cols':cols}]
        for x in fistarlist
        ]

    # fire up the pool of workers
    results = pool.map(parallel_anet_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    # this is the return dictionary
    returndict = {x:y for (x,y) in results}
    return returndict




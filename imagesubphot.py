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


#################
## DEFINITIONS ##
#################

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


##################################
## ASTROMETRIC REFERENCE FRAMES ##
##################################

def select_astromref_frame(fitsdir,
                           fitsglob,
                           srclistdir=None,
                           srclistext='.fistar',
                           photdir=None,
                           photext='.fiphot'):
    '''
    This picks an astrometric reference frame.

    We're looking for (in order):

    - highest median S (smallest FWHM)
    - median D value closest to zero (roundest stars)
    - lowest median background
    - largest number of sources with good extractions


    '''

    # first, get the frames
    fitslist = glob.glob(os.path.join(fitsdir, fitsglob))

    if not srclistdir:
        srclistdir = fitsdir
    if not photdir:
        photdir = fitsdir

    print('%sZ: %s FITS files found in %s matching glob %s, '
          'finding photometry and source lists...' %
          (datetime.utcnow().isoformat(),
           len(fitslist), fitsdir, fitsglob))

    goodframes = []
    goodphots = []
    goodsrclists = []

    # associate the frames with their fiphot files
    for fits in fitslist:

        photpath = os.path.join(
            photdir,
            os.path.basename(fits).strip('.fits.fz') + photext
            )
        srclistpath = os.path.join(
            srclistdir,
            os.path.basename(fits).strip('.fits.fz') + srclistext
            )

        if os.path.exists(photpath) and os.path.exists(srclistpath):
            goodframes.append(fits)
            goodphots.append(photpath)
            goodsrclists.append(srlistpath)


    # we only work on goodframes now
    print('%sZ: selecting an astrometric reference frame...' %
          (datetime.utcnow().isoformat(),))


    median_sval = []
    median_dval = []
    median_background = []
    good_detections = []

    # go through all the frames and find their properties
    for frame, phot, srclist in zip(goodframes, goodphots, goodsrclists):

        # decide if the phot file is binary or not. read the first 600
        # bytes and look for the '--binary-output' text
        with open(phot,'rb') as photf:
            header = photf.read(600)

        if '--binary-output' in header and HAVEBINPHOT:

            photdata_f = read_fiphot(phot)
            photdata = {
                'mag':np.array(photdata_f['per aperture'][2]['mag']),
                'err':np.array(photdata_f['per aperture'][2]['mag err']),
                'flag':np.array(
                    photdata_f['per aperture'][2]['status flag']
                    )
                }
            del photdata_f

        elif '--binary-output' in header and not HAVEBINPHOT:

            print('%sZ: %s is a binary phot file, '
                  'but no binary phot reader is present, skipping...' %
                  (datetime.utcnow().isoformat(), phot))
            continue

        else:

            # read in the phot file
            photdata = np.genfromtxt(
                phot,
                usecols=(12,13,14),
                dtype='f8,f8,S5',
                names=['mag','err','flag']
                )

        # now, get the data from the associated fistar file
        srcdata = np.genfromtxt(srclist,
                                usecols=(3,5,6),
                                dtype='f8,f8,f8',
                                names=['background',
                                       'svalue',
                                       'dvalue'])

        # number of good detections
        good_detind = (photdata['flag'] == 'G') | (photdata['flag'] == '0')
        good_detections.append(len(photdata['mags'][good_detind]))

        # median background, d, and s
        median_background.append(np.nanmedian(srcdata['background']))
        median_dval.append(np.nanmedian(srcdata['dvalue']))
        median_sval.append(np.nanmedian(srcdata['svalue']))

    # now find the best astrometric reference frame

    # to np.arrays first
    median_sval = np.array(median_sval)
    median_dval = np.array(median_dval)
    median_background = np.array(median_background)
    good_detections = np.array(good_detections)

    # get the best S --> largest S at the top
    median_sval_ind = np.argsort(median_sval)[::-1]

    # here, we want the values closest to zero to be at the top
    median_dval_ind = np.argsort(np.fabs(median_dval))

    # want the smallest background
    median_background_ind = np.argsort(median_background)

    # and the most number of detections
    good_detections_ind = np.argsort(good_detections)[::-1]

    # get the top 200 of each index
    median_sval_ind = median_sval_ind[:200]
    median_dval_ind = median_dval_ind[:200]
    median_background_ind = median_background_ind[:200]
    good_detections_ind = good_detections_ind[:200]

    # now intersect all of these arrays to find the best candidates for the
    # astrometric reference frame
    best_frame_ind = np.intersect1d(
        np.intersect1d(median_sval_ind,
                       median_dval_ind,
                       assume_unique=True),
        np.intersect1d(median_background_ind,
                       good_detections_ind,
                       assume_unique=True),
        assume_unique=True
        )

    if len(best_frame_ind) > 0:
        goodframes = np.array(goodframes)
        selectedreference = goodframes[best_frame_ind[0]]

        print('%sZ: selected best astrometric reference frame is %s' %
              (datetime.utcnow().isoformat(), selectedreference))

        return selectedreference

    else:

        print('ERR! %sZ: could not select a good astrometric reference frame!' %
              (datetime.utcnow().isoformat(), ))

        return



def get_smoothed_xysdk_coeffs(fistardir):
    '''
    This generates smoothed xy and sdk coefficents for use with iphot later
    (these go into the output photometry file).

    grtrans ${APPHOT}/$base.${EXT_FISTAR} --col-xy 2,3 --col-fit 6,7,8 \
                  --col-weight 10 --order 4 \
                  --iterations 3 --rejection-level 3 \
                  --comment --output-transformation ${IPHOT}/$base.xysdk

    '''


def get_astromref_shifts():
    '''
    This gets shifts between the astrometric reference frame and all other
    frames.

    '''


def transform_frames_to_astromref():
    '''
    This shifts all frames to the astrometric reference.

    '''


##################################
## PHOTOMETRIC REFERENCE FRAMES ##
##################################

def select_photref_frames(photdir,
                          minframes=80):
    '''
    This selects a group of photometric reference frames that will later be
    stacked and medianed to form the single photometric reference frame.

    Should be similar to the aperturephot version.

    - best median scatter of photometry
    - lowest median error in photometry
    - lowest median background measurement
    - low zenith distance
    - high moon and sun distance
    - large number of stars detected
    - fewest number of failed source photometry extractions

    '''



def stack_photref_frames(framelist):
    '''
    This stacks the photometric reference frames in median mode.

    - first, transform all of them to the astrometric reference frame

    - then use fiarith in median mode to stack them

    - (or actually, use scipy to do this instead, making sure to copy over the
      first frame's header + adding a new keyword listing all the component
      frames)

    '''


def photometry_on_stacked_photref(
    stackedphotref,
    apertures='2.95:7.0:6.0,3.35:7.0:6.0,3.95:7.0:6.0'
    ):
    '''
    This runs fiphot in the special iphot mode on the stacked photometric
    reference frame. See cmrawphot.sh for the correct commandline to use.

    '''

def generate_photref_registration(stackedphotref):
    '''
    This generates a registration file for use with convolution later. The
    convolution and subtraction step below needs this.

    '''


#################################
## CONVOLUTION AND SUBTRACTION ##
#################################

def convolve_and_subtract_frames(transformedframelist,
                                 stackedphotref,
                                 stackedphotrefregistration,
                                 kernelspec='b/4;i/4;d=4/4'):
    '''
    This convolves the photometric reference to each frame, using the specified
    kernel, then subtracts the frame from the photometric reference to produce
    the subtracted frames.

    '''

def photometry_on_subtracted_frames(subtractedframelist,
                                    photrefrawphot):
    '''
    This runs photometry on the subtracted frames and finally produces the ISM
    magnitudes.

    See run_iphot.py and IMG-3-PHOT_st5.sh for what this is supposed to do.

    '''


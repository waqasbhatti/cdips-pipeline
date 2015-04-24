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

def select_astrometric_reference_frame(photdir):
    '''
    This picks an astrometric reference frame.

    We're looking for:

    - highest median S (smallest FWHM)
    - lowest median background
    - smallest Z

    TODO:


    '''


def get_smoothed_xysdk_coeffs(fistardir):
    '''
    This generates smoothed xy and sdk coefficents for use with iphot later
    (these go into the output photometry file).

    grtrans ${APPHOT}/$base.${EXT_FISTAR} --col-xy 2,3 --col-fit 6,7,8 \
                  --col-weight 10 --order 4 \
                  --iterations 3 --rejection-level 3 \
                  --comment --output-transformation ${IPHOT}/$base.xysdk

    '''


def get_reference_frame_shifts():
    '''
    This gets shifts between the astrometric reference frame and all other
    frames.

    '''


def transform_frames_to_astrometric_reference():
    '''
    This shifts all frames to the astrometric reference.

    '''


##################################
## PHOTOMETRIC REFERENCE FRAMES ##
##################################

def select_photometric_reference_frames(photdir,
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



def stack_photometric_reference_frames(framelist):
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
    This generates a registration file for use with convolution later.

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


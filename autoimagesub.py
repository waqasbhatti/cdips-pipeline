#!/usr/bin/env python

'''autoimagesub.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - 09/2016

This contains functions to run image subtraction photometry continuously.

TODO:

- functions for running imagesub steps on lists of reduced FITS
- functions that search for reference images
- functions that search for reference image photometry catalogs
- function that knows how to time-out a subprocess call to ficonv

'''

#############
## IMPORTS ##
#############

import os
import os.path
import glob
import multiprocessing as mp
import subprocess
from subprocess import check_output
import shlex
from datetime import datetime
import re
import json
import shutil
import random
import cPickle as pickle


import aperturephot as ap
import imagesubphot as ism


############
## CONFIG ##
############

DEBUG = False

# used to get the station ID, frame number, and CCD number from a FITS filename
FRAMEREGEX = re.compile(r'(\d{1})\-(\d{6}\w{0,1})_(\d{1})')

# this defines the field string and CCDs
FIELD_REGEX = re.compile('^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$')
FIELD_CCDS = [5,6,7,8]

# defines where the reference frames go
REFBASEDIR = '/P/HP0/BASE/reference-frames'
REFINFO = os.path.join(FIELD_REFBASEDIR,'refinfo.sqlite')


######################
## REFERENCE FRAMES ##
######################

def generate_astromref(fitsfiles,
                       field=None,
                       ccd=None,
                       projectid=None,
                       refdir=REFBASEDIR,
                       refinfo=REFINFO):
    '''
    This chooses an astrometry reference frame from the frames in fitfiles.

    writes the frame to refdir.

    ref frames have the following filename pattern:

    proj{projectid}-ccd{ccd}-{field}-astromref.fits

    '''




def find_astromref(field, refinfo=REFINFO):
    '''
    This finds the reference frame for the field.

    reference frames are found in:

    {REFBASEDIR}/{field}



    '''

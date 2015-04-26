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

# get fiphot binary reader
try:
    from HATpipepy.Common.BinPhot import read_fiphot
    HAVEBINPHOT = True
except:
    print("can't import binary fiphot reading functions from "
          "HATpipe, binary fiphot files will be unreadable!")
    HAVEBINPHOT = False


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


#######################
## COMMAND TEMPLATES ##
#######################

XYSDKCMD = ('grtrans {fistarfile} --col-xy 2,3 --col-fit 6,7,8 '
            '--col-weight 10 --order 4 '
            '--iterations 3 --rejection-level 3 '
            '--comment --output-transformation {xysdkfile}')

FRAMESHIFTCALCCMD = (
    'grmatch --match-points '
    '-r {astromref} --col-ref 2,3 --col-ref-ordering +9 '
    '-i {fistartoshift} --col-inp 2,3 --col-inp-ordering +9 '
    '--weight reference,column=9 '
    '--triangulation maxinp=5000,maxref=5000,conformable,auto,unitarity=0.01 '
    '--order 4 --max-distance 1 --comment '
    '--output-transformation {outtransfile} '
    '--output /dev/null'
    )

FRAMETRANSFORMCMD = ('fitrans {frametoshift} -k '
                     '--input-transformation {itransfile} '
                     '--reverse -o {outtransframe}')



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
            goodsrclists.append(srclistpath)


    # we only work on goodframes now
    print('%sZ: selecting an astrometric reference frame...' %
          (datetime.utcnow().isoformat(),))


    median_sval = []
    median_dval = []
    median_background = []
    good_detections = []

    # go through all the frames and find their properties
    for frame, phot, srclist in zip(goodframes, goodphots, goodsrclists):

        if DEBUG:
            print('working on frame %s' % frame)

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

        # find good frames
        if '--binary-output' in header:
            goodind = np.where(photdata['flag'] == 0)
        else:
            goodind = np.where(photdata['flag'] == 'G')
        # number of good detections
        good_detections.append(len(photdata['mag'][goodind]))

        # median background, d, and s
        median_background.append(np.nanmedian(srcdata['background']))
        median_dval.append(np.nanmedian(srcdata['dvalue']))
        median_sval.append(np.nanmedian(srcdata['svalue']))

    #
    # now find the best astrometric reference frame
    #

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

    sd_ind =  np.intersect1d(median_sval_ind,
                             median_dval_ind,
                             assume_unique=True)


    best_frame_ind = np.intersect1d(
        sd_ind,
        np.intersect1d(median_background_ind,
                       good_detections_ind,
                       assume_unique=True),
        assume_unique=True
        )

    sdndet_ind = np.intersect1d(sd_ind,
                                good_detections_ind,
                                assume_unique=True)

    # pick a good astrometric reference frame
    goodframes = np.array(goodframes)


    # if all selectors produced a result, use that one
    if len(best_frame_ind) > 0:

        selectedreference = goodframes[best_frame_ind[0]]

        print('%sZ: selected best astrometric reference frame is %s' %
              (datetime.utcnow().isoformat(), selectedreference))

        return selectedreference

    # otherwise, fall back to to the frames with the best values of S, D and
    # a large number of detections
    elif len(sdndet_ind) > 0:

        selectedreference = goodframes[sdndet_ind[0]]

        print('WRN! %sZ: selected best astrometric reference frame '
              '(using S, D, and ndet only) is %s' %
              (datetime.utcnow().isoformat(), selectedreference))

        return selectedreference


    # otherwise, fall back to to the frames with the best values of S and D
    elif len(sd_ind) > 0:

        selectedreference = goodframes[sd_ind[0]]

        print('WRN! %sZ: selected best astrometric reference frame '
              '(using S and D only) is %s' %
              (datetime.utcnow().isoformat(), selectedreference))

        return selectedreference

    # if that fails, fail back to the best S value frame
    elif len(median_sval_ind) > 0:

        selectedreference = goodframes[median_sval_ind[0]]

        print('WRN! %sZ: selected best astrometric reference frame '
              '(using S only) is %s' %
              (datetime.utcnow().isoformat(), selectedreference))

        return selectedreference

    else:

        print('ERR! %sZ: could not select a good astrometric reference frame!' %
              (datetime.utcnow().isoformat(), ))

        return



def xysdk_coeffs_worker(task):
    '''
    This is a parallel worker to run the xysdk coeff operation.

    '''

    fistar, fistarglob = task

    outxysdk = fistar.replace(fistarglob.split('.')[-1], 'xysdk')

    cmdtorun = XYSDKCMD.format(fistarfile=fistar,
                               xysdkfile=outxysdk)

    returncode = os.system(cmdtorun)

    if returncode == 0:
        print('%sZ: XYSDK coeffs OK: %s -> %s' %
              (datetime.utcnow().isoformat(), fistar, outxysdk))
        return fistar, outxysdk
    else:
        print('ERR! %sZ: XYSDK coeffs failed for %s' %
              (datetime.utcnow().isoformat(), fistar))
        if os.path.exists(outxysdk):
            os.remove(outxysdk)
        return fistar, None


def get_smoothed_xysdk_coeffs(fistardir,
                              fistarglob='*.fistar',
                              nworkers=16,
                              maxworkertasks=1000):
    '''
    This generates smoothed xy and sdk coefficents for use with iphot later
    (these go into the output photometry file or something).

    grtrans ${APPHOT}/$base.${EXT_FISTAR} --col-xy 2,3 --col-fit 6,7,8 \
                  --col-weight 10 --order 4 \
                  --iterations 3 --rejection-level 3 \
                  --comment --output-transformation ${IPHOT}/$base.xysdk

    '''

    fistarlist = glob.glob(os.path.join(os.path.abspath(fistardir), fistarglob))

    print('%sZ: %s files to process in %s' %
          (datetime.utcnow().isoformat(), len(fistarlist), fistardir))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, fistarglob) for x in fistarlist]

    # fire up the pool of workers
    results = pool.map(xysdk_coeffs_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {x:y for (x,y) in results}



def astromref_shift_worker(task):
    '''
    This is a parallel worker for getting the shifts between the astromref frame
    and the frame in the task definition.

    task[0] = target frame fistar
    task[1] = astromref frame fistar
    task[2] = outdir

    '''

    fistartoshift, astromref, outdir = task

    if outdir:
        outfile = os.path.join(
            os.path.abspath(outdir),
            os.path.basename(fistartoshift).replace('fistar','itrans')
            )

    else:
        outfile = fistartoshift.replace('fistar','itrans')

    cmdtorun = FRAMESHIFTCALCCMD.format(astromref=astromref,
                                        fistartoshift=fistartoshift,
                                        outtransfile=outfile)

    returncode = os.system(cmdtorun)

    if returncode == 0:
        print('%sZ: shift transform calc OK: %s -> %s' %
              (datetime.utcnow().isoformat(), fistartoshift, outfile))
        return fistartoshift, outfile
    else:
        print('ERR! %sZ: shift transform calc failed for %s' %
              (datetime.utcnow().isoformat(), fistartoshift))
        if os.path.exists(outfile):
            os.remove(outfile)
        return fistartoshift, None



def get_astromref_shifts(fistardir,
                         astromrefframe,
                         fistarglob='*.fistar',
                         outdir=None,
                         nworkers=16,
                         maxworkertasks=1000):
    '''
    This gets shifts between the astrometric reference frame and all other
    frames.

    '''

    fistarlist = glob.glob(os.path.join(os.path.abspath(fistardir), fistarglob))

    print('%sZ: %s files to process in %s' %
          (datetime.utcnow().isoformat(), len(fistarlist), fistardir))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, astromrefframe, outdir) for x in fistarlist]

    # fire up the pool of workers
    results = pool.map(astromref_shift_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {x:y for (x,y) in results}



def frame_to_astromref_worker(task):
    '''
    This is a parallel worker for the frame shift to astromref frame operation.

    task[0] = FITS frame to shift
    task[1] = directory where transform files are
    task[2] = output directory

    '''

    frametoshift, transdir, outdir = task

    # figure out the transfile path
    if transdir:
        itransfile = os.path.join(
            os.path.join(
                os.path.abspath(transpath),
                os.path.basename(frametoshift).replace('.fits','.itrans')
            )
        )

    else:
        itransfile = frametoshift.replace('.fits','.itrans')

    # make sure the itransfile for this frame exists before we proceed
    if not os.path.exists(itransfile):
        print('ERR! %sZ: frame transform to astromref failed for %s, '
              'no itrans file found' %
              (datetime.utcnow().isoformat(), frametoshift))
        return frametoshift, None

    # figure out the output path
    if outdir:
        outtransframe = os.path.join(
            os.path.abspath(outdir),
            os.path.basename(frametoshift).replace('.fits','-xtrns.fits')
            )

    else:
        outtransframe = frametoshift.replace('.fits','-xtrns.fits')

    cmdtorun = FRAMETRANSFORMCMD.format(itransfile=itransfile,
                                        frametoshift=frametoshift,
                                        outtransframe=outtransframe)

    returncode = os.system(cmdtorun)

    if returncode == 0:
        print('%sZ: transform to astromref OK: %s -> %s' %
              (datetime.utcnow().isoformat(), frametoshift, outtransframe))
        return frametoshift, outtransframe
    else:
        print('ERR! %sZ: transform to astromref failed for %s' %
              (datetime.utcnow().isoformat(), frametoshift))
        if os.path.exists(outtransframe):
            os.remove(outtransframe)
        return frametoshift, None



def transform_frames_to_astromref(fitsdir,
                                  fitsglob='*.fits',
                                  itransdir=None,
                                  outdir=None,
                                  nworkers=16,
                                  maxworkertasks=1000):
    '''
    This shifts all frames to the astrometric reference.

    '''

    fitslist = glob.glob(os.path.join(os.path.abspath(fitsdir), fitsglob))

    print('%sZ: %s files to process in %s' %
          (datetime.utcnow().isoformat(), len(fitslist), fitsdir))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, itransdir, outdir) for x in fitslist]

    # fire up the pool of workers
    results = pool.map(frame_to_astromref_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {x:y for (x,y) in results}


def generate_astromref_registration_info(astromrefsrclist,
                                         outfile,
                                         xycols=(2,3)):
    '''This generates a registration information file using the astrometry
    reference frame. This file is then used by the convolution step somehow to
    figure out the convolution kernel? In any case, it's needed for:

    - generating convolved reference frames to be ultimately stacked into a
      single photometric reference frame

    - do the convolution of the reference frame to each -xtrns target frame when
      doing the image subtraction

    '''

    # get the x and y coordinate columns from the source list (fistar)
    srcxy = np.genfromtxt(astromrefsrclist,
                          usecols=xycols,
                          dtype='f8,f8',
                          names=['x','y'])

    # set up the grid (this weirdness is transcribed directly from Chelsea's
    # regslct.py) TODO: figure out WTF this does

    BX = 30.; BY = 30.
    mx = np.zeros(BX*BY)-1
    my = np.zeros(BX*BY)-1
    ma = np.zeros(BX*BY)
    xsize = 2048.
    ysize = 2048.
    bx = (srcxy['x']*BX/xsize).astype(int)
    by = (srcxy['y']*BY/ysize).astype(int)
    mx[by*bx+bx] = srcxy['x']
    my[by*bx+bx] = srcxy['y']

    outf = open(outfile,'wb')

    for i in xrange(int(BX*BY)):
        outf.write("%8.0f %8.0f %8.0f\n" % (mx[i],my[i],20))

    outf.close()


##################################
## PHOTOMETRIC REFERENCE FRAMES ##
##################################

def select_photref_frames(fitsdir,
                          fitsglob='*-xtrns.fits',
                          photdir=None,
                          photext='.fiphot',
                          srclistdir=None,
                          srclistext='.fistar',
                          minframes=80,
                          forcecollectinfo=False):
    '''This selects a group of photometric reference frames that will later be
    stacked and medianed to form the single photometric reference frame.

    0. this is run on the transformed frames (the ones with -xtrns.fits)

    1. returns a list that is at least minframes long of frames suitable for
    combining into a median reference frame, using the following list of
    selectors (on the original versions (?; probably)):

    - best median scatter of photometry
    - lowest median error in photometry
    - lowest median background measurement
    - low zenith distance
    - high moon distance and lowest moon elevation
    - hour angle with +/- 3 hours
    - large number of stars detected with good flags

    for all selected frames, we will get the median of the background values
    near the center 512x512 pixels. then, we'll enforce that the median of
    background values of each frames be within some delta of the overall
    median. this is basically a slacker way to get rid of cloudy nights.

    2. we convolve all of these to the astrometric reference frame's PSF
    (they're already in the same coordinates as the astromref).

    3. now that all the frames are in the same coordinate system, and have been
    convolved to the same PSF, we can median-stack them (using scipy or
    ficombine)

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

        fitsbase = os.path.splitext(os.path.basename(fits))[0]

        # if the xtrns files are passed in, make sure we look at the
        # right fistar and fiphot files
        if '-xtrns' in fitsbase:
            fitsbase = fitsbase.rstrip('-xtrns')

        photpath = os.path.join(
            photdir,
            fitsbase + photext
            )
        srclistpath = os.path.join(
            srclistdir,
            fitsbase + srclistext
            )

        if os.path.exists(photpath) and os.path.exists(srclistpath):
            goodframes.append(fits)
            goodphots.append(photpath)
            goodsrclists.append(srclistpath)

    # we only work on goodframes now
    print('%sZ: %s good frames found in %s, '
          'now selecting photometric reference frames...' %
          (datetime.utcnow().isoformat(), len(goodframes), fitsdir))

    # things we need to worry about
    # best median scatter of photometry
    # lowest median error in photometry
    # lowest median background measurement
    # low zenith distance
    # high moon and sun distance
    # large number of stars detected

    if (not os.path.exists(os.path.join(fitsdir,
                                         'TM-imagesub-photref.pkl')) or
        forcecollectinfo):

        # from the FITS
        zenithdist, moondist, moonelev, moonphase, hourangle = [], [], [], [], []

        # from the fiphot files
        ngoodobjects, medmagerr, magerrmad, medsrcbg = [], [], [], []


        for frame, phot, srclist in zip(goodframes, goodphots, goodsrclists):

            if DEBUG:
                print('working on frame %s' % frame)

            # 1. get the data from FITS header
            headerdata = imageutils.get_header_keyword_list(
                frame,
                ['Z','MOONDIST','MOONELEV','MOONPH','HA']
                )

            # 2. get the data from the fiphot file

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

                print('WRN! %sZ: %s is a binary phot file, '
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

            # 3. get the data fro mthe fistar file
            srcdata = np.genfromtxt(srclist,
                                    usecols=(3,5,6),
                                    dtype='f8,f8,f8',
                                    names=['background',
                                           'svalue',
                                           'dvalue'])

            # now we have headerdata, photdata, and srcdata, fill in the lists

            # header data
            if 'Z' in headerdata:
                zenithdist.append(headerdata['Z'])
            else:
                zenithdist.append(np.nan)

            if 'MOONDIST' in headerdata:
                moondist.append(headerdata['MOONDIST'])
            else:
                moondist.append(np.nan)

            if 'MOONELEV' in headerdata:
                moonelev.append(headerdata['MOONELEV'])
            else:
                moonelev.append(np.nan)

            if 'MOONPH' in headerdata:
                moonphase.append(headerdata['MOONPH'])
            else:
                moonphase.append(np.nan)

            if 'HA' in headerdata:
                hourangle.append(headerdata['HA'])
            else:
                hourangle.append(np.nan)

            # fiphot data
            if '--binary-output' in header:
                goodind = np.where(photdata['flag'] == 0)
            else:
                goodind = np.where(photdata['flag'] == 'G')

            median_mag = np.median(photdata['mag'][goodind])

            # these are the quantities we're interested in
            ngood = len(goodind[0])
            median_magerr = np.nanmedian(photdata['err'][goodind])
            medabsdev_mag = np.nanmedian(
                np.abs(photdata['mag'][goodind] - median_mag)
                )

            # put these in the lists
            ngoodobjects.append(ngood)
            medmagerr.append(median_magerr)
            magerrmad.append(medabsdev_mag)

            # fistar data
            medsrcbg.append(np.nanmedian(srcdata['background']))

        #
        # done with collecting data, choose the best photometric reference frames
        #

        # convert all lists to np.arrays first
        zenithdist = np.array(zenithdist)
        moondist = np.array(moondist)
        moonelev = np.array(moonelev)
        moonphase = np.array(moonphase)
        hourangle = np.array(hourangle)

        ngoodobjects = np.array(ngoodobjects)
        medmagerr = np.array(medmagerr)
        magerrmad = np.array(magerrmad)

        medsrcbg = np.array(medsrcbg)

        goodframes = np.array(goodframes)

        infodict = {
            'frames':goodframes,
            'zenithdist':zenithdist,
            'moondist':moondist,
            'moonelev':moonelev,
            'moonphase':moonphase,
            'hourangle':hourangle,
            'ngoodobjs':ngoodobjects,
            'medmagerr':medmagerr,
            'magerrmad':magerrmad,
            'medsrcbkg':medsrcbg
        }

        # write this info dict to a file so we can quickly load it later
        outpf = open(os.path.join(fitsdir, 'TM-imagesub-photref.pkl'), 'wb')
        pickle.dump(infodict, outpf, pickle.HIGHEST_PROTOCOL)
        outpf.close()

    # if the imagesub photref info file exists already, load it up
    else:

        print('%sZ: wrote the photref select info to %s' %
              (datetime.utcnow().isoformat(),
               os.path.join(fitsdir, 'TM-imagesub-photref.pkl')))

        inpf = open(os.path.join(fitsdir, 'TM-imagesub-photref.pkl'), 'rb')
        infodict = pickle.load(inpf)
        inpf.close()

    #
    # now do the filtering
    #

    # filter on hour angle
    haind = np.where(np.fabs(infodict['hourangle']) < 3.0)
    haind = haind[0]

    # get dark nights
    moonind = np.where((np.fabs(infodict['moonphase']) < 25.0) |
                       (moonelev < -15.0))
    moonind = moonind[0]

    selectind = haind & moonind

    print selectind, len(selectind)
    return infodict





def convolve_photref_frames(photreflist,
                            targetframe,
                            outdir=None):
    '''This convolves all photref frames to the targetframe. See getref() in
    run_ficonv.py.

    '''



def stack_photref_frames(framelist):
    '''
    This stacks the photometric reference frames in median mode.

    - first, transform all of them to the astrometric reference frame

    - then use ficombine in median mode to stack them

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

    fiphot --input $Refframe --input-list $Inlist --col-id 1 --col-xy 5,6
           --gain 2.725 --mag-flux $zp,30
           --apertures '2.5:7.0:6.0,2.95:7.0:6.0,3.35:7.0:6.0,3.95:7.0:6.0'
           --sky-fit 'mode,sigma=3,iterations=2' --disjoint-radius 2
           --serial $seriel --format 'IXY-----,sMm' --nan-string 'NaN'
           --aperture-mask-ignore 'saturated' --comment '--comment'
           --single-background 3 -op $output -k

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

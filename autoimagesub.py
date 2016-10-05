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
import sqlite3
import time
from hashlib import md5, sha256
import gzip
from traceback import format_exc

import numpy as np
import psycopg2 as pg

import aperturephot as ap
import imagesubphot as ism
from imageutils import get_header_keyword, get_header_keyword_list, \
    fits_to_full_jpeg, check_frame_warping

from numpy import nan

# get fiphot binary reader
try:
    from HATpipepy.Common.BinPhot import read_fiphot
    HAVEBINPHOT = True
except:
    print("can't import binary fiphot reading functions from "
          "HATpipe, binary fiphot files will be unreadable!")
    HAVEBINPHOT = False


############
## CONFIG ##
############

DEBUG = False

# used to get the station ID, frame number, and CCD number from a FITS filename
FRAMEREGEX = re.compile(r'(\d{1})\-(\d{6}\w{0,1})_(\d{1})')

# used to get the station ID, frame number, subframe, and CCD number
FRAMESUBREGEX = re.compile(r'(\d{1})\-(\d{6})(\w{0,1})_(\d{1})')

# this defines the field string and CCDs
FIELD_REGEX = re.compile('^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$')
FIELD_CCDS = [5,6,7,8]

# defines where the reference frames go
REFBASEDIR = '/P/HP0/BASE/reference-frames'
REFINFO = os.path.join(REFBASEDIR,'TM-refinfo.sqlite')

# define where the frameinfo cache is
FRAMEINFOCACHEDIR = '/P/HP0/BASE/frameinfo-cache'

# these define the field catalog location and properties
FIELDCAT_DIR = '/P/HP0/BASE/field-catalogs'

# these define the postgres database credentials
PGPASSFILE = '/home/hatuser/.pgpass'
PGUSER = 'hpx'
PGDATABASE = 'hpx'
PGHOST = 'localhost'

# these define the light curve directory
LCBASEDIR = '/P/LC'

# this is to recognize a HATID
HATIDREGEX = re.compile(r'^HAT\-\d{3}\-\d{7}$')


with open(PGPASSFILE) as infd:
    pgpass_contents = infd.readlines()
    pgpass_contents = [x.split(':') for x in pgpass_contents]
    PGPASSWORD = [x[-1] for x in pgpass_contents
                  if (x[0] == PGHOST and x[2] == PGDATABASE and x[3] == PGUSER)]
    PGPASSWORD = PGPASSWORD[0].strip('\n')



###############
## UTILITIES ##
###############

def smartcast(castee, caster, subval=None):
    '''
    This just tries to apply the caster function to castee.

    Returns None on failure.

    '''

    try:
        return caster(castee)
    except Exception as e:
        if caster is float or caster is int:
            return nan
        elif caster is str:
            return ''
        else:
            return subval



def fits_fieldprojectidccd_worker(frame):
    '''
    This is a worker for the two functions below.

    '''

    try:

        # first, figure out the input frame's projid, field, and ccd
        frameelems = get_header_keyword_list(frame,
                                             ['object',
                                              'projid'])
        felems = FRAMEREGEX.findall(
            os.path.basename(frame)
        )
        field, ccd, projectid = (frameelems['object'],
                                 felems[0][2],
                                 frameelems['projid'])

        return frame, (field, projectid, int(ccd))

    except Exception as e:

        print('ERR! %sZ: could get info from frame %s, error was: %s' %
              (datetime.utcnow().isoformat(), frame, e))
        return frame, None



def find_original_fits_in_database(field, projectid, ccd,
                                   enforceok=True,
                                   database=None):
    '''
    This finds the FITS matching the field-projectid-ccd combo in the DB.

    '''




def find_original_fits_fieldprojectidccd(dirlist,
                                         field,
                                         projectid,
                                         ccd,
                                         fglob='?-???????_?.fits',
                                         nworkers=8,
                                         maxworkertasks=1000):
    '''This searches in dirlist for all original FITS files matching the specified
    projectid, field, and ccd combination.

    Returns a flat list of matching FITS, and list of all fits + their info.

    '''

    # first, go through the directories and get all the original FITS files
    print('%sZ: finding frames matching %s...' %
          (datetime.utcnow().isoformat(), fglob))
    fitslist = []
    for fdir in dirlist:
        fitslist.extend(glob.glob(os.path.join(fdir, fglob)))
    fitslist = sorted(fitslist)

    # next, run through all these files and get the info needed
    print('%sZ: %s frames found, getting info...' %
          (datetime.utcnow().isoformat(), len(fitslist)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = fitslist

    # fire up the pool of workers
    results = pool.map(fits_fieldprojectidccd_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    # now filter the results based on the requested field, projectid, and ccd
    matchingframes = []

    for elem in results:

        if (elem[1] and
            elem[1][0] == field and
            elem[1][1] == projectid and
            elem[1][2] == ccd):
            matchingframes.append(elem[0])

    print('%sZ: %s frames with field = %s, projectid = %s, and ccd = %s' %
          (datetime.utcnow().isoformat(),
           len(matchingframes),
           field, projectid, ccd))

    return matchingframes, results



def find_arefshifted_fits_fieldprojectidccd(dirlist,
                                            field,
                                            projectid,
                                            ccd,
                                            fglob='?-???????_?-xtrns.fits',
                                            nworkers=8,
                                            maxworkertasks=1000):
    '''This searches for all astromref-shifted FITS files matching the
    specified projectid, field, and ccd combination.

    Returns a flat list of matching FITS, and list of all fits + their info.

    '''

    return find_original_fits_fieldprojectidccd(dirlist,
                                                field,
                                                projectid,
                                                ccd,
                                                fglob=fglob,
                                                nworkers=nworkers,
                                                maxworkertasks=maxworkertasks)



def find_subtracted_fits_fieldprojectidccd(
        dirlist,
        field,
        projectid,
        ccd,
        subtracttype,
        photreftype,
        nworkers=8,
        maxworkertasks=1000
):
    '''This searches for all subtracted FITS files matching the specified
    projectid, field, and ccd combination.

    Returns a flat list of matching FITS, and list of all fits + their info.

    '''

    fglob= '%s-%s-?-???????_?-xtrns.fits' % (subtracttype, photreftype)

    return find_original_fits_fieldprojectidccd(dirlist,
                                                field,
                                                projectid,
                                                ccd,
                                                fglob=fglob,
                                                nworkers=nworkers,
                                                maxworkertasks=maxworkertasks)

########################
## GETTING FRAME INFO ##
########################

def get_frame_info(frame):
    '''
    This gets the needed info from a frame for selecting photref candidates.

    '''

    try:

        # find the frame's fistar and fiphot files
        fitsdir = os.path.dirname(frame)
        fitsbase = os.path.splitext(os.path.basename(frame))[0]

        # if the xtrns files are passed in, make sure we look at the
        # right fistar and fiphot files
        if '-xtrns' in fitsbase:
            fitsbase = fitsbase.rstrip('-xtrns')

        photpath = os.path.join(
            fitsdir,
            fitsbase + '.fiphot'
            )
        srclistpath = os.path.join(
            fitsdir,
            fitsbase + '.fistar'
            )

        if not (os.path.exists(frame) and
                os.path.exists(photpath) and
                os.path.exists(srclistpath)):

            print('ERR! %sZ: frame %s has missing fiphot/fistar, '
                  'so no phot info...' %
                  (datetime.utcnow().isoformat(), frame))
            return frame, None



        # 1. get the data from FITS header
        headerdata = get_header_keyword_list(
            frame,
            ['Z','MOONDIST','MOONELEV','MOONPH','HA']
        )

        # 2. get the data from the fiphot file

        # decide if the phot file is binary or not. read the first 600
        # bytes and look for the '--binary-output' text
        with open(photpath,'rb') as photf:
            header = photf.read(1000)

        if '--binary-output' in header and HAVEBINPHOT:

            photdata_f = read_fiphot(photpath)
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
                  (datetime.utcnow().isoformat(), photpath))
            return frame, None

        else:

            # read in the phot file
            photdata = np.genfromtxt(
                photpath,
                usecols=(12,13,14),
                dtype='f8,f8,S5',
                names=['mag','err','flag']
                )

        # 3. get the data from the fistar file
        srcdata = np.genfromtxt(srclistpath,
                                usecols=(3,5,6),
                                dtype='f8,f8,f8',
                                names=['background',
                                       'svalue',
                                       'dvalue'])

        # process fiphot data
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


        # now consolidate all the data
        frameinfo = {
            'zenithdist':(headerdata['Z']
                          if 'Z' in headerdata else np.nan),
            'moondist':(headerdata['MOONDIST']
                        if 'MOONDIST' in headerdata else np.nan),
            'moonelev':(headerdata['MOONELEV']
                        if 'MOONELEV' in headerdata else np.nan),
            'moonphase':(headerdata['MOONPH']
                         if 'MOONPH' in headerdata else np.nan),
            'hourangle':(headerdata['HA']
                         if 'HA' in headerdata else np.nan),
            'ngoodobjects':ngood,
            'medmagerr':median_magerr,
            'magerrmad':medabsdev_mag,
            'medsrcbgv':np.nanmedian(srcdata['background']),
            'stdsrcbgv':np.nanstd(srcdata['background']),
            'medsval':np.nanmedian(srcdata['svalue']),
            'meddval':np.nanmedian(srcdata['dvalue']),
        }

        return frame, frameinfo

    except Exception as e:

        print('ERR! %sZ: could not get info from frame %s, error was: %s' %
              (datetime.utcnow().isoformat(), frame, e))
        return frame, None



def fitslist_frameinfo(fitslist,
                       forcecollectinfo=False,
                       nworkers=8,
                       maxworkertasks=1000):
    '''
    This runs a parallel get_frame_info job.

    '''

    # check if we have it in the frameinfo cache, if so, use it. if not, redo
    # the info collection, and then write it back to the cache.
    cachefile = os.path.join(FRAMEINFOCACHEDIR,
                             ('TM-frameinfo-%s.pkl.gz' %
                              md5(repr(fitslist)).hexdigest()))

    if os.path.exists(cachefile) and not forcecollectinfo:

        with gzip.open(cachefile,'rb') as infd:
            frameinfo = pickle.load(infd)

        print('%sZ: frameinfo found in cache file: %s' %
              (datetime.utcnow().isoformat(), cachefile))

        return frameinfo

    # if the cache doesn't exist, we'll run the frameinfo procedure and write
    # the results back to the cache
    else:

        print('%sZ: getting frameinfo for %s frames...' %
              (datetime.utcnow().isoformat(), len(fitslist)))

        pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

        tasks = fitslist

        # fire up the pool of workers
        results = pool.map(get_frame_info, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()

        # now turn everything into ndarrays
        frames = np.array([x[0] for x in results])
        zenithdist = np.array([(x[1]['zenithdist']
                                if x[1] else np.nan) for x in results])
        moondist = np.array([(x[1]['moondist']
                                if x[1] else np.nan) for x in results])
        moonelev = np.array([(x[1]['moonelev']
                                if x[1] else np.nan) for x in results])
        moonphase = np.array([(x[1]['moonphase']
                                if x[1] else np.nan) for x in results])
        hourangle = np.array([(x[1]['hourangle']
                                if x[1] else np.nan) for x in results])
        ngoodobjects = np.array([(x[1]['ngoodobjects']
                                if x[1] else np.nan) for x in results])
        medmagerr = np.array([(x[1]['medmagerr']
                                if x[1] else np.nan) for x in results])
        magerrmad = np.array([(x[1]['magerrmad']
                                if x[1] else np.nan) for x in results])
        medsrcbgv = np.array([(x[1]['medsrcbgv']
                                if x[1] else np.nan) for x in results])
        stdsrcbgv = np.array([(x[1]['stdsrcbgv']
                                if x[1] else np.nan) for x in results])
        medsval = np.array([(x[1]['medsval']
                                if x[1] else np.nan) for x in results])
        meddval = np.array([(x[1]['meddval']
                                if x[1] else np.nan) for x in results])


        outdict = {'frames':frames,
                   'zenithdist':zenithdist,
                   'moondist':moondist,
                   'moonelev':moonelev,
                   'moonphase':moonphase,
                   'hourangle':hourangle,
                   'ngoodobjects':ngoodobjects,
                   'medmagerr':medmagerr,
                   'magerrmad':magerrmad,
                   'medsrcbgv':medsrcbgv,
                   'stdsrcbgv':stdsrcbgv,
                   'medsval':medsval,
                   'meddval':meddval}

        with gzip.open(cachefile,'wb') as outfd:
            pickle.dump(outdict, outfd, pickle.HIGHEST_PROTOCOL)

        print('%sZ: wrote frameinfo to cache file: %s' %
              (datetime.utcnow().isoformat(), cachefile))

        return outdict



######################################
## PUTTING FRAMES INTO THE DATABASE ##
######################################

def calibrated_frame_to_database(fitsfile,
                                 network='HP',
                                 overwrite=False,
                                 badframetag='badframes',
                                 database=None):
    '''This puts an original fully calibrated FITS into the database.

    Requires that the frame have gone through the full
    cron_transfer_calibrate.autocal_fitsdir function process, so that it has an
    extracted source list (fistar), astrometry attempted on it (wcs), and basic
    aperture photometry carried out (fiphot).

    Searches the same directory as the FITS for fistar, fiphot, and wcs
    files. Should be used with a "find {fitsdir} -type f -name '?-*_?.fits'"
    command to generate the list of the FITS files, so that it can work on bad
    frames that have been moved to the {badframedir} subdirectory of {fitsdir}
    and import these into the database as well. Bad frames will be tagged with
    frameisok = false in the database. If a frame doesn't have an accompanying
    wcs or fiphot, it will be tagged with frameisok = false in the database as
    well.

    '''
    # open a database connection
    if database:
        cursor = database.cursor()
        closedb = False
    else:
        database = pg.connect(user=PGUSER,
                              password=PGPASSWORD,
                              database=PGDATABASE,
                              host=PGHOST)
        cursor = database.cursor()
        closedb = True

    # start work here
    try:

        # get the frame filename elements
        felems = FRAMESUBREGEX.findall(
            os.path.basename(fitsfile)
        )
        if not felems:
            message = ('failed to get basic frameinfo from %s, skipping...'
                       % fitsfile)
            print('ERR! %sZ: %s' %
                   (datetime.utcnow().isoformat(), message) )
            return fitsfile, False


        # get a big list of header keywords
        headerdata = get_header_keyword_list(
            fitsfile,
            ['PROJID','OBJECT','JD','EXPTIME',
             'STVER','MTID','MTVER',
             'CMID','CMVER','TELID','TELVER','FOV',
             'BIASVER','DARKVER','FLATVER','MNTSTATE','PTVER',
             'FILID','IMAGETYP','RA','DEC',
             'MOONDIST','MOONELEV','MOONPH',
             'HA','Z','MGENSTAT','WIND','HUMIDITY','SKYTDIFF','AMBTEMP',
             'DEWPT','FOCUS', 'COMMENT'],
        )

        # sort out the commentary card
        headercomments = repr(headerdata['COMMENT'])
        headercomments = headercomments.lower()

        # get the frame's projectid, stationid, observed field, rjd, centerra,
        # centerrdec, widthdeg, frameisok values
        if 'PROJID' in headerdata and headerdata['PROJID']:
            projectid = headerdata['PROJID']
        else:
            projectid = None
        if 'OBJECT' in headerdata and headerdata['OBJECT']:
            obsfield = headerdata['OBJECT']
        else:
            obsfield = None
        if 'JD' in headerdata and headerdata['JD']:
            framerjd = headerdata['JD']
        else:
            framerjd = None

        # convert the RA [h] to RA (deg)
        if 'RA' in headerdata and headerdata['RA'] is not None:
            centerra = (headerdata['RA']/24.0)*360.0
        else:
            centerra = None
        if 'DEC' in headerdata and headerdata['DEC'] is not None:
            centerdec = headerdata['DEC']
        else:
            centerdec = None

        if 'FOV' in headerdata and headerdata['FOV']:
            fovdeg = headerdata['FOV']
        else:
            fovdeg = None

        # get the frame's framenum, ccdnum, subframe, frametype
        stationid, cfn, cfs, ccd = felems[0]
        if 'IMAGETYP' in headerdata and headerdata['IMAGETYP']:
            frt = headerdata['IMAGETYP']
        else:
            frt = None

        # find the frame's fistar, wcs, and fiphot
        prospective_fistar = fitsfile.replace('.fits','.fistar')
        prospective_wcs = fitsfile.replace('.fits','.wcs')
        prospective_fiphot = fitsfile.replace('.fits','.fiphot')

        # check if this is a bad frame and set stuff accordingly
        if (os.path.exists(prospective_fistar) and
            os.path.exists(prospective_wcs) and
            os.path.exists(prospective_fiphot) and
            (badframetag not in os.path.abspath(fitsfile))):

            frameisok = True
            fits = os.path.abspath(fitsfile)
            fistar = os.path.abspath(prospective_fistar)
            wcs = os.path.abspath(prospective_wcs)
            fiphot = os.path.abspath(prospective_fiphot)

        else:

            frameisok = False

            fits = os.path.abspath(fitsfile)

            if os.path.exists(prospective_fistar):
                fistar = os.path.abspath(prospective_fistar)
            else:
                fistar = None
            if os.path.exists(prospective_wcs):
                wcs = os.path.abspath(prospective_wcs)
            else:
                wcs = None
            if os.path.exists(prospective_fiphot):
                fiphot = os.path.abspath(prospective_fiphot)
            else:
                fiphot = None


        # get the frame's filter ID and version
        if 'FILID' in headerdata and headerdata['FILID'] is not None:
            flt = headerdata['FILID']
        else:
            flt = None
        if 'filters.flver' in headercomments:
            flvind = headercomments.find('filters.flver')
            flv = headercomments[flvind:flvind+16]
            flv = int(flv.split('=')[-1].rstrip("'\""))
        else:
            flv = 0

        # get the frame's camera config
        if 'CMID' in headerdata and headerdata['CMID'] is not None:
            cid = headerdata['CMID']
        else:
            cid = None
        if 'CMVER' in headerdata and headerdata['CMVER'] is not None:
            cvn = headerdata['CMVER']
        else:
            cvn = 0
        if 'BIASVER' in headerdata and headerdata['BIASVER'] is not None:
            cbv = headerdata['BIASVER']
        else:
            cbv = 0
        if 'DARKVER' in headerdata and headerdata['DARKVER'] is not None:
            cdv = headerdata['DARKVER']
        else:
            cdv = 0
        if 'FLATVER' in headerdata and headerdata['FLATVER'] is not None:
            cfv = headerdata['FLATVER']
        else:
            cfv = 0
        if 'EXPTIME' in headerdata and headerdata['EXPTIME'] is not None:
            exp = headerdata['EXPTIME']
        else:
            exp = None

        # get the frame's telescope config
        if 'TELID' in headerdata and headerdata['TELID'] is not None:
            tid = headerdata['TELID']
        else:
            tid = None
        if 'TELVER' in headerdata and headerdata['TELVER'] is not None:
            tvn = headerdata['TELVER']
        else:
            tvn = 0
        if 'FOCUS' in headerdata and headerdata['FOCUS'] is not None:
            tfs = headerdata['FOCUS']
        else:
            tfs = None
        if 'MNTSTATE' in headerdata and headerdata['MNTSTATE'] is not None:
            tms = headerdata['MNTSTATE']
        else:
            tms = None
        if 'MTID' in headerdata and headerdata['MTID'] is not None:
            tmi = headerdata['MTID']
        else:
            tmi = None
        if 'MTVER' in headerdata and headerdata['MTVER'] is not None:
            tmv = headerdata['MTVER']
        else:
            tmv = 0
        if 'MGENSTAT' in headerdata and headerdata['MGENSTAT'] is not None:
            tgs = headerdata['MGENSTAT']
        else:
            tgs = None

        # get the frame's photometry info (useful for selecting refs)
        photfits, photinfo = get_frame_info(fitsfile)

        # make sure the frame actually has photinfo before inserting
        if photinfo:
            ngo = photinfo['ngoodobjects']
            mme = photinfo['medmagerr']
            mem = photinfo['magerrmad']
            mbg = photinfo['medsrcbgv']
            sbg = photinfo['stdsrcbgv']
            mfs = photinfo['medsval']
            mfd = photinfo['meddval']
        else:
            ngo = None
            mme = None
            mem = None
            mbg = None
            sbg = None
            mfs = None
            mfd = None
            frameisok = False


        # get the frame's environmental info
        if 'MOONPH' in headerdata and headerdata['MOONPH'] is not None:
            mph = headerdata['MOONPH']
        else:
            mph = None
        if 'MOONPH' in headerdata and headerdata['MOONPH'] is not None:
            mph = headerdata['MOONPH']
        else:
            mph = None
        if 'MOONDIST' in headerdata and headerdata['MOONDIST'] is not None:
            mds = headerdata['MOONDIST']
        else:
            mds = None
        if 'MOONELEV' in headerdata and headerdata['MOONELEV'] is not None:
            mel = headerdata['MOONELEV']
        else:
            mel = None
        if 'HA' in headerdata and headerdata['HA'] is not None:
            iha = headerdata['HA']
        else:
            iha = None
        if 'Z' in headerdata and headerdata['Z'] is not None:
            izd = headerdata['Z']
        else:
            izd = None
        if 'WIND' in headerdata and headerdata['WIND'] is not None:
            wis = headerdata['WIND']
        else:
            wis = None
        if 'HUMIDITY' in headerdata and headerdata['HUMIDITY'] is not None:
            hum = headerdata['HUMIDITY']
        else:
            hum = None
        if 'SKYTDIFF' in headerdata and headerdata['SKYTDIFF'] is not None:
            skt = headerdata['SKYTDIFF']
        else:
            skt = None
        if 'AMBTEMP' in headerdata and headerdata['AMBTEMP'] is not None:
            amt = headerdata['AMBTEMP']
        else:
            amt = None
        if 'DEWPT' in headerdata and headerdata['DEWPT'] is not None:
            dew = headerdata['DEWPT']
        else:
            dew = None

        # put together the query and execute it, inserting the object into the
        # database and overwriting if told to do so
        if overwrite:

            query = ("insert into calibratedframes ("
                     "network, projectid, stationid, obsfield, "
                     "framerjd, centerra, centerdec, fovdeg, frameisok, "
                     "fits, fistar, fiphot, wcs, "
                     "cfn, cfs, ccd, frt, "
                     "flt, flv, "
                     "cid, cvn, cbv, cdv, cfv, exp, "
                     "tid, tvn, tfs, tms, tmi, tmv, tgs, "
                     "ngo, mme, mem, mbg, sbg, mfs, mfd, "
                     "mph, mds, mel, iha, izd, wis, hum, skt, amt, dew"
                     ") values ("
                     "%s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, "
                     "%s, %s, %s, %s, "
                     "%s, %s, "
                     "%s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                     ") "
                     "on conflict "
                     "(network, projectid, stationid, obsfield, "
                     "cfn, cfs, ccd, frt) "
                     "do update set "
                     "framerjd = %s, centerra = %s, centerdec = %s, "
                     "fovdeg = %s, frameisok = %s, fits = %s, fistar = %s, "
                     "fiphot = %s, wcs = %s, flt = %s, flv = %s, "
                     "cid = %s, cvn = %s, cbv = %s, cdv = %s, cfv = %s, "
                     "exp = %s, tid = %s, tvn = %s, tfs = %s, tms = %s, "
                     "tmi = %s, tmv = %s, tgs = %s, ngo = %s, mme = %s, "
                     "mem = %s, mbg = %s, sbg = %s, mfs = %s, mfd = %s, "
                     "mph = %s, mds = %s, mel = %s, iha = %s, izd = %s, "
                     "wis = %s, hum = %s, skt = %s, amt = %s, dew = %s, "
                     "entryts = current_timestamp")
            params = (network, projectid, stationid, obsfield,
                      framerjd, centerra, centerdec, fovdeg, frameisok,
                      fits, fistar, fiphot, wcs,
                      cfn, cfs, ccd, frt,
                      flt, flv,
                      cid, cvn, cbv, cdv, cfv, exp,
                      tid, tvn, tfs, tms, tmi, tmv, tgs,
                      ngo, mme, mem, mbg, sbg, mfs, mfd,
                      mph, mds, mel, iha, izd, wis, hum, skt, amt, dew,
                      framerjd, centerra, centerdec, fovdeg, frameisok,
                      fits, fistar, fiphot, wcs, flt, flv,
                      cid, cvn, cbv, cdv, cfv,
                      exp, tid, tvn, tfs, tms,
                      tmi, tmv, tgs, ngo, mme,
                      mem, mbg, sbg, mfs, mfd,
                      mph, mds, mel, iha, izd,
                      wis, hum, skt, amt, dew)

        else:

            query = ("insert into calibratedframes ("
                     "network, projectid, stationid, obsfield, "
                     "framerjd, centerra, centerdec, fovdeg, frameisok, "
                     "fits, fistar, fiphot, wcs, "
                     "cfn, cfs, ccd, frt, "
                     "flt, flv, "
                     "cid, cvn, cbv, cdv, cfv, exp, "
                     "tid, tvn, tfs, tms, tmi, tmv, tgs, "
                     "ngo, mme, mem, mbg, sbg, mfs, mfd, "
                     "mph, mds, mel, iha, izd, wis, hum, skt, amt, dew"
                     ") values ("
                     "%s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, "
                     "%s, %s, %s, %s, "
                     "%s, %s, "
                     "%s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                     ") ")
            params = (network, projectid, stationid, obsfield,
                      framerjd, centerra, centerdec, fovdeg, frameisok,
                      fits, fistar, fiphot, wcs,
                      cfn, cfs, ccd, frt,
                      flt, flv,
                      cid, cvn, cbv, cdv, cfv, exp,
                      tid, tvn, tfs, tms, tmi, tmv, tgs,
                      ngo, mme, mem, mbg, sbg, mfs, mfd,
                      mph, mds, mel, iha, izd, wis, hum, skt, amt, dew)

        # execute the query and commit
        cursor.execute(query, params)
        database.commit()

        message = 'inserted %s into DB OK' % fits
        print('%sZ: %s' %
              (datetime.utcnow().isoformat(), message) )
        returnval = (fits, True)

    # catch the overwrite = False scenario
    except pg.IntegrityError as e:

        database.rollback()

        message = ('failed to insert %s '
                   'into DB because it exists already '
                   'and overwrite = False'
                   % fits)
        print('EXC! %sZ: %s\n%s' %
               (datetime.utcnow().isoformat(), message, format_exc()) )
        returnval = (fits, False)


    # if everything goes wrong, exit cleanly
    except Exception as e:

        database.rollback()

        message = 'failed to insert %s into DB' % fits
        print('EXC! %sZ: %s\nexception was: %s' %
               (datetime.utcnow().isoformat(),
                message, format_exc()) )
        returnval = (fits, False)

        # TEMPORARY
        # raise


    finally:

        cursor.close()
        if closedb:
            database.close()

    return returnval



def calframe_to_db_worker(task):
    '''
    This wraps calibrated_frames_to_databse for the parallel driver below.

    '''

    fitsfile, kwargs = task
    return calibrated_frame_to_database(fitsfile, **kwargs)



def parallel_calibrated_frames_to_database(fitsbasedir,
                                           fitsglob='1-???????_?.fits',
                                           network='HP',
                                           overwrite=False,
                                           badframetag='badframes',
                                           nworkers=16,
                                           maxworkertasks=1000):
    '''
    This runs a DB ingest on all FITS located in fitsbasedir and subdirs.

    Runs a 'find' subprocess to find all the FITS to process.

    '''

    print('%sZ: finding all FITS frames matching %s starting in %s' %
          (datetime.utcnow().isoformat(), fitsglob, fitsbasedir))

    # find all the FITS files
    try:

        findcmd = "find {fitsbasedir} -type f -name '{fitsglob}' -print"
        findcmd = findcmd.format(fitsbasedir=fitsbasedir,
                                 fitsglob=fitsglob)
        fitslist = subprocess.check_output(findcmd,shell=True)
        fitslist = fitslist.split('\n')
        fitslist = sorted(fitslist[:-1])

        # generate the task list
        tasks = [(x, {'network':network,
                      'overwrite':overwrite,
                      'badframetag':badframetag}) for x in fitslist]

        print('%sZ: %s files to process' %
              (datetime.utcnow().isoformat(), len(tasks)))

        if len(tasks) > 0:

            pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

            # fire up the pool of workers
            results = pool.map(calframe_to_db_worker, tasks)

            # wait for the processes to complete work
            pool.close()
            pool.join()

            return {x:y for (x,y) in results}

        else:

            print('ERR! %sZ: none of the files specified exist, bailing out...' %
                  (datetime.utcnow().isoformat(),))
            return

    except subprocess.CalledProcessError:

        print('ERR! %sZ: no files in directory %s, bailing out...' %
              (datetime.utcnow().isoformat(), fitsbasedir))
        return




##################################
## ASTROMETRIC REFERENCE FRAMES ##
##################################

def dbgen_astromref_projectidfieldccd(projectid,
                                      field,
                                      ccd,
                                      makeactive=True,
                                      overwrite=False,
                                      refdir=REFBASEDIR,
                                      database=None):
    '''
    This gets all the frame info from the DB and finds a good astromref.

    '''
    # open a database connection
    if database:
        cursor = database.cursor()
        closedb = False
    else:
        database = pg.connect(user=PGUSER,
                              password=PGPASSWORD,
                              database=PGDATABASE,
                              host=PGHOST)
        database.autocommit = True
        cursor = database.cursor()
        closedb = True

    try:

        # get all the frame info for the requested projectid-field-ccd
        query = ("select framekey, fits, fistar, fiphot, mfs, mfd, mbg, ngo "
                 "from calibratedframes where "
                 "(projectid = %s) and (obsfield = %s) and (ccd = %s) and "
                 "frameisok = true")
        params = (str(projectid), field, ccd)

        cursor.execute(query, params)
        rows = cursor.fetchall()

        # if we're successful
        if rows and len(rows) > 0:

            # get the frame info

            framekey = np.array([x[0] for x in rows])
            fits = np.array([x[1] for x in rows])
            fistar = np.array([x[2] for x in rows])
            fiphot = np.array([x[3] for x in rows])

            mfs = np.array([x[4] for x in rows])
            mfd = np.array([x[5] for x in rows])
            mbg = np.array([x[6] for x in rows])
            ngo = np.array([x[7] for x in rows])

            #
            # now, find the best astrometric reference frame
            #

            # get the best S --> largest S at the top
            median_sval_ind = np.argsort(mfs)[::-1]

            # here, we want the values closest to zero to be at the top
            median_dval_ind = np.argsort(np.fabs(mfd))

            # want the smallest background
            median_background_ind = np.argsort(mbg)

            # and the most number of detections
            good_detections_ind = np.argsort(ngo)[::-1]

            # get the top 500 of each index
            median_sval_ind = median_sval_ind[:500]
            median_dval_ind = median_dval_ind[:500]
            median_background_ind = median_background_ind[:500]
            good_detections_ind = good_detections_ind[:500]

            # now intersect all of these arrays to find the best candidates for
            # the astrometric reference frame

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

            # if all selectors produced a result, use that one
            if len(best_frame_ind) > 0:

                selectedreference = fits[best_frame_ind[0]]

                print('%sZ: selected best astrometric reference frame is %s' %
                      (datetime.utcnow().isoformat(), selectedreference))

                framejpg = fits_to_full_jpeg(
                    selectedreference,
                    out_fname=os.path.join(
                            os.path.dirname(selectedreference),
                        ('JPEG-ASTROMREF-%s.jpg' %
                         os.path.basename(selectedreference).strip('.fits.fz'))
                    )
                )

                arefinfo = {
                    'fits':selectedreference,
                    'jpg':framejpg,
                    'fistar':fistar[best_frame_ind[0]],
                    'fiphot':fiphot[best_frame_ind[0]],
                    'framekey':framekey[best_frame_ind[0]],
                    'sval':mfs[best_frame_ind[0]],
                    'dval':mfd[best_frame_ind[0]],
                    'bgv':mbg[best_frame_ind[0]],
                    'ndet':ngo[best_frame_ind[0]],
                    'comment':'astromref chosen using sval, dval, bgval, ndet'
                }

            # otherwise, fall back to frames with the best values of S, D, and a
            # large number of detections
            elif len(sdndet_ind) > 0:

                selectedreference = fits[sdndet_ind[0]]

                print('WRN! %sZ: selected best astrometric reference frame '
                      '(using S, D, and ndet only) is %s' %
                      (datetime.utcnow().isoformat(), selectedreference))


                framejpg = fits_to_full_jpeg(
                    selectedreference,
                    out_fname=os.path.join(
                            os.path.dirname(selectedreference),
                        ('JPEG-ASTROMREF-%s.jpg' %
                         os.path.basename(selectedreference).strip('.fits.fz'))
                    )
                )

                arefinfo = {
                    'fits':selectedreference,
                    'jpg':framejpg,
                    'fistar':fistar[sdndet_ind[0]],
                    'fiphot':fiphot[sdndet_ind[0]],
                    'framekey':framekey[sdndet_ind[0]],
                    'sval':mfs[sdndet_ind[0]],
                    'dval':mfd[sdndet_ind[0]],
                    'bgv':mbg[sdndet_ind[0]],
                    'ndet':ngo[sdndet_ind[0]],
                    'comment':'astromref chosen using sval, dval, ndet'
                }

            # otherwise, fall back to the frames with best values of S and D
            elif len(sd_ind) > 0:

                selectedreference = fits[sd_ind[0]]

                print('WRN! %sZ: selected best astrometric reference frame '
                      '(using S and D only) is %s' %
                      (datetime.utcnow().isoformat(), selectedreference))

                framejpg = fits_to_full_jpeg(
                    selectedreference,
                    out_fname=os.path.join(
                            os.path.dirname(selectedreference),
                        ('JPEG-ASTROMREF-%s.jpg' %
                         os.path.basename(selectedreference).strip('.fits.fz'))
                    )
                )

                arefinfo = {
                    'fits':selectedreference,
                    'jpg':framejpg,
                    'fistar':fistar[sd_ind[0]],
                    'fiphot':fiphot[sd_ind[0]],
                    'framekey':framekey[sd_ind[0]],
                    'sval':mfs[sd_ind[0]],
                    'dval':mfd[sd_ind[0]],
                    'bgv':mbg[sd_ind[0]],
                    'ndet':ngo[sd_ind[0]],
                    'comment':'astromref chosen using sval, dval'
                }

            # if that also fails, fall back to the best S value frame
            elif len(median_sval_ind) > 0:

                selectedreference = fits[median_sval_ind[0]]

                print('WRN! %sZ: selected best astrometric reference frame '
                      '(using S only) is %s' %
                      (datetime.utcnow().isoformat(), selectedreference))

                framejpg = fits_to_full_jpeg(
                    selectedreference,
                    out_fname=os.path.join(
                            os.path.dirname(selectedreference),
                        ('JPEG-ASTROMREF-%s.jpg' %
                         os.path.basename(selectedreference).strip('.fits.fz'))
                    )
                )

                arefinfo = {
                    'fits':selectedreference,
                    'jpg':framejpg,
                    'fistar':fistar[median_sval_ind[0]],
                    'fiphot':fiphot[median_sval_ind[0]],
                    'framekey':framekey[median_sval_ind[0]],
                    'sval':mfs[median_sval_ind[0]],
                    'dval':mfd[median_sval_ind[0]],
                    'bgv':mbg[median_sval_ind[0]],
                    'ndet':ngo[median_sval_ind[0]],
                    'comment':'astromref chosen using sval'
                }

            # if everything fails, do nothing
            else:

                print('ERR! %sZ: could not select a '
                      'good astrometric reference frame! all tests failed!' %
                      (datetime.utcnow().isoformat(), ))
                arefinfo = None


            # update the astromrefs table in the database if we found an
            # astromref frame, and copy it to the reference-frames directory
            # with the appropriate filename
            if arefinfo:

                # now that we have the astromref frame, copy it over to the
                # system-wide reference-images directory along with its JPEG
                # snapshot
                areftargetfits = (
                    'proj{projectid}-{field}-'
                    'ccd{ccd}-astromref-{origfname}.fits').format(
                        projectid=projectid,
                        field=field,
                        ccd=ccd,
                        origfname=os.path.splitext(
                            os.path.basename(arefinfo['fits'])
                        )[0]
                    )
                areftargetjpeg = areftargetfits.replace('.fits','.jpg')
                areftargetfistar = areftargetfits.replace('.fits','.fistar')
                areftargetfiphot = areftargetfits.replace('.fits','.fiphot')

                # copy the frame, jpeg, and fistar to the reference-frames dir
                shutil.copy(arefinfo['fits'],os.path.join(REFBASEDIR,
                                                                areftargetfits))
                shutil.copy(arefinfo['jpg'],os.path.join(REFBASEDIR,
                                                                areftargetjpeg))
                shutil.copy(arefinfo['fistar'],
                            os.path.join(REFBASEDIR, areftargetfistar))
                shutil.copy(arefinfo['fiphot'],
                            os.path.join(REFBASEDIR, areftargetfiphot))

                # now, update the astomrefs table in the database
                if overwrite:

                    query = (
                        "insert into astromrefs ("
                        "projectid, field, ccd, isactive, entryts, "
                        "framekey, fits, fistar, fiphot, jpeg, "
                        "mediansval, mediandval, medianbgv, ngoodobj, "
                        "comment"
                        ") values ("
                        "%s, %s, %s, %s, current_timestamp, "
                        "%s, %s, %s, %s, %s, "
                        "%s, %s, %s, %s, "
                        "%s"
                        ") on conflict on constraint astromrefs_pkey "
                        "do update set "
                        "entryts = current_timestamp, "
                        "framekey = %s, fits = %s, fistar = %s, fiphot = %s, "
                        "jpeg = %s, mediansval = %s, mediandval = %s, "
                        "medianbgv = %s, ngoodobj = %s, comment = %s"
                    )
                    params = (
                        projectid, field, ccd, makeactive,
                        arefinfo['framekey'], arefinfo['fits'],
                        arefinfo['fistar'], arefinfo['fiphot'],
                        arefinfo['jpg'], arefinfo['sval'], arefinfo['dval'],
                        arefinfo['bgv'],arefinfo['ndet'], arefinfo['comment'],
                        arefinfo['framekey'], arefinfo['fits'],
                        arefinfo['fistar'], arefinfo['fiphot'],
                        arefinfo['jpg'], arefinfo['sval'], arefinfo['dval'],
                        arefinfo['bgv'],arefinfo['ndet'], arefinfo['comment'],
                    )

                else:

                    query = (
                        "insert into astromrefs ("
                        "projectid, field, ccd, isactive, entryts, "
                        "framekey, fits, fistar, fiphot, jpeg, "
                        "mediansval, mediandval, medianbgv, ngoodobj, "
                        "comment"
                        ") values ("
                        "%s, %s, %s, %s, current_timestamp, "
                        "%s, %s, %s, %s, %s, "
                        "%s, %s, %s, %s, "
                        "%s"
                        ")"
                    )
                    params = (
                        projectid, field, ccd, makeactive,
                        arefinfo['framekey'], arefinfo['fits'],
                        arefinfo['fistar'], arefinfo['fiphot'],
                        arefinfo['jpg'], arefinfo['sval'], arefinfo['dval'],
                        arefinfo['bgv'],arefinfo['ndet'], arefinfo['comment']
                    )

                # execute the query to insert the astromref into the DB
                cursor.execute(query, params)
                database.commit()

                returnval = arefinfo

            # if we failed to find an astromref, do nothing
            else:

                returnval = None

        # if we didn't find any frames for this projectid-field-ccd combo,
        # complain and drop out
        else:

            print('ERR! %sZ: could not select a '
                  'good astrometric reference frame for %s, %s, %s!'
                  ' no frames exist' %
                  (datetime.utcnow().isoformat(), projectid, field, ccd))
            returnval = None

    except Exception as e:

        database.rollback()

        message = 'failed to get astromref for %s, %s, %s' % (projectid,
                                                              field,
                                                              ccd)
        print('EXC! %sZ: %s\nexception was: %s' %
               (datetime.utcnow().isoformat(),
                message, format_exc()) )
        returnval = None

    finally:

        cursor.close()
        if closedb:
            database.close()

    return returnval





def generate_astromref(fitsfiles,
                       makeactive=True,
                       field=None,
                       ccd=None,
                       projectid=None,
                       refdir=REFBASEDIR,
                       refinfo=REFINFO):

    '''This chooses an astrometry reference frame from the frames in fitfiles.

    writes the frame to refdir.

    ref frames have the following filename pattern:

    proj{projectid}-ccd{ccd}-{field}-astromref-{origfname}.fits

    if field, ccd, or projectid are None, these values are taken from the FITS
    file headers.

    updates the refinfo database.

    '''

    goodfits = [x for x in fitsfiles if os.path.exists(x)]

    if not goodfits:
        print('ERR! %sZ: no good FITS files found in input list' %
              (datetime.utcnow().isoformat(),))
        return

    # find the astromref
    astromref = ism.select_astromref_frame(
        fitsfiles,
        '1-*.fits',
    )

    # if an astromref was successfully found, then add its info to the DB
    if astromref:

        if field and ccd and projectid:

            frameinfo = {'field':field,
                         'ccd':ccd,
                         'projectid':projectid}

        else:

            # get the frame info
            frameelems = get_header_keyword_list(astromref['astromref'],
                                                 ['object',
                                                  'projid'])

            felems = FRAMEREGEX.findall(
                os.path.basename(astromref['astromref'])
            )

            if felems and felems[0]:

                ccd = felems[0][2]
                frameinfo = {'field':frameelems['object'],
                             'ccd':ccd,
                             'projectid':frameelems['projid']}

            else:

                print('ERR! %sZ: could not figure out CCD for astromref: %s' %
                      (datetime.utcnow().isoformat(), astromref['astromref']))
                return

            # now that we have the astromref frame, copy it over to the
            # system-wide reference-images directory along with its JPEG
            # snapshot
            areftargetfits = ('proj{projectid}-{field}-'
                              'ccd{ccd}-astromref-{origfname}.fits').format(
                                  projectid=frameinfo['projectid'],
                                  field=frameinfo['field'],
                                  ccd=frameinfo['ccd'],
                                  origfname=os.path.splitext(
                                      os.path.basename(astromref['astromref'])
                                  )[0]
                               )
            areftargetjpeg = areftargetfits.replace('.fits','.jpg')
            areftargetfistar = areftargetfits.replace('.fits','.fistar')

            # copy the frame, jpeg, and fistar to the reference-frames dir
            shutil.copy(astromref['astromref'],os.path.join(REFBASEDIR,
                                                            areftargetfits))
            shutil.copy(astromref['framejpg'],os.path.join(REFBASEDIR,
                                                            areftargetjpeg))
            shutil.copy(astromref['astromref'].replace('.fits','.fistar'),
                        os.path.join(REFBASEDIR, areftargetfistar))

            # now, put together the information and write to the refinfo sqlite

            query = ("insert into astromrefs "
                     "(field, projectid, ccd, isactive, unixtime, "
                     "framepath, jpegpath, sval, dval, bgv, ndet, "
                     "comment) values "
                     "(?, ?, ?, ?, ?, "
                     "?, ?, ?, ?, ?, ?, "
                     "?)")
            params = (frameinfo['field'],
                      frameinfo['projectid'],
                      frameinfo['ccd'],
                      1 if makeactive else 0,
                      time.time(),

                      os.path.join(REFBASEDIR,areftargetfits),
                      os.path.join(REFBASEDIR,areftargetjpeg),
                      astromref['sval'],
                      astromref['dval'],
                      astromref['bgv'],
                      astromref['ndet'],

                      (astromref['comment'] +
                       '; original: %s' % astromref['astromref']))

            db = sqlite3.connect(
                refinfo,
                detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES
            )
            cur = db.cursor()

            try:

                astromref.update(frameinfo)
                cur.execute(query, params)
                db.commit()

                print('%sZ: using astromref %s for '
                      'field %s, ccd %s, project id %s, database updated.' %
                      (datetime.utcnow().isoformat(),
                       astromref['astromref'],
                       astromref['field'],
                       astromref['ccd'],
                       astromref['projectid']))

                returnval = astromref

            except Exception as e:

                print('ERR! %sZ: could not update refinfo DB! error was: %s' %
                      (datetime.utcnow().isoformat(), e))
                returnval = None
                db.rollback()

            db.close()

    # otherwise, do nothing
    else:

        print('ERR! %sZ: could not find an astromref frame' %
              (datetime.utcnow().isoformat(),))
        returnval = None


    return returnval




def get_astromref(projectid, field, ccd, refinfo=REFINFO):
    '''This finds the reference frame for the field, projectid, and ccd
    combination using the TM-refinfo.sqlite database.


    '''

    db = sqlite3.connect(
        refinfo,
        detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES
    )
    cur = db.cursor()

    query = ('select field, projectid, ccd, unixtime, '
             'framepath, jpegpath, sval, dval, bgv, ndet, comment '
             'from astromrefs where '
             'projectid = ? and field = ? and ccd = ? and '
             'isactive = 1')
    params = (projectid, field, ccd)

    try:

        cur.execute(query, params)
        rows = cur.fetchone()

        astromref = {x:y for (x,y) in zip(('field','projectid','ccd',
                                           'unixtime','framepath','jpegpath',
                                           'sval','dval','bgv',
                                           'ndet','comment'),rows)}

        returnval = astromref

    except Exception as e:

        print('ERR! %sZ: could not get astromref info '
              'from DB! error was: %s' %
              (datetime.utcnow().isoformat(), e))
        returnval = None

    db.close()

    return returnval



def frames_astromref_worker(task):
    '''
    This is the parallel worker for frames_to_astromref.

    task[0] = fits file
    task[1] = outdir
    task[2] = refinfo
    task[3] = warpcheck
    task[4] = warpcheck kwargs {'threshold', 'margins'}

    '''

    try:

        frame, outdir, refinfo, warpcheck, warpcheckkwargs = task

        # figure out this frame's field, ccd, and projectid
        frameelems = get_header_keyword_list(frame,
                                             ['object',
                                              'projid'])

        felems = FRAMEREGEX.findall(
            os.path.basename(frame)
        )

        framefistar = frame.replace('.fits','.fistar')

        if felems and felems[0] and os.path.exists(framefistar):

            ccd = felems[0][2]
            frameinfo = {'field':frameelems['object'],
                         'ccd':ccd,
                         'projectid':frameelems['projid']}


            # find this frame's associated active astromref
            framearef = get_astromref(frameinfo['projectid'],
                                      frameinfo['field'],
                                      frameinfo['ccd'],
                                      refinfo=refinfo)
            areffistar = framearef['framepath'].replace('.fits','.fistar')

            # calculate the shift and write the itrans back to the frame's
            # directory
            shifted_fistar, shifted_itrans = ism.astromref_shift_worker(
                (framefistar, areffistar, outdir)
            )

            # if the shift calculation is successful, shift the image itself
            if shifted_itrans and os.path.exists(shifted_itrans):

                frame_to_shift, shifted_frame = ism.frame_to_astromref_worker(
                    (frame, None, None)
                )

                if shifted_frame and os.path.exists(shifted_frame):

                    # check if the frame has warped too much after the shift,
                    # these frames look like they're folding into/out of the
                    # z-direction. we need to throw these away.
                    if warpcheck:

                        notwarped, warpinfo = check_frame_warping(
                            shifted_frame,
                            **warpcheckkwargs
                        )

                        # if image is OK, return it
                        if notwarped:

                            print('%sZ: SHIFT OK %s -> %s' %
                                  (datetime.utcnow().isoformat(),
                                   frame, shifted_frame))

                            return frame, shifted_frame

                        # otherwise, move it to the badframes subdir and mark it
                        # as warped
                        else:

                            badframesdir = os.path.join(os.path.dirname(frame),
                                                        'badframes')
                            if not os.path.exists(badframesdir):
                                os.mkdir(badframesdir)

                            # find all the components of this frame and move
                            # them to the badframes subdir
                            badframeglob = glob.glob(
                                os.path.join(
                                    os.path.dirname(shifted_frame),
                                    '*%s*.*' % (
                                        os.path.splitext(
                                            os.path.basename(frame)
                                        )[0]
                                    )
                                )
                            )

                            for x in badframeglob:
                                shutil.move(x, badframesdir)

                            print('WRN! %sZ: SHIFT HAS WARPED '
                                  'IMAGE, moved %s and metadata to %s' %
                                  (datetime.utcnow().isoformat(),
                                   frame, badframesdir))

                            return frame, None

                    # if we're not checking for warps, just check if the image
                    # was shifted fine
                    else:

                        print('%sZ: SHIFT OK %s -> %s' %
                              (datetime.utcnow().isoformat(),
                               frame, shifted_frame))

                        return frame, shifted_frame


                else:

                    print('ERR! %sZ: SHIFT OPERATION FAILED for %s' %
                          (datetime.utcnow().isoformat(), frame))
                    return frame, None

            else:

                print('ERR! %sZ: SHIFT CALCULATION FAILED for %s' %
                      (datetime.utcnow().isoformat(), frame))
                return frame, None

        else:

            print('ERR! %sZ: could not figure out '
                  'CCD info or fistar for frame: %s' %
                  (datetime.utcnow().isoformat(), frame))
            return frame, None

    except Exception as e:

        print('ERR! %sZ: could not shift frame %s to astromref, error was: %s' %
              (datetime.utcnow().isoformat(), frame, e))
        return frame, None



def framelist_make_xtrnsfits(fitsfiles,
                             outdir=None,
                             refinfo=REFINFO,
                             warpcheck=True,
                             warpthreshold=2000.0,
                             warpmargins=100,
                             nworkers=16,
                             maxworkertasks=1000):
    '''This calculates the shifts between frames in fitsfiles and the appropriate
    astromref for the projectid, field and CCD, then shifts each frame to the
    astromref's coordinate system, generating -xtrns.fits files.

    '''

    print('%sZ: %s files to process' %
          (datetime.utcnow().isoformat(), len(fitsfiles)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, outdir, refinfo, warpcheck,
              {'threshold':warpthreshold, 'margins':warpmargins})
             for x in fitsfiles if os.path.exists(x)]

    # fire up the pool of workers
    results = pool.map(frames_astromref_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {x:y for (x,y) in results}



##################################
## PHOTOMETRIC REFERENCE FRAMES ##
##################################

def generate_photref_candidates_from_xtrns(fitsfiles,
                                           minframes=50,
                                           maxhourangle=3.0,
                                           maxmoonphase=25.0,
                                           maxmoonelev=0.0,
                                           maxzenithdist=30.0,
                                           maxbackgroundstdev=10.0,
                                           maxbackgroundmedian=1000.0,
                                           forcecollectinfo=False,
                                           nworkers=8,
                                           maxworkertasks=1000):
    '''This uses ism.select_photref_frames run on fitsfiles to get photref
    candidates.

    fitsfiles must be a list of frames, which have been already transformed to
    the astromref, and are all from a single projectid, ccd, field combination
    for this operation to make sense.

    '''

    # first, get all the info from these fits files.
    frameinfo = fitslist_frameinfo(fitsfiles,
                                   forcecollectinfo=False,
                                   nworkers=nworkers,
                                   maxworkertasks=maxworkertasks)

    # this is the cachekey used to store the photref selection info
    cachekey = '%s-%i-%.1f-%.1f-%.1f-%.1f-%.1f-%.1f' % (repr(fitsfiles),
                                                        minframes,
                                                        maxhourangle,
                                                        maxmoonphase,
                                                        maxmoonelev,
                                                        maxzenithdist,
                                                        maxbackgroundstdev,
                                                        maxbackgroundmedian)
    cachekey = md5(cachekey).hexdigest()
    cachedir = os.path.join(FRAMEINFOCACHEDIR,'TM-photref-%s' % cachekey)
    cacheinfofile = os.path.join(cachedir, 'selection-info.pkl.gz')

    # get the data from the cache if it exists and we're allowed to use it
    if ((not forcecollectinfo) and
        os.path.exists(cachedir) and
        os.path.exists(cacheinfofile)):

        with gzip.open(cacheinfofile) as infd:
            photrefinfo = pickle.load(infd)

        print('%sZ: candidate photref JPEGs in: %s, '
              'cached photrefinfo from: %s' %
              (datetime.utcnow().isoformat(), cachedir, cacheinfofile))

        return photrefinfo

    ## OTHERWISE, RUN THE FULL PROCESS ##

    # then, apply our conditions to these fits files to generate a list of
    # photref candidates
    # filter on hour angle
    haind = np.fabs(frameinfo['hourangle']) < maxhourangle
    print('%sZ: %s frames with hour angle < %s' %
          (datetime.utcnow().isoformat(),
           len(np.where(haind)[0]),
           maxhourangle))

    # get dark nights
    moonind = ((np.fabs(frameinfo['moonphase']) < maxmoonphase) |
               (frameinfo['moonelev'] < maxmoonelev))
    print('%sZ: %s frames with moon phase < %s or moon elev < %s' %
          (datetime.utcnow().isoformat(),
           len(np.where(moonind)[0]),
           maxmoonphase,
           maxmoonelev))

    # get low zenith distance nights
    zenithind = frameinfo['zenithdist'] < maxzenithdist
    print('%sZ: %s frames with zenith distance < %s' %
          (datetime.utcnow().isoformat(),
           len(np.where(zenithind)[0]),
           maxzenithdist))

    # get nights with background stdev < max_bgv_stdev (to possibly remove
    # cloudy nights)
    backgroundstdevind = frameinfo['stdsrcbgv'] < maxbackgroundstdev
    print('%sZ: %s frames with background stdev < %s' %
          (datetime.utcnow().isoformat(),
           len(np.where(backgroundstdevind)[0]),
           maxbackgroundstdev))

    # get nights with background median < maxbackgroundmedian (to possibly
    # remove cloudy nights)
    backgroundmedind = frameinfo['medsrcbgv'] < maxbackgroundmedian
    print('%sZ: %s frames with background median < %s' %
          (datetime.utcnow().isoformat(),
           len(np.where(backgroundmedind)[0]),
           maxbackgroundmedian))

    # this is the final operating set of frames that will be sorted for the
    # following tests
    selectind = haind & moonind & zenithind & backgroundstdevind

    selected_frames = frameinfo['frames'][selectind]
    selected_ngoodobj = frameinfo['ngoodobjects'][selectind]

    selected_medmagerr = frameinfo['medmagerr'][selectind]
    selected_magerrmad = frameinfo['magerrmad'][selectind]

    selected_medsrcbgv = frameinfo['medsrcbgv'][selectind]
    selected_stdsrcbgv = frameinfo['stdsrcbgv'][selectind]

    selected_medsvalue = frameinfo['medsval'][selectind]
    selected_meddvalue = frameinfo['meddval'][selectind]

    print('\n%sZ: selected %s frames with acceptable '
          'HA, Z, moon phase, background, and elevation '
          'for further filtering...\n' %
          (datetime.utcnow().isoformat(), len(selected_frames)))

    # we select in the following order
    # 1. D closest to 0
    # 2. largest S

    # then we get filter out any images that have background >
    # maxbackgroundmedian and backgroundstdev > maxbackgroundstdev

    # first sort selector
    stage1_sort_ind = (np.argsort(selected_medsvalue))[::-1]

    stage1_frames = selected_frames[stage1_sort_ind[:2*minframes]]
    stage1_median_bgv = selected_medsrcbgv[stage1_sort_ind[:2*minframes]]
    stage1_stdev_bgv = selected_stdsrcbgv[stage1_sort_ind[:2*minframes]]
    stage1_svalue = selected_medsvalue[stage1_sort_ind[:2*minframes]]
    stage1_dvalue = selected_meddvalue[stage1_sort_ind[:2*minframes]]

    # next, sort by roundest stars
    stage2_sort_ind = (np.argsort(np.fabs(stage1_dvalue)))

    stage2_frames = stage1_frames[stage2_sort_ind]
    stage2_median_bgv = stage1_median_bgv[stage2_sort_ind]
    stage2_stdev_bgv = stage1_stdev_bgv[stage2_sort_ind]
    stage2_svalue = stage1_svalue[stage2_sort_ind]
    stage2_dvalue = stage1_dvalue[stage2_sort_ind]

    final_bgvmed_ind = stage2_median_bgv < maxbackgroundmedian
    final_bgvstd_ind = stage2_stdev_bgv < maxbackgroundstdev
    final_selector_ind = final_bgvmed_ind & final_bgvstd_ind

    final_frames = stage2_frames[final_selector_ind][:minframes]
    final_median_bgv = stage2_median_bgv[final_selector_ind][:minframes]
    final_stdev_bgv = stage2_stdev_bgv[final_selector_ind][:minframes]
    final_svalues = stage2_svalue[final_selector_ind][:minframes]
    final_dvalues = stage2_dvalue[final_selector_ind][:minframes]

    print('%sZ: selected %s final frames as photref' %
          (datetime.utcnow().isoformat(), len(final_frames)))

    # the master photref is the frame we'll convolve all of the rest of the
    # photrefs to. it's the softest of these frames
    try:

        candidate_master_photref = final_frames[np.nanargmin(final_svalues)]
        final_jpegs = []

        # make JPEGs of the selected photref frames and copy them to the
        # cachedir
        if not os.path.exists(cachedir):
            print('WRN! %sZ: making new photref cache directory: %s' %
                  (datetime.utcnow().isoformat(), cachedir))
            os.mkdir(cachedir)

        for final_frame in final_frames:

            framejpg = fits_to_full_jpeg(
                final_frame,
                out_fname=os.path.join(
                    cachedir,
                    ('JPEG-PHOTREF-%s.jpg' %
                     os.path.basename(final_frame).rstrip('.fits.fz'))
                    )
                )
            final_jpegs.append(framejpg)

        photrefinfo = {
            'framelist':fitsfiles,
            'frameinfo':frameinfo,
            'cachekey':cachekey,
            'minframes':minframes,
            'maxhourangle':maxhourangle,
            'maxmoonphase':maxmoonphase,
            'maxmoonelev':maxmoonelev,
            'maxzenithdist':maxzenithdist,
            'maxbackgroundstdev':maxbackgroundstdev,
            'maxbackgroundmedian':maxbackgroundmedian,
            'masterphotref':os.path.abspath(candidate_master_photref),
            'photrefs':[os.path.abspath(x) for x in final_frames],
            'photrefjpegs':final_jpegs
        }

        # dump the photrefinfo to a pickle
        with gzip.open(cacheinfofile,'wb') as outfd:
            pickle.dump(photrefinfo, outfd, pickle.HIGHEST_PROTOCOL)

        print('%sZ: candidate photref JPEGs in: %s, photrefinfo dumped to: %s' %
              (datetime.utcnow().isoformat(), cachedir, cacheinfofile))

        return photrefinfo

    except Exception as e:

        print('ERR! %sZ: selection failed, some criteria '
              'may be too strict for this frame list' %
              (datetime.utcnow().isoformat()))

        return {'framelist':fitsfiles,
                'frameinfo':frameinfo,
                'cachekey':cachekey,
                'maxhourangle':maxhourangle,
                'maxmoonphase':maxmoonphase,
                'maxmoonelev':maxmoonelev,
                'maxzenithdist':maxzenithdist,
                'maxbackgroundstdev':maxbackgroundstdev,
                'maxbackgroundmedian':maxbackgroundmedian,
                'masterphotref':None,
                'photrefs':None,
                'photrefjpegs':None}



def amend_candidate_photrefs(photrefinfo):
    '''This is an interactive way to update masterphotref, photrefs, and
    photrefjpegs after reviewing them.

    This will automatically update the photrefinfo cache.

    '''

    cachekey = photrefinfo['cachekey']
    cachedir = os.path.join(FRAMEINFOCACHEDIR,'TM-photref-%s' % cachekey)
    cacheinfofile = os.path.join(cachedir, 'selection-info.pkl.gz')

    print('reviewing photrefinfo for %s\n' % cachedir)

    # now deal with the photrefs:
    print('-- CANDIDATE PHOTREFS --\n')

    initialphotrefs = sorted(photrefinfo['photrefs'][::])
    initialphotrefjpegs = sorted(photrefinfo['photrefjpegs'][::])

    for frame, jpeg in zip(initialphotrefs, initialphotrefjpegs):

        breakloop = False

        photref_prompt = (
            'photref = %s, jpeg = %s\n'
            '[ENTER] to keep this, or [x] to remove: ' %
            (frame, os.path.basename(jpeg))
        )

        while not breakloop:

            photref_check = raw_input(photref_prompt)

            if photref_check and photref_check == 'x':

                photrefinfo['photrefs'].remove(frame)
                photrefinfo['photrefjpegs'].remove(jpeg)
                os.remove(jpeg)

                print('REMOVED photref %s' % frame)
                breakloop = True

            elif not photref_check:
                breakloop = True

    print('\nfinal photrefs set to:')
    for frame in photrefinfo['photrefs']:
        print(frame)

    # next, update the masterphotref
    masterphotref_prompt = (
        'current masterphotref = %s\n'
        '[ENTER] to keep this, or new masterphot: ' %
        photrefinfo['masterphotref']
    )

    breakloop = False

    print('\n-- MASTERPHOTREF --\n')

    # loop until masterphotref is satisfied
    while not breakloop:

        masterphotref_amendment = raw_input(masterphotref_prompt)

        if masterphotref_amendment and os.path.exists(masterphotref_amendment):

            photrefinfo['masterphotref'] = masterphotref_amendment[::]

            masterphotref_prompt = (
                'new masterphotref = %s\n'
                '[ENTER] to keep this, or new masterphot: ' %
                photrefinfo['masterphotref']
            )

        elif masterphotref_amendment and not os.path.exists(masterphotref_amendment):

            masterphotref_prompt = (
                'masterphotref = %s does not exist\n'
                'new masterphot: ' %
                masterphotref_amendment
            )

        elif not masterphotref_amendment:
            breakloop = True

    print('\nmasterphotref set to %s' % photrefinfo['masterphotref'])

    # update the cache info file
    print('\nupdating photref cached selection-info pickle...')

    # dump the photrefinfo to a pickle
    with gzip.open(cacheinfofile,'wb') as outfd:
        pickle.dump(photrefinfo, outfd, pickle.HIGHEST_PROTOCOL)

    print('%sZ: candidate photref JPEGs in: %s, photrefinfo dumped to: %s' %
          (datetime.utcnow().isoformat(), cachedir, cacheinfofile))

    return photrefinfo



def generate_combined_photref(
        photrefinfo,
        photreftype,
        makeactive=True,
        field=None,
        ccd=None,
        projectid=None,
        refdir=REFBASEDIR,
        refinfo=REFINFO,
        fovcatdir=FIELDCAT_DIR,
        combinemethod='median',
        kernelspec='b/4;i/4;d=4/4',
        ccdgain=None,
        zeropoint=None,
        ccdexptime=None,
        extractsources=True,
        astrometrysrcthreshold=25000,
        apertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
        framewidth=None,
        searchradius=8.0,
        nworkers=8,
        maxworkertasks=1000
):
    '''This generates a combined photref from photref target and candidates and
    updates the TM-refinfo.sqlite database.

    Use this after reviewing the results from
    generate_photref_candidates_from_xtrns function above. Amend the
    photrefinfo['masterphotref'], photrefinfo['photrefs'], and
    photrefinfo['photrefjpegs'] arrays as needed using the
    amend_candidate_photrefs function above.

    photreftype is the type of the combined photref produced. it must be one of
    the following strings:

    'oneframe' -> single HATPI frame
    'onehour' -> up to 120 HATPI frames
    'onenight' -> up to 960 HATPI frames

    updates photrefinfo with the following dict and keys:

    'combinedphotref':{'frame': -> combined photref frame path
                       'jpeg' -> combined photref jpeg path
                       'cmrawphot' -> cmrawphot file path
                       'regfile' -> convolution registration file path
                       'combinemethod'- > combine type
                       'reftype' -> combined photref type
                       'phottype' -> either 're-extracted' or 'cat-projected'
                       'photaps' -> photometry apertures for combined photref
                       'fovcat' -> fovcat file used for photometry
                       'kernelspec' -> convolution kernel specs}

    and updates the cached selection-info pickle file as well.

    the output combined photref frame, jpeg, cmrawphot (and byproducts) go to
    the REFBASEDIR, using the following prototype for the filename:

    {REFBASEDIR}/proj{projid}-{field}-ccd{ccd}-combinedphotref-{photreftype}.XXX

    '''

    # get the field, ccd, projectid first (from the convolvetarget =
    # masterphotref)

    masterphotref = photrefinfo['masterphotref']

    frameelems = get_header_keyword_list(masterphotref,
                                         ['object','projid'])

    felems = FRAMEREGEX.findall(
        os.path.basename(masterphotref)
    )

    if felems and felems[0]:

        ccd = felems[0][2]
        masterphotrefinfo = {'field':frameelems['object'],
                             'ccd':int(ccd),
                             'projectid':frameelems['projid']}

    else:

        print('ERR! %sZ: could not figure out CCD for masterphotref: %s' %
              (datetime.utcnow().isoformat(), masterphotref))
        return

    # make the convolution registration file

    photreffname = ('proj{projid}-{field}-ccd{ccd}'
                    '-combinedphotref-{photreftype}.{fileext}')

    regfpath = os.path.join(
        refdir,
        photreffname.format(
            projid=masterphotrefinfo['projectid'],
            field=masterphotrefinfo['field'],
            ccd=masterphotrefinfo['ccd'],
            photreftype=photreftype,
            fileext='reg'
        )
    )

    masterphotref_fistar = masterphotref.replace('-xtrns.fits','.fistar')

    if not os.path.exists(masterphotref_fistar):
        print('ERR! %sZ: no fistar available for masterphotref: %s' %
              (datetime.utcnow().isoformat(), masterphotref))
        return

    # conv registration file
    ism.genreg(masterphotref_fistar,
               regfpath)

    if not os.path.exists(regfpath):
        print('ERR! %sZ: could not make regfile for masterphotref: %s' %
              (datetime.utcnow().isoformat(), masterphotref))
        return

    # convolve all candidate photrefs to the masterphotref
    convresult = ism.convolve_photref_frames(photrefinfo['photrefs'],
                                             masterphotref,
                                             regfpath,
                                             kernelspec=kernelspec,
                                             nworkers=nworkers,
                                             maxworkertasks=maxworkertasks)

    # get the convolved photref frames
    convphotrefs = [convresult[x] for x in convresult
                    if os.path.exists(convresult[x])]

    if len(convphotrefs) == 0:
        print('ERR! %sZ: convolution of photrefs to masterphotref: %s failed' %
              (datetime.utcnow().isoformat(), masterphotref))
        return

    # combine all the convolved photrefs into a single combinedphotref, using
    # combinemethod

    # the output combinedphotref path
    combinedphotrefpath = os.path.join(
        refdir,
        photreffname.format(
            projid=masterphotrefinfo['projectid'],
            field=masterphotrefinfo['field'],
            ccd=masterphotrefinfo['ccd'],
            photreftype=photreftype,
            fileext='fits'
        )
    )

    combinedphotref = ism.combine_frames(convphotrefs,
                                         combinedphotrefpath,
                                         combinemethod=combinemethod)

    if not (combinedphotref[1] and os.path.exists(combinedphotref[1])):
        print('ERR! %sZ: combining conv photrefs '
              'into masterphotref: %s failed' %
              (datetime.utcnow().isoformat(), masterphotref))
        return

    # rearrange the returned combinedphotref filename
    combinedphotref = combinedphotref[1]

    # find the fovcat file for the field, ccd, projectid, photreftype combo
    # photreftype = 'oneframe' -> default field-gri.catalog
    # photreftype = 'onehour' -> default field-gri-18.0.catalog
    # photreftype = 'onenight' -> default field-gri-20.0.catalog

    fovcat_template = '{field}{bandspec}{magspec}.catalog'

    if photreftype == 'oneframe':
        photref_fovcatpath = os.path.join(
            fovcatdir,
            fovcat_template.format(
                field=masterphotrefinfo['field'],
                bandspec='-gri',
                magspec=''
                )
            )
    elif photreftype == 'onehour':
        photref_fovcatpath = os.path.join(
            fovcatdir,
            fovcat_template.format(
                field=masterphotrefinfo['field'],
                bandspec='-gri',
                magspec='-18.0'
                )
            )
    elif photreftype == 'onenight':
        photref_fovcatpath = os.path.join(
            fovcatdir,
            fovcat_template.format(
                field=masterphotrefinfo['field'],
                bandspec='-gri',
                magspec='-20.0'
                )
            )
    else:
        print('ERR! %sZ: unknown photreftype: %s specified '
              'can\'t continue...' %
              (datetime.utcnow().isoformat(), photreftype))
        return

    if not os.path.exists(photref_fovcatpath):
        print('ERR! %sZ: no FOV catalog available '
              'for field %s, photreftype %s, '
              'can\'t do photometry on combinedphotref %s' %
              (datetime.utcnow().isoformat(),
               masterphotref, photreftype, combinedphotref))
        return


    # run photometry on the combinedphotref and generate a cmrawphot file
    cphotref_photometry = ism.photometry_on_combined_photref(
        combinedphotref,
        photref_fovcatpath,
        masterphotrefinfo['ccd'],
        ccdgain=ccdgain,
        zeropoint=zeropoint,
        ccdexptime=ccdexptime,
        extractsources=extractsources,
        apertures=apertures,
        framewidth=framewidth,
        searchradius=searchradius,
        astrometrysourcethreshold=astrometrysrcthreshold,
    )

    if not (cphotref_photometry and
            cphotref_photometry[1] and
            os.path.exists(cphotref_photometry[1])):
        print('ERR! %sZ: photometry failed for combinedphotref %s '
              'can\'t continue...' %
              (datetime.utcnow().isoformat(), combinedphotref))
        return

    # update the cache photref selection-info.pkl.gz file
    combinedphotrefinfo = {
        'reftype':photreftype,
        'frame':combinedphotref,
        'jpeg':combinedphotref.replace('.fits',
                                       '.jpg').replace('proj',
                                                       'JPEG-COMBINED-proj'),
        'cmrawphot':cphotref_photometry[1],
        'regfile':regfpath,
        'combinemethod':combinemethod,
        'kernelspec':kernelspec,
        'phottype':'re-extracted' if extractsources else 'cat-projected',
        'photaps':apertures,
        'fovcat':photref_fovcatpath,
    }
    photrefinfo['combinedphotref'] = combinedphotrefinfo

    cachekey = photrefinfo['cachekey']
    cachedir = os.path.join(FRAMEINFOCACHEDIR,'TM-photref-%s' % cachekey)
    cacheinfofile = os.path.join(cachedir, 'selection-info.pkl.gz')

    with gzip.open(cacheinfofile, 'wb') as outfd:
        print('%sZ: combined photref JPEG: %s, photrefinfo updated: %s' %
              (datetime.utcnow().isoformat(),
               photrefinfo['combinedphotref']['jpeg'],
               cacheinfofile))
        pickle.dump(photrefinfo, outfd, pickle.HIGHEST_PROTOCOL)


    # update the TM-refinfo.sqlite database

    # first, get the frame info from the combinedphotref
    _, photref_frameinfo = get_frame_info(combinedphotref)

    if not photref_frameinfo:
        print('ERR! %sZ: could not extract frame info from combinedphotref %s' %
              (datetime.utcnow().isoformat(), combinedphotref))
        return


    query = ("insert into photrefs "
             "(field, projectid, ccd, photreftype, isactive, unixtime, "
             "framepath, jpegpath, "
             "convolvetarget, convolveregpath, cmrawphotpath, "
             "target_zenithdist, target_moondist, target_moonelev, "
             "target_moonphase, target_hourangle, target_ndet, "
             "target_medmagerr, target_magerrmad, target_medsrcbgv, "
             "target_stdsrcbgv, target_medsval, target_meddval, "
             "photrefinfo) values "
             "(?, ?, ?, ?, ?, ?, "
             "?, ?, "
             "?, ?, ?, "
             "?, ?, ?, "
             "?, ?, ?, "
             "?, ?, ?, "
             "?, ?, ?, "
             "?)")
    params = (
        masterphotrefinfo['field'],
        masterphotrefinfo['projectid'],
        masterphotrefinfo['ccd'],
        photreftype,
        1 if makeactive else 0,
        time.time(),

        photrefinfo['combinedphotref']['frame'],
        photrefinfo['combinedphotref']['jpeg'],

        masterphotref,
        photrefinfo['combinedphotref']['regfile'],
        photrefinfo['combinedphotref']['cmrawphot'],

        photref_frameinfo['zenithdist'],
        photref_frameinfo['moondist'],
        photref_frameinfo['moonelev'],

        photref_frameinfo['moonphase'],
        photref_frameinfo['hourangle'],
        photref_frameinfo['ngoodobjects'],

        photref_frameinfo['medmagerr'],
        photref_frameinfo['magerrmad'],
        photref_frameinfo['medsrcbgv'],

        photref_frameinfo['stdsrcbgv'],
        photref_frameinfo['medsval'],
        photref_frameinfo['meddval'],

        json.dumps(photrefinfo['combinedphotref'],ensure_ascii=True)
    )

    db = sqlite3.connect(
        refinfo,
        detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES
    )
    cur = db.cursor()

    try:

        cur.execute(query, params)
        db.commit()

        print('%sZ: will use combinedphotref %s for '
              'field %s, ccd %s, project id %s, database updated.' %
              (datetime.utcnow().isoformat(),
               combinedphotref,
               masterphotrefinfo['field'],
               masterphotrefinfo['ccd'],
               masterphotrefinfo['projectid']))

    except Exception as e:

        print('ERR! %sZ: could not update refinfo DB! error was: %s' %
              (datetime.utcnow().isoformat(), e))
        db.rollback()

    db.close()

    # return the updated photrefinfo dict
    return photrefinfo



def get_combined_photref(projectid,
                         field,
                         ccd,
                         photreftype,
                         refinfo=REFINFO):
    '''This gets the combined photref for the given combo of projid, field, ccd.

    Used for the convsubphot functions below.

    '''

    db = sqlite3.connect(
        refinfo,
        detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES
    )
    cur = db.cursor()

    query = (
        'select field,projectid,ccd,photreftype,unixtime,'
        'framepath,jpegpath,convolvetarget,convolveregpath,'
        'cmrawphotpath,target_zenithdist,target_moondist,'
        'target_moonelev,target_moonphase,target_hourangle,'
        'target_ndet,target_medmagerr,target_magerrmad,'
        'target_medsrcbgv,target_stdsrcbgv,target_medsval,'
        'target_meddval,photrefinfo from photrefs where '
        '(isactive = 1) and '
        '(projectid = ?) and '
        '(ccd = ?) and '
        '(field = ?) and '
        '(photreftype = ?)'
    )
    params = (projectid, ccd, field, photreftype)

    cur.execute(query, params)

    try:

        cur.execute(query, params)
        rows = cur.fetchone()

        cphotref = {x:y for (x,y) in zip(('field','projectid','ccd',
                                          'photreftype','unixtime',
                                          'framepath','jpegpath',
                                          'convolvetarget','convolveregpath',
                                          'cmrawphotpath',
                                          'target_zenithdist',
                                          'target_moondist',
                                          'target_moonelev',
                                          'target_moonphase',
                                          'target_hourangle',
                                          'target_ndet',
                                          'target_medmagerr',
                                          'target_magerrmad',
                                          'target_medsrcbgv',
                                          'target_stdsrcbgv',
                                          'target_medsval',
                                          'target_meddval',
                                          'photrefinfo'),rows)}

        # load the JSON string for photrefinfo
        cphotref['photrefinfo'] = json.loads(cphotref['photrefinfo'])

        returnval = cphotref

    except Exception as e:

        print('ERR! %sZ: could not get combinedphotref info '
              'from DB! error was: %s' %
              (datetime.utcnow().isoformat(), e))
        returnval = None

        raise


    db.close()
    return returnval



##################################
## IMAGE SUBTRACTION PHOTOMETRY ##
##################################

def xtrnsfits_convsubphot_worker(task):
    '''This is a parallel worker for framelist_convsubphot_photref below.

    task[0] = xtrnsfits file
    task[1] = photreftype to use <"oneframe"|"onehour"|"onenight">
    task[2] = outdir
    task[3] = kernelspec
    task[4] = reversesubtract boolean
    task[5] = findnewobjects boolean
    task[6] = photdisjointradius
    task[7] = refinfo
    task[8] = lcapertures

    TODO: think about combining this with the database import step below.
          - write the iphot to /dev/shm
          - use pg's \copy protocol to quickly import it into the database
          - once committed, remove the /dev/shm file

    TODO: we're now planning on writing photometry directly to the pg
    database. we'll use partitioned tables created automatically per week to
    hold the photometry. the partition key will be the framejd converted to a pg
    UTC timestamp.

    '''

    (frame, photreftype, outdir,
     kernelspec, reversesubtract, findnewobjects, photdisjointradius,
     refinfo, lcapertures) = task

    try:

        # first, figure out the input frame's projid, field, and ccd
        frameelems = get_header_keyword_list(frame,
                                             ['object',
                                              'projid'])
        felems = FRAMEREGEX.findall(
            os.path.basename(frame)
        )
        field, ccd, projectid = (frameelems['object'],
                                 int(felems[0][2]),
                                 frameelems['projid'])

        # then, find the associated combined photref frame, regfile, cmrawphot
        cphotref = get_combined_photref(projectid, field, ccd, photreftype,
                                        refinfo=refinfo)
        cphotref_frame = cphotref['framepath']
        cphotref_reg = cphotref['convolveregpath']
        cphotref_cmrawphot = cphotref['cmrawphotpath']

        # do the subtraction (take care of reversesubtract here)
        _, convsub = ism.subframe_convolution_worker(
            (frame, cphotref_frame, cphotref_reg,
             kernelspec, outdir, reversesubtract, photreftype)
        )

        if not (convsub and os.path.exists(convsub)):
            print('ERR! %sZ: convulution and subtraction failed on frame %s' %
                  (datetime.utcnow().isoformat(), frame))
            return frame, None

        # find new objects in the subtracted frame if told to do so
        if findnewobjects:
            pass

        # find matching kernel, itrans, and xysdk files for each subtracted
        # frame

        # generate the convsubfits hash
        convsubhash = get_convsubfits_hash(
            photreftype,
            ('reverse' if os.path.basename(convsub).startswith('rsub')
             else 'normal'),
            kernelspec
        )

        frameinfo = FRAMEREGEX.findall(os.path.basename(convsub))
        kernel = '%s-%s-%s_%s-xtrns.fits-kernel' % (photreftype,
                                                    convsubhash,
                                                    frameinfo[0][0],
                                                    frameinfo[0][1],
                                                    frameinfo[0][2])
        kernel = os.path.abspath(os.path.join(os.path.dirname(convsub),kernel))

        itrans = '%s-%s_%s.itrans' % (frameinfo[0][0],
                                      frameinfo[0][1],
                                      frameinfo[0][2])
        itrans = os.path.abspath(os.path.join(os.path.dirname(convsub),itrans))

        xysdk = '%s-%s_%s.xysdk' % (frameinfo[0][0],
                                    frameinfo[0][1],
                                    frameinfo[0][2])
        xysdk = os.path.abspath(os.path.join(os.path.dirname(convsub),xysdk))


        # write the photometry file to /dev/shm by default
        if outdir is None:
            outdir = '/dev/shm'

        _, subphot = ism.subframe_photometry_worker(
            (convsub, cphotref_cmrawphot, photdisjointradius,
             kernel, itrans, xysdk, outdir,
             photreftype, kernelspec, lcapertures)
        )

        if subphot and os.path.exists(subphot):

            print('%sZ: CONVSUBPHOT OK: frame %s, '
                  'subtracted frame %s, photometry file %s' %
                  (datetime.utcnow().isoformat(), frame, convsub, subphot))

            return frame, (convsub, subphot)

        else:

            print('%sZ: CONVSUBPHOT FAILED: frame %s' %
                  (datetime.utcnow().isoformat(), frame))
            return frame, None

    except Exception as e:

        print('ERR! %sZ: could not do convsubphot on frame %s, error was: %s' %
              (datetime.utcnow().isoformat(), frame, e))

        return frame, None



def xtrnsfits_convsubphot(xtrnsfits,
                          photreftype,
                          outdir=None,
                          refinfo=REFINFO,
                          reversesubtract=True,
                          kernelspec='b/4;i/4;d=4/4',
                          lcapertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                          photdisjointradius=2,
                          findnewobjects=False,
                          nworkers=16,
                          maxworkertasks=1000):
    '''This convolves, subtracts, and does photometry of known photref sources for
    all FITS files in the xtrnsfits list of FITS transformed to astromrefs.

    If findnewobjects is True, this will run source extraction on each
    subtracted frame, remove all known sources from the photref, see if there
    are any new sources, add them to the source catalog as HAT-999 objects if
    there are no matches to them within catmatcharcsec arcseconds, regenerate
    the cmrawphot for the combined photref, and then run aperturephot on them.

    '''

    tasks = [(x, photreftype, outdir, kernelspec,
              reversesubtract, findnewobjects,
              photdisjointradius, refinfo, lcapertures)
             for x in xtrnsfits if os.path.exists(x)]

    print('%sZ: %s files to process' %
          (datetime.utcnow().isoformat(), len(tasks)))

    if len(tasks) > 0:

        pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

        # fire up the pool of workers
        results = pool.map(xtrnsfits_convsubphot_worker, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()

        return {x:y for (x,y) in results}

    else:

        print('ERR! %sZ: none of the files specified exist, bailing out...' %
              (datetime.utcnow().isoformat(),))
        return



#########################
## PHOTOMETRY DATABASE ##
#########################

def parse_iphot_line(iphotline):
    '''This parses the iphot line and returns a row of formatted elems.

    These can be then attached to other metadata and written directly to the
    database.

    '''

    photelemline = iphotline.rstrip(' \n')
    photelem = photelemline.split()

    if len(photelem) > 0:

        # get the station, framenumber, subframe, and ccd
        framekeyelems = FRAMESUBREGEX.findall(
            os.path.basename(iphot)
        )
        stf, cfn, cfs, fccd = (
            framekeyelems[0][0],
            framekeyelems[0][1],
            framekeyelems[0][2],
            framekeyelems[0][3]
        )
        net = 'HP' # this will be 'HP' for HATPI

        # these are things we get from the iphot line
        # hat, xcc, ycc, xic, yic
        # fsv, fdv, fkv, bgv, bge
        # ifl1, ife1, irm1, ire1, irq1
        # ifl2, ife2, irm2, ire2, irq2
        # ifl2, ife2, irm2, ire2, irq2

        objectid, xcc, ycc, xic, yic = photelem[:5]
        fsv, fdv, fkv, bgv, bge = photelem[5:10]
        ifl1, ife1, irm1, ire1, irq1 = photelem[10:15]
        ifl2, ife2, irm2, ire2, irq2 = photelem[15:20]
        ifl3, ife3, irm3, ire3, irq3 = photelem[20:]

        # format these as needed
        xcc = smartcast(xcc, float)
        ycc = smartcast(ycc, float)
        xic = smartcast(xic, float)
        yic = smartcast(yic, float)

        fsv = smartcast(fsv, float)
        fdv = smartcast(fdv, float)
        fkv = smartcast(fkv, float)
        bgv = smartcast(bgv, float)
        bge = smartcast(bge, float)

        ifl_000 = smartcast(ifl1, float)
        ife_000 = smartcast(ife1, float)
        irm_000 = smartcast(irm1, float)
        ire_000 = smartcast(ire1, float)
        irq_000 = smartcast(irq1, str)

        ifl_001 = smartcast(ifl2, float)
        ife_001 = smartcast(ife2, float)
        irm_001 = smartcast(irm2, float)
        ire_001 = smartcast(ire2, float)
        irq_001 = smartcast(irq2, str)

        ifl_002 = smartcast(ifl3, float)
        ife_002 = smartcast(ife3, float)
        irm_002 = smartcast(irm3, float)
        ire_002 = smartcast(ire3, float)
        irq_002 = smartcast(irq3, str)

        return {'objectid':objectid,
                'xcc':xcc,
                'ycc':ycc,
                'xic':xic,
                'yic':yic,
                'fsv':fsv,
                'fdv':fdv,
                'fkv':fkv,
                'bgv':bgv,
                'bge':bge,
                'ifl_000':ifl_000,
                'ife_000':ife_000,
                'irm_000':irm_000,
                'ire_000':ire_000,
                'irq_000':irq_000,
                'ifl_001':ifl_001,
                'ife_001':ife_001,
                'irm_001':irm_001,
                'ire_001':ire_001,
                'irq_001':irq_001,
                'ifl_002':ifl_002,
                'ife_002':ife_002,
                'irm_002':irm_002,
                'ire_002':ire_002,
                'irq_002':irq_002}

    else:

        return None




def convsub_photometry_to_ismphot_database(convsubfits,
                                           convsubphot,
                                           photreftype,
                                           subtracttype,
                                           kernelspec,
                                           lcapertures,
                                           projectid=None,
                                           field=None,
                                           ccd=None,
                                           overwrite=False,
                                           database=None):

    '''This inserts the ISM photometry from a single convsub FITS into the DB.

    If projectid, field, ccd are not provided, gets them from the FITS
    file. Also gets the photreftype from the filename of the
    convolved-subtracted photometry iphot file.

    '''

    # open a database connection
    if database:
        cursor = database.cursor()
        closedb = False
    else:
        database = pg.connect(user=PGUSER,
                              password=PGPASSWORD,
                              database=PGDATABASE,
                              host=PGHOST)
        cursor = database.cursor()
        closedb = True

    # start work here
    try:

        # figure out the projectid, field, ccd, photreftype
        # first, figure out the input frame's projid, field, and ccd
        frameelems = get_header_keyword_list(convsubfits,
                                             ['object',
                                              'projid'])
        felems = FRAMEREGEX.findall(
            os.path.basename(convsubfits)
        )

        if not (projectid and field and ccd):

            field, ccd, projectid = (frameelems['object'],
                                     int(felems[0][2]),
                                     frameelems['projid'])

        convsubdir = os.path.abspath(os.path.dirname(convsubfits))

        if not os.path.exists(iphotpath):
            print('ERR! %sZ: expected iphot %s for '
                  'convsub FITS %s does not exist, '
                  'not processing...' %
                  (datetime.utcnow().isoformat(), iphotpath, convsubfits))
            return (convsubfits, False)

        # find the frame's original FITS file (unsubtracted calibrated frame)
        originalfitsbasename = '%s-%s_%s.fits' % (felems[0][0],
                                                  felems[0][1],
                                                  felems[0][2])
        originalfitspath = os.path.join(convsubdir, originalfitsbasename)

        if not os.path.exists(originalfitspath):
            print('%ERR! sZ: expected original FITS %s '
                  'for convsub FITS %s does not exist, '
                  'not processing...' %
                  (datetime.utcnow().isoformat(),
                   originalfitspath, convsubfits))
            return (convsubfits, False)

        # figure out the frame's info from the original frame's header
        framerjd = get_header_keyword_list(originalfitspath,
                                           ['JD',''])

        # also get some metadata from the frameheader


        # now open the accompanying iphot file, and stream the photometry to the
        # database
        with open(convsubphot,'rb') as infd:

            # prepare the statement
            query = ("insert into ism_photometry ("
                     "frameutcdt, objectid, framekey, photkey, "
                     "xcc, ycc, xic, yic, bgv, bge, fsv, fdv, fkv, "
                     "ifl_000, ife_000, irm_000, ire_000, irq_000, "
                     "ifl_001, ife_001, irm_001, ire_001, irq_001, "
                     "ifl_002, ife_002, irm_002, ire_002, irq_002"
                     ") values ("
                     "%s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s, %s, %s"
                     ")")

            # prepare the input params



            for line in infd:

                parsedline = parse_iphot_line(line)



        # update the iphotfiles table file with all of this info. if there's a
        # uniqueness conflict, i.e. this same combination exists, then overwrite
        # if told to do so
        if overwrite:

            print('WRN! %sZ: overwriting existing photometry info in DB for %s'
                  %
                  (datetime.utcnow().isoformat(), convsubfits))

            query = ("insert into iphotfiles "
                     "(projectid, field, ccd, photreftype, convsubtype, "
                     "isactive, iphotfilepath, framerjd, framefilepath) "
                     "values ("
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s"
                     ") on conflict on constraint iphotfiles_pkey "
                     "do update "
                     "set projectid = %s, field = %s, ccd = %s, "
                     "photreftype = %s, convsubtype = %s, "
                     "isactive = %s, iphotfilepath = %s, framerjd = %s, "
                     "framefilepath = %s, entrytimestamp = current_timestamp")

            params = (projectid, field, ccd, photreftype, subtractiontype,
                      True, iphotpath, framerjd, originalfitspath,
                      projectid, field, ccd, photreftype, subtractiontype,
                      True, iphotpath, framerjd, originalfitspath)

        else:

            query = ("insert into iphotfiles "
                     "(projectid, field, ccd, photreftype, convsubtype, "
                     "isactive, iphotfilepath, framerjd, framefilepath) "
                     "values ("
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s"
                     ")")
            params = (projectid, field, ccd, photreftype, subtractiontype,
                      True, iphotpath, framerjd, originalfitspath)


        # execute the query to insert the object
        cursor.execute(query, params)
        database.commit()

        # update the iphotobjects table with all of these objects. if there's a
        # uniqueness conflict, i.e. this same combination exists, then overwrite
        # if told to do so

        if overwrite:

            query = ("insert into iphotobjects "
                     "(projectid, field, ccd, photreftype, convsubtype, "
                     "isactive, objectid, iphotfilepath, iphotfileline) "
                     "values ("
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s"
                     ") on conflict on constraint iphotobjects_pkey "
                     "do update set "
                     "projectid = %s, field = %s, ccd = %s, photreftype = %s, "
                     "convsubtype = %s, isactive = %s, objectid = %s, "
                     "iphotfilepath = %s, iphotfileline = %s, "
                     "entrytimestamp = current_timestamp")

        else:

            query = ("insert into iphotobjects "
                     "(projectid, field, ccd, photreftype, convsubtype, "
                     "isactive, objectid, iphotfilepath, iphotfileline) "
                     "values ("
                     "%s, %s, %s, %s, %s, "
                     "%s, %s, %s, %s"
                     ")")

        # execute statements for all of the iphot objects
        for ind, objectid in enumerate(iphotobjects):

            if overwrite:
                params = (projectid, field, ccd, photreftype, subtractiontype,
                          True, objectid, iphotpath, ind,
                          projectid, field, ccd, photreftype, subtractiontype,
                          True, objectid, iphotpath, ind,)
            else:
                params = (projectid, field, ccd, photreftype, subtractiontype,
                          True, objectid, iphotpath, ind)

            cursor.execute(query, params)

        database.commit()

        print('%sZ: convsub FITS %s with iphot %s and %s objects '
              'inserted into DB OK' %
              (datetime.utcnow().isoformat(),
               convsubfits,
               iphotpath,
               len(iphotobjects)) )

        # return True if everything succeeded
        returnval = (convsubfits, True)


    # catch the overwrite = False scenario
    except pg.IntegrityError as e:

        database.rollback()

        message = ('failed to insert photometry from %s '
                   'into DB because it exists already '
                   'and overwrite = False'
                   % convsubfits)
        print('EXC! %sZ: %s\n%s' %
               (datetime.utcnow().isoformat(), message, format_exc()) )
        returnval = (convsubfits, False)


    # if everything goes wrong, exit cleanly
    except Exception as e:

        database.rollback()

        message = 'failed to insert photometry from %s into DB' % convsubfits
        print('EXC! %sZ: %s\nexception was: %s' %
               (datetime.utcnow().isoformat(),
                message, format_exc()) )
        returnval = (convsubfits, False)
        raise


    finally:

        cursor.close()
        if closedb:
            database.close()

    return returnval



def parallel_convsubphotdb_worker(task):
    '''This wraps the function above for use with the parallel driver below.

    task[0] = convsubfits
    task[1] = {'projectid', 'field', 'ccd', 'overwrite'}

    '''

    convsubfits = task[0]
    kwargs = task[1]

    return convsub_photometry_to_ismphot_databse(convsubphots,**kwargs)



def parallel_convsubphot_to_db(convsubfitslist,
                               projectid=None,
                               field=None,
                               ccd=None,
                               overwrite=False,
                               nworkers=16,
                               maxworkertasks=1000):
    '''This runs a convsubphot ingest in parallel.

    '''

    tasks = [(x, {'projectid':projectid, 'field':field,
                  'ccd':ccd, 'overwrite':overwrite})
             for x in convsubfitslist if os.path.exists(x)]

    print('%sZ: %s files to process' %
          (datetime.utcnow().isoformat(), len(tasks)))

    if len(tasks) > 0:

        pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)


        # fire up the pool of workers
        results = pool.map(parallel_convsubphotdb_worker, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()

        return {x:y for (x,y) in results}

    else:

        print('ERR! %sZ: none of the files specified exist, bailing out...' %
              (datetime.utcnow().isoformat(),))
        return




############################
## LIGHT CURVE PRODUCTION ##
############################

# we'll make hatlc.sqlite type files, collecting them in /P/LC, under the
# following organization:
# {primary_field}/{hatid}-DR{datarelease}-V{lcversion}.hatlc.sqlite(.gz)
# we'll collect all photometry across observed fields and CCDs in the same file

def collect_lightcurve(objectid,
                       datarelease,
                       lcversion,
                       objectfield=None,
                       lcbasedir=LCBASEDIR,
                       updateifexists=True,
                       database=None):
    '''This collects lightcurves for objectid into a hatlc.sqlite file.

    We'll collect all photometry across observed fields and CCDs in the same
    file, organized as below:

    {LCBASEDIR} / DR{datarelease} / {primary_field} /
    {hatid}-DR{datarelease}-V{lcversion}.hatlc.sqlite

    If an LC doesn't exist for this objectid/datarelease/lcversion, then creates
    a new sqlitecurve. If one already exists and updateifexists is True, updates
    the light curve with new observations.

    FIXME: if we're going to run realtime imagesubphot, this will require that a
    collect_lightcurve is run after every image is taken. This will probably be
    stupidly slow...

    '''

    # open a database connection
    if database:
        database.autocommit = True # this is a readonly connection
        cursor = database.cursor()
        closedb = False
    else:
        database = pg.connect(user=PGUSER,
                              password=PGPASSWORD,
                              database=PGDATABASE,
                              host=PGHOST)
        database.autocommit = True
        cursor = database.cursor()
        closedb = True

    # start work here
    try:

        # find the photometry for the specified objectid
        query = ("select a.projectid, a.field, a.ccd, a.photreftype, "
                 "a.convsubtype, a.iphotfilepath, a.framefilepath, "
                 "a.framerjd, b.iphotfileline "
                 "from iphotfiles a join iphotobjects b on "
                 "(a.iphotfilepath = b.iphotfilepath) "
                 "where "
                 "(isactive = true) and "
                 "(b.objectid = %s) order by a.framerjd asc")
        params = (objectid, )

        cursor.execute(query, params)
        rows = cursor.fetchall()

        # if there are any rows, then there are some detections of this object
        if rows and len(rows) > 0:

            if HATIDREGEX.match(objectid):
                primary_field = objectid.split('-')[1]
            elif objectfield:
                primary_field = objectfield
            else:
                primary_field = 'nofield'

            # prepare the output directory
            lcdir = os.path.join(LCBASEDIR,
                                 'DR%s' % datarelease,
                                 '%s' % primary_field)

            # make the light curve directory if it doesn't exist
            if not os.path.exists(lcdir):
                os.mkdirs(lcdir)

            # this is the path to the light curve file itself
            lcfbasename = '{objectid}-DR{datarelease}-V{lcversion}-hatlc.sqlite'
            lcfpath = os.path.join(lcdir, lcfbasename)

            # if it doesn't exist, we can make a new file
            if not os.path.exists(lcfpath):

                sqlcdb = sqlite3.connect(lcfpath)
                sqlcur = sqlcdb.cursor()

                # read in the hatlc.sql file and execute it to generate the
                # tables for the hatlc
                hatlcsqlf = os.path.join(os.path.dirname(__file__),
                                         'hatlc.sql')
                with open(hatlcsqlf,'rb') as infd:
                    hatlcsql = infd.read()

                sqlcur.executescript(hatlcsql)
                # the lightcurve is now ready to use, but has no information in
                # it, we'll add this at a later step


            # check if the path to it exists and if updateifexisting is True
            # if so, we can update it
            if updateifexisting:

                print('WRN! %sZ: objectid %s has an '
                      'existing LC %s, updating...' %
                      (datetime.utcnow().isoformat(), objectid, lcfpath) )

            # if we're not allowed to update it, fail out
            else:

                print('ERR! %sZ: objectid %s has an existing LC %s '
                      'and updateifexisting = False, skipping this object...' %
                      (datetime.utcnow().isoformat(), objectid, lcfpath) )
                return objectid, None

            # now read in the iphot files one-by-one and write their photometry
            # lines to the sqlite
            for row in rows:

                (projectid, ofield, ccd, photreftype,
                 convsubtype, iphot, frame, rjd, iphotline) = row

                # open the iphot file
                with open(iphot,'rb') as infd:
                    for ind, line in enumerate(infd):
                        if ind == iphotline:
                            photelemline = line.rstrip(' \n')
                            break

                photelem = photelemline.split()

                if len(photelem) > 0:

                    # get the station, framenumber, subframe, and ccd
                    framekeyelems = FRAMESUBREGEX.findall(
                        os.path.basename(iphot)
                    )
                    stf, cfn, cfs, fccd = (
                        framekeyelems[0][0],
                        framekeyelems[0][1],
                        framekeyelems[0][2],
                        framekeyelems[0][3]
                    )
                    net = 'HP' # this will be 'HP' for HATPI

                    # these are things we get from the iphot line
                    # hat, xcc, ycc, xic, yic
                    # fsv, fdv, fkv, bgv, bge
                    # ifl1, ife1, irm1, ire1, irq1
                    # ifl2, ife2, irm2, ire2, irq2
                    # ifl2, ife2, irm2, ire2, irq2

                    hat, xcc, ycc, xic, yic = photelem[:5]
                    fsv, fdv, fkv, bgv, bge = photelem[5:10]
                    ifl1, ife1, irm1, ire1, irq1 = photelem[10:15]
                    ifl2, ife2, irm2, ire2, irq2 = photelem[15:20]
                    ifl3, ife3, irm3, ire3, irq3 = photelem[20:]

                    # format these as needed
                    xcc = smartcast(xcc, float)
                    ycc = smartcast(ycc, float)
                    xic = smartcast(xic, float)
                    yic = smartcast(yic, float)

                    fsv = smartcast(fsv)
                    fdv = smartcast(fdv)
                    fkv = smartcast(fkv)
                    bgv = smartcast(bgv)
                    bge = smartcast(bge)




    # if everything goes wrong, exit cleanly
    except Exception as e:

        database.rollback()

        message = 'failed to collect light curve for %s, DR%s, V%s' % (
            objectid, datarelease, lcversion
        )
        print('EXC! %sZ: %s\nexception was: %s' %
               (datetime.utcnow().isoformat(),
                message, format_exc()) )
        returnval = (objectid, None)
        raise


    finally:

        cursor.close()
        if closedb:
            database.close()

    return returnval


#############################
## LIGHT CURVE EPD AND TFA ##
#############################



###########################
## CRON ROLLUP FUNCTIONS ##
###########################

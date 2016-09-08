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
from hashlib import md5
import gzip

import numpy as np


import aperturephot as ap
import imagesubphot as ism
from imageutils import get_header_keyword_list

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

# this defines the field string and CCDs
FIELD_REGEX = re.compile('^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$')
FIELD_CCDS = [5,6,7,8]

# defines where the reference frames go
REFBASEDIR = '/P/HP0/BASE/reference-frames'
REFINFO = os.path.join(REFBASEDIR,'TM-refinfo.sqlite')

# define where the frameinfo cache is
FRAMEINFOCACHEDIR = '/P/HP0/BASE/frameinfo-cache'

###############
## UTILITIES ##
###############

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
    '''This searches in dirlist for all astromref-shifted FITS files matching the
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


##################################
## ASTROMETRIC REFERENCE FRAMES ##
##################################

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
            areftargetjpeg = areftargetfits.replace('.fits','.jpeg')
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

            db = sqlite3.connect(refinfo)
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

    db = sqlite3.connect(REFINFO)
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

    '''

    try:

        frame, outdir, refinfo = task

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

                    print('%sZ: SHIFT OK %s -> %s' %
                          (datetime.utcnow().isoformat(), frame, shifted_frame))

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
                             nworkers=16,
                             maxworkertasks=1000):
    '''This calculates the shifts between frames in fitsfiles and the appropriate
    astromref for the projectid, field and CCD, then shifts each frame to the
    astromref's coordinate system, generating -xtrns.fits files.

    '''

    print('%sZ: %s files to process' %
          (datetime.utcnow().isoformat(), len(fitsfiles)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, outdir, refinfo) for x in fitsfiles if os.path.exists(x)]

    # fire up the pool of workers
    results = pool.map(frames_astromref_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {x:y for (x,y) in results}


##################################
## PHOTOMETRIC REFERENCE FRAMES ##
##################################

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

            print('ERR! %sZ: frame %s has missing fiphot/fistar, skipping...' %
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



def generate_photref_candidates_from_xtrns(fitsfiles,
                                           makeactive=True,
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

    # then, apply our conditions to these fits files to generate a list of
    # photref candidates





def generate_combined_photref(photreftarget,
                              photrefcandidates,
                              makeactive=True,
                              field=None,
                              ccd=None,
                              projectid=None,
                              refdir=REFBASEDIR,
                              refinfo=REFINFO,
                              combinetype='median',
                              kernelspec='b/4;i/4;d=4/4',
                              ccdgain=None,
                              zeropoint=None,
                              ccdexptime=None,
                              extractsources=True,
                              apertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                              framewidth=None,
                              searchradius=8.0,
                              nworkers=16,
                              maxworkertasks=1000):
    '''This generates a combined photref from photref target and candidates and
    updates the TM-refinfo.sqlite database.

    '''



def get_combined_photref(projectid, field, ccd, refinfo=REFINFO):
    '''
    This gets the combined photref for the given combo of projid, field, ccd.

    '''



##################################
## IMAGE SUBTRACTION PHOTOMETRY ##
##################################

def find_new_sources(subtractedframe,
                     fluxthreshold=1000,
                     catmatcharcsec=3.0,
                     refinfo=REFINFO):
    '''
    This finds new sources in the subtracted frame.

    '''


def xtrsfits_convsubphot_worker(task):
    '''
    This is a parallel worker for framelist_convsubphot_photref below.

    task[0] = xtrnsfits file
    task[1] = outdir
    task[2] = kernelspec
    task[3] = reversesubtract boolean
    task[4] = findnewobjects boolean

    '''

    frame, outdir, kernelspec, reversesubtract, findnewobjects = task

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

        # then, find the associated combined photref frame

        # then, find the associated combined photref registration file

        # then, find the associated combined photref's cmrawphot for photometry

        # generate the output subtracted frame filename and kernel filename

        # do the subtraction (take care of reversesubtract here)

        # find new objects in the subtracted frame if told to do so

        # then do photometry on the subtracted frame


    except Exception as e:

        print('ERR! %sZ: could not do convsubphot on frame %s, error was: %s' %
              (datetime.utcnow().isoformat(), frame, e))
        return frame, None




def xtrnsfits_convsubphot(xtrnsfits,
                          outdir=None,
                          refinfo=REFINFO,
                          kernelspec='b/4;i/4;d=4/4',
                          findnewobjects=True,
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


#########################
## PHOTOMETRY DATABASE ##
#########################



############################
## LIGHT CURVE PRODUCTION ##
############################



#############################
## LIGHT CURVE EPD AND TFA ##
#############################



###########################
## CRON ROLLUP FUNCTIONS ##
###########################

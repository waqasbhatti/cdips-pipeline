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


import aperturephot as ap
import imagesubphot as ism
from imageutils import get_header_keyword_list


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

###############
## UTILITIES ##
###############

def find_original_fits_fieldprojectidccd(dirlist,
                                         field,
                                         projectid,
                                         ccd):
    '''This searches in dirlist for all original FITS files matching the specified
    projectid, field, and ccd combination.

    Returns a flat list of FITS.

    '''



def find_arefshifted_fits_fieldprojectidccd(dirlist,
                                            field,
                                            projectid,
                                            ccd):
    '''This searches in dirlist for all astromref-shifted FITS files matching the
    specified projectid, field, and ccd combination.

    Returns a flat list of FITS.

    '''


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



def framelist_to_astromref(fitsfiles,
                           outdir=None,
                           refinfo=REFINFO,
                           nworkers=16,
                           maxworkertasks=1000):
    '''This calculates the shifts between frames in fitsfiles and the appropriate
    astromref for the projectid, field and CCD, then shifts each frame to the
    astromref's coordinate system.

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

def generate_photref_candidates(fitsfiles,
                                makeactive=True,
                                minframes=50,
                                maxhourangle=3.0,
                                maxmoonphase=25.0,
                                maxmoonelev=0.0,
                                maxzenithdist=30.0,
                                maxbackgroundstdev=10.0,
                                maxbackgroundmedian=1000.0,
                                forcecollectinfo=False):
    '''This uses ism.select_photref_frames run on fitsfiles to get photref
    candidates.

    fitsfiles must be a list of frames, which have been already transformed to
    the astromref, and are all from a single projectid, ccd, field combination
    for this operation to make sense.

    '''

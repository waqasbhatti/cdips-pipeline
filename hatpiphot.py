#!/usr/bin/env python

'''
hatpiphot.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Dec 2014

Contains aperture photometry routines for HATPI.

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
import ConfigParser
import json

import numpy as np
import pyfits

from scipy.spatial import cKDTree as kdtree

import imageutils

# get fiphot binary reader
try:
    from HATpipepy.Common.BinPhot import read_fiphot
    HAVEBINPHOT = True
except:
    print("can't import binary fiphot reading functions from "
          "HATpipe, binary fiphot files will be unreadable!")
    HAVEBINPHOT = False

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
             'path':'/S/CAT/2MASS/2MASS_JH_AP/data'},
    'UCAC4':{'cmd':UCAC4READCMD,
             'path':'/S/CAT/UCAC4'}
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
FIPHOTCMD = ("fiphot --input {fits} --input-list {sourcelist} "
             "--col-id 1 --col-xy {xycols} --gain {ccdgain:f} "
             "--mag-flux {zeropoint:f},{ccdexptime:f} "
             "--apertures {aperturelist} "
             "--sky-fit 'mode,sigma=3,iterations=2' --disjoint-radius 2 "
             "--serial {fitsbase} "
             "--format 'ISXY,BbMms' --nan-string 'NaN' "
             "--aperture-mask-ignore 'saturated' --comment '--comment' "
             "--single-background 3 {binaryout} --output {outfile} -k")


# MagnitudeFitting.py commandline
# {magfitexec}:         path to the executable to run for magfit (usually
#                       MagnitudeFitting.py)
# {network}:            HATNet or HATSouth
# {fit_type}:           single for single reference mag fit, master for
#                       master mag fit
# {sphotref_frame}:     FITS to use as the single photometric reference
# {sphotref_phot} :     fiphot file for the single photometric reference
# {nprocs}        :     number of processors to use
# {magfit_config_file}: path to the magfit config file
# {magfit_frame_list}:  path to the magfit frame list produced by
#                       get_magfit_frames
MAGFITCMD = ("python {magfitexec} {network} {fit_type} "
             "{singlephotref_frame} {singlephotref_phot} "
             "-p {nprocs} --config-file={magfit_config_file} "
             "--manual-frame-list={magfit_frame_list} --stat")

# command to generate master photometric reference
# {mphotrefexec}:        path to executable to run (usually do_masterphotref.py)
# {network}:             HATNet or HATSouth
# {sphotref_frame}:      FITS used for the single photometric reference
# {fiphot_list}:         list of fiphot files containing all photometry
# {magfit_config_file}:  path to the magfit config file
MPHOTREFCMD = ("python {mphotrefexec} {network} {sphotref_frame} "
               "--manual-frame-list={fiphot_list} "
               "--config-file={magfit_config_file} --nostat")



##########################
## PHOTOMETRY FUNCTIONS ##
##########################

def make_fov_catalog(ra=None, dec=None, size=None,
                     brightrmag=6.0,
                     faintrmag=13.0,
                     fits=None,
                     outfile=None,
                     catalog='2MASS',
                     columns=None):
    '''
    This function gets all the sources in the field of view of the frame, given
    its central pointing coordinates and plate-scale from either 2MASS or
    UCAC4. Makes a catalog file that can then be used as input to
    make_source_list below.

    catalog = 'UCAC4' or '2MASS'

    if ra, dec, size are None, fits must not be None. fits is the filename of
    the FITS file to get the center RA, DEC, and platescale values from.

    Returns the path of the catalog file produced.

    '''

    if ra and dec and size:

        catra, catdec, catbox = ra, dec, size

    elif fits:
        raise NotImplementedError('not done yet!')

    else:
        print('%sZ: need a FITS file to work on, or center coords and size' %
              (datetime.utcnow().isoformat(),))
        return


    if not outfile:
        outfile = '%s-RA%s-DEC%s-SIZE%s.catalog' % (catalog,
                                                    catra,
                                                    catdec,
                                                    catsize)

    print('%sZ: making FOV catalog for '
          'center RA, DEC = %.5f, %.5f with size = %.5f deg' %
          (datetime.utcnow().isoformat(),
           catra, catdec, catbox))


    catalogcmd = CATALOGS[catalog]['cmd'].format(
        ra=catra,
        dec=catdec,
        boxlen=catbox,
        catalogpath=CATALOGS[catalog]['path'],
        brightrmag=brightrmag,
        faintrmag=faintrmag,
        outfile=outfile
        )

    if DEBUG:
        print(catalogcmd)

    # execute the cataloger shell command
    catalogproc = subprocess.Popen(shlex.split(catalogcmd),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

    # get results
    catalog_stdout, catalog_stderr = catalogproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if catalogproc.returncode == 0:

        print('%sZ: FOV catalog %s generated for '
              'center RA, DEC = %.5f, %.5f with size = %.5f deg' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(outfile), catra, catdec, catbox))

        return os.path.abspath(outfile)

    else:

        print('%sZ: FOV catalog generation failed for '
              'center RA, DEC = %.5f, %.5f with size = %.5f deg!' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(outfile), catra, catdec, catbox))

        return None



def reform_fov_catalog(incat,
                       outcat,
                       columns='id,ra,dec,xi,eta,J,K,qlt,I,r'):
    '''
    This converts the full output catalog from 2massread, etc. to the format
    required for magfit. Also useful for general reforming of the columns.

    columns is a CSV string containing columns needed from allcolumns below.

    '''

    allcolumns = ['id','ra','dec','xi','eta','arcdis','J',
                  'Junc','H','Hunc','K','Kunc','qlt','B',
                  'V','R','I','u','g','r','i','z','field','num']


    columns = columns.split(',')
    colstoget = [allcolumns.index(x) for x in columns]

    inf = open(incat,'rb')
    outf = open(outcat,'wb')

    for line in inf:
        if '#' not in line:
            sline = line.split()
            outcols = [sline[x] for x in colstoget]
            outline = ' '.join(outcols)
            outf.write('%s\n' % outline)



def extract_frame_sources(fits,
                          outfile,
                          fistarexec='fistar',
                          ccdextent='0:0,2048:2048',
                          ccdgain=2.725,
                          fluxthreshold=1000,
                          zeropoint=17.11,
                          exptime=30.0):
    '''
    This uses fistar to extract sources from the image.

    fistar -i 1-377741e_5.fits -o test.fistar --model elliptic --iterations symmetric=4,general=2 --algorithm uplink --format id,x,y,bg,amp,s,d,k,flux,s/n -g 2.725 --mag-flux 17,30 --sort flux --flux-threshold 1000 --section 0:0,2048:2048
    '''

    if not os.path.exists(fits):

        print('%sZ: no FITS file to work on!' %
              (datetime.utcnow().isoformat(),))
        return None

    if not outfile:
        outfile = fits.strip('.fits.fz') + '.fistar'

    fistarcmd = FISTARCMD.format(
        fistarexec=fistarexec, # assuming fistar is in the path
        frame=fits,
        extractedlist=outfile,
        ccdgain=ccdgain,
        zeropoint=zeropoint,
        exptime=exptime,
        fluxthreshold=fluxthreshold,
        ccdsection=ccdextent
        )

    if DEBUG:
        print(fistarcmd)

    print('%sZ: starting fistar for %s...' %
          (datetime.utcnow().isoformat(), fits))

    # execute the fistar shell command
    fistarproc = subprocess.Popen(shlex.split(fistarcmd),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    # get results
    fistar_stdout, fistar_stderr = fistarproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if fistarproc.returncode == 0:

        print('%sZ: fistar completed for %s: %s' %
              (datetime.utcnow().isoformat(),fits, outfile))

        return outfile

    else:

        print('%sZ: fistar failed for %s!' %
              (datetime.utcnow().isoformat(), fits))

        return None


def parallel_sourceextract_worker(task):
    '''
    This expands the task arg into the args and kwargs necessary for
    extract_frame_sources.

    '''

    return (task[0][0], extract_frame_sources(*task[0],**task[1]))



def parallel_extract_sources(fitsdir,
                             outdir,
                             nworkers=8,
                             maxtasksperworker=1000,
                             fistarexec='fistar',
                             ccdextent='0:0,2048:2048',
                             ccdgain=2.725,
                             fluxthreshold=1000,
                             zeropoint=17.11,
                             exptime=30.0):
    '''
    This does parallel source extraction from all FITS in fitsdir, and puts the
    results in outdir.

    '''

    # get a list of all fits files in the directory
    fitslist = glob.glob(os.path.join(fitsdir,'*.fits'))

    print('%sZ: found %s FITS files in %s, starting source extraction...' %
          (datetime.utcnow().isoformat(),
           len(fitslist), fitsdir))

    if outdir and not os.path.exists(outdir):

        print('%sZ: making new output directory %s' %
              (datetime.utcnow().isoformat(),
               outdir))
        os.mkdir(outdir)

    pool = mp.Pool(nworkers, maxtasksperchild=maxtasksperworker)

    tasks = [
        [(x, os.path.join(outdir,
                          os.path.basename(x.strip('.fits.fz') +
                                           '.fistar'))),
         {'fistarexec':fistarexec,
          'ccdextent':ccdextent,
          'ccdgain':ccdgain,
          'fluxthreshold':fluxthreshold,
          'zeropoint':zeropoint,
          'exptime':exptime,}]
        for x in fitslist
        ]

    # fire up the pool of workers
    results = pool.map(parallel_sourceextract_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    # this is the return dictionary
    returndict = {x:y for (x,y) in results}
    return returndict



def match_fovcatalog_framesources(frame_extracted_sourcelist,
                                  frame_projected_fovcatalog,
                                  outfile,
                                  srclist_cols=(0,1,2,5,6,7),
                                  fovcat_cols=(0,1,2,12,13),
                                  match_pixel_distance=1.0):
    '''Does frame_projected_fovcatalog and frame_extracted_sourcelist matching.

    This matches the fovcatalog transformed to pixel coordinates to the
    extracted source list pixel coordinates and gets the IDs of the sources to
    be later used when creating the fiphot file.

    The procedure is:

    1. read frame sourcelist xy cols into ndarray
    2. read projected fovcatalog xy cols into ndarray
    3. vstack both ndarrays, combined = [sourcelist, fovcatalog]
    4. create a kd-Tree using tree = cKDTree(combined)
    5. find pairs using tree.query_pairs(match_pixel_distance, p=2) [p=2 ->
       standard euclidean distance for the Minkowski norm]
    6. convert result to list of index pairs
    7. get indices of fovcatalog for each pair and use to get ID from fovcatalog
    8. get indices of framelist for each pair and use to get S, D, K info
    9. make sure fovcatalog IDs are unique
    10. write to a combined fiphot file.

    '''

    srclist = np.genfromtxt(frame_extracted_sourcelist,
                           dtype='S17,f8,f8,f8,f8,f8',
                           names=['id','x','y','s','d','k'],
                           usecols=srclist_cols)

    fovcat = np.genfromtxt(frame_projected_fovcatalog,
                            dtype='S17,f8,f8,f8,f8',
                            names=['id','ra','dec','x','y'],
                            usecols=fovcat_cols)

    srclist_xys = np.column_stack((srclist['x'],srclist['y']))
    fovcat_xys = np.column_stack((fovcat['x'],fovcat['y']))

    fovcat_tree = kdtree(fovcat_xys)

    distances, fovcat_indices = fovcat_tree.query(srclist_xys)

    # now that we have matches, put the two files together based on the indices
    # the output format will be:
    # srcext = source extracted frame list from fistar
    # fovcat = frame projected FOV catalog object list from
    # make_frame_sourcelist
    # fovcat HAT-ID, fovcat RA, fovcat DEC, fovcat x, fovcat y,
    # srcext id, srcext x, srcext y, srcext S, srcext D, srcext K
    if not outfile:
        outfile = (frame_extracted_sourcelist.strip('.fits.fz') +
                   '.matched-sources')


    outf = open(outfile,'wb')

    outlinestr = '%s %.5f %.5f %.5f %.5f %s %.5f %.5f %.5f %.5f %.5f\n'

    for dist, fovcat_ind, srclist_ind in zip(distances,
                                             fovcat_indices,
                                             range(len(srclist))):

        if dist < match_pixel_distance:
            outf.write(outlinestr %
                       (fovcat['id'][fovcat_ind],
                        fovcat['ra'][fovcat_ind],
                        fovcat['dec'][fovcat_ind],
                        fovcat['x'][fovcat_ind],
                        fovcat['y'][fovcat_ind],
                        srclist['id'][srclist_ind],
                        srclist['x'][srclist_ind],
                        srclist['y'][srclist_ind],
                        srclist['s'][srclist_ind],
                        srclist['d'][srclist_ind],
                        srclist['k'][srclist_ind]))

    outf.close()

    return outfile



def make_frameprojected_catalog(fits,
                                catalog,
                                catalogxycols=(12,13),
                                transformer='anrd2xy',
                                wcs=None,
                                ccdextent=None,
                                out=None,
                                removetemp=True,
                                pixborders=0.0):

    '''
    This makes the projected catalog for the frame to use with fiphot using the
    anet rdtoxy transform stored in the filename wcs and projects the catalog
    into pixel coordinates of the image. If wcs is None, a file with .wcs
    extension in the same directory as fits will be used as input. catalog is
    the path to the FOV catalog generated by 2massread or ucac4read. out is the
    name of the output sourcelist. Uses the executable defined in the
    transformer kwarg (usually a variant of the rd2xy binary from
    astrometry.net). ccdextent is a dictionary like CCDEXTENT above noting the
    extent of the CCD and tells us which objects are outside the FOV so they're
    removed from the output source list.

    Returns the path of the source list file produced if successful, otherwise
    returns None.

    '''

    fitspath = os.path.abspath(fits)

    if wcs:
        framewcsfile = wcs
    else:
        wcspath = fitspath.rstrip('.fits.fz')
        wcspath = wcspath + '.wcs'
        framewcsfile = wcspath

    if out:
        outfile = out
        temppath = out + '.projcattemp'
    else:
        outpath = fitspath.rstrip('.fits.fz')
        temppath = outpath + '.projcattemp'
        outpath = outpath + '.projcatalog'
        outfile = outpath

    # format the transformer shell command
    transformcmd = TRANSFORMCMD.format(transformer=transformer,
                                       framewcsfile=framewcsfile,
                                       catalogsourcelist=catalog,
                                       outputfile=temppath)

    if DEBUG:
        print(transformcmd)

    # execute the transformer shell command
    transformproc = subprocess.Popen(shlex.split(transformcmd),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

    # get results
    transform_stdout, transform_stderr = transformproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if transformproc.returncode == 0:

        # now we need to take out sources outside the CCD extent
        sourcelist_x, sourcelist_y = np.loadtxt(temppath,
                                                usecols=catalogxycols,
                                                unpack=True)

        # get the extent of the CCD
        if not ccdextent:
            ccdextent = CCDEXTENT

        # get indices for the lines to be kept
        # optionally, remove all sources within pixborders pixels of the edges
        # of the image.
        keep_ind = np.where(
            (sourcelist_x > (ccdextent['x'][0] + pixborders)) &
            (sourcelist_x < (ccdextent['x'][1] - pixborders)) &
            (sourcelist_y > (ccdextent['y'][0] + pixborders)) &
            (sourcelist_y < (ccdextent['y'][1] - pixborders))
        )[0].tolist()

        # output the lines to be kept
        outf = open(outfile, 'wb')
        with open(temppath,'rb') as tempf:

            templines = tempf.readlines()
            templines = [x for x in templines if '#' not in x]
            for ind in keep_ind:
                outf.write(templines[ind])

        outf.close()

        if removetemp:
            os.remove(temppath)

        print('%sZ: frame source list generation OK for %s' %
              (datetime.utcnow().isoformat(),
               fits))

        return outfile
    else:
        print('%sZ: frame source list generation '
              'failed for %s: error was %s' %
              (datetime.utcnow().isoformat(),
               fits,
               transform_stderr))
        if removetemp:
            try:
                os.remove(temppath)
            except:
                pass
        return None



def run_fiphot(fits,
               sourcelist=None,
               xycols='25,26', # set for full 2MASS fov catalog output format
               ccdgain=None,
               zeropoint=None,
               ccdexptime=None,
               aperturelist='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
               outfile=None,
               removesourcelist=False,
               binaryoutput=True):
    '''
    Thus runs fiphot for a single frame. Only the fits filename is required. If
    other parameters are not provided, they will be obtained from the image
    header and defaults.

    Returns the path of the .fiphot file produced if successful, otherwise
    returns None.

    '''

    # get the required header keywords from the FITS file
    header = imageutils.get_header_keyword_list(fits,
                                                ['GAIN',
                                                 'GAIN1',
                                                 'GAIN2',
                                                 'EXPTIME'])

    # handle the gain and exptime parameters
    if not ccdgain:

        # FIXME: is this right? should we use separate gain values for each side
        # of the CCD? what about stars that straddle the middle?
        if 'GAIN1' in header and 'GAIN2' in header:
            ccdgain = (header['GAIN1'] + header['GAIN2'])/2.0
        elif 'GAIN' in header:
            ccdgain = header['GAIN']
        else:
            ccdgain = None

    if not ccdexptime:
        ccdexptime = header['EXPTIME'] if 'EXPTIME' in header else None

    if not (ccdgain or ccdexptime):
        print('%sZ: no GAIN or EXPTIME defined for %s' %
              (datetime.utcnow().isoformat(),
               fits))
        return None

    # figure out the fitsbase from the fits filename
    fitsbase = os.path.basename(fits).strip('.fits.fz')

    # handle the zeropoints
    if not zeropoint:

        # if the zeropoint isn't provided and if this is a HAT frame, the ccd
        # number will get us the zeropoint in the ZEROPOINTS dictionary
        frameinfo = FRAMEREGEX.findall(fits)
        if frameinfo:
            zeropoint = ZEROPOINTS[int(frameinfo[0][-1])]
        else:
            print('%sZ: no zeropoint magnitude defined for %s' %
                  (datetime.utcnow().isoformat(),
                   fits))
            return None


    # figure out the output path
    if not outfile:
        outfile = os.path.abspath(fits.strip('.fits.fz') + '.fiphot')

    # figure out the sourcelist path
    if not sourcelist:
        sourcelist = os.path.abspath(fits.strip('.fits.fz') + '.sourcelist')
        if not os.path.exists(sourcelist):

            print("%sZ: can't find a source list for %s" %
                  (datetime.utcnow().isoformat(),
                   fits))
            return None

    # figure out if we want binary output or not
    if binaryoutput:
        binaryout = '--binary-output'
    else:
        binaryout = ''


    # assemble the fiphot command string
    fiphotcmd = FIPHOTCMD.format(fits=fits,
                                 sourcelist=sourcelist,
                                 xycols=xycols,
                                 zeropoint=zeropoint,
                                 ccdgain=ccdgain,
                                 ccdexptime=ccdexptime,
                                 aperturelist=aperturelist,
                                 fitsbase=fitsbase,
                                 binaryout=binaryout,
                                 outfile=outfile)

    if DEBUG:
        print(fiphotcmd)

    # execute the fiphot command
    fiphotproc = subprocess.Popen(shlex.split(fiphotcmd),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    # get results
    fiphot_stdout, fiphot_stderr = fiphotproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if fiphotproc.returncode == 0:

        print('%sZ: photometry using fiphot OK for %s' %
              (datetime.utcnow().isoformat(),
               fits))

        if removesourcelist:
            os.remove(sourcelist)

        return outfile


    else:
        print('%sZ: fiphot failed for %s: error was %s' %
              (datetime.utcnow().isoformat(),
               fits,
               fitphot_stderr))
        return None



def do_photometry(fits,
                  fovcatalog,
                  extractsources=True,
                  fovcat_xycols=(12,13),
                  fiphot_xycols='7,8', # set for matched source list
                  outdir=None,
                  ccdextent=None,
                  removesourcetemp=True,
                  pixborders=0.0,
                  aperturelist='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                  removesourcelist=False,
                  binaryoutput=True):
    '''This rolls up the sourcelist and fiphot functions above.

    Runs both stages on fits, and puts the output in outdir if it exists. If it
    doesn't or is None, then puts the output in the same directory as fits.

    fovcatalog is the path to the 2MASS/UCAC4 catalog for all sources in the
    observed field.

    '''

    outprojcat = os.path.basename(fits).strip('.fits.fz') + '.projcatalog'
    outsourcelist = os.path.basename(fits).strip('.fits.fz') + '.sourcelist'
    outfiphot = os.path.basename(fits).strip('.fits.fz') + '.fiphot'

    if outdir and os.path.exists(outdir):

        outprojcat = os.path.join(outdir, outprojcat)
        outsourcelist = os.path.join(outdir,outsourcelist)
        outfiphot = os.path.join(outdir, outfiphot)

    elif outdir and not os.path.exists(outdir):

        print('%sZ: making new output directory %s' %
              (datetime.utcnow().isoformat(),
               outdir))
        os.mkdir(outdir)
        outprojcat = os.path.join(outdir, outprojcat)
        outsourcelist = os.path.join(outdir,outsourcelist)
        outfiphot = os.path.join(outdir, outfiphot)

    else:

        outprojcat, outsourcelist, outfiphot = None, None, None

    if DEBUG:
        print('output projected catalog will be %s' % outprojcat)
        print('output source list will be %s' % outsourcelist)
        print('output fiphot will be %s' % outfiphot)

    # get the frame project catalog
    projcatfile = make_frameprojected_catalog(fits,
                                              fovcatalog,
                                              ccdextent=ccdextent,
                                              out=outprojcat,
                                              removetemp=removesourcetemp,
                                              pixborders=pixborders)

    if projcatfile:

        # if we're supposed to extract sources and run photometry on them
        # instead of just the sources in the projected fovcatalog, do so
        if extractsources:

            # extract sources
            framesources = extract_frame_sources(
                fits,
                os.path.join(
                    outdir,
                    os.path.basename(fits).strip('.fits.fz') + '.fistar'
                    )
                )

            if framesources:

                # match these to the projected fovcatalog
                matchedsources = match_fovcatalog_framesources(
                    framesources,
                    projcatfile,
                    outsourcelist
                    )

                fiphot_xycols = '7,8'

            else:

                print('%sZ: extracting sources failed for %s!' %
                      (datetime.utcnow().isoformat(), fits))
                return None


        # run fiphot on the source list
        fiphotfile = run_fiphot(fits,
                                sourcelist=outsourcelist,
                                aperturelist=aperturelist,
                                outfile=outfiphot,
                                xycols=fiphot_xycols,
                                removesourcelist=removesourcelist,
                                binaryoutput=binaryoutput)

        if fiphotfile:

            return fiphotfile

        else:

            print('%sZ: photometry failed for %s!' %
                  (datetime.utcnow().isoformat(), fits))

    else:

        print('%sZ: creating a projected source catalog failed for %s!' %
              (datetime.utcnow().isoformat(), fits))
        return None



def fitsdir_photometry(fitsdir,
                       outdir,
                       fovcatalog,
                       ccdextent=None,
                       pixborders=0.0,
                       aperturelist='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                       removesourcetemp=True,
                       removesourcelist=False,
                       binaryoutput=True):
    '''This does photometry for all FITS file in a directory.

    Puts the results in outdir. Returns a dictionary with photometry output for
    each file.

    '''

    # get a list of all fits files in the directory
    fitslist = glob.glob(os.path.join(fitsdir,'*.fits'))
    resultdict = {}

    print('%sZ: found %s FITS files in %s, starting photometry...' %
          (datetime.utcnow().isoformat(),
           len(fitslist), fitsdir))

    for i, fits in enumerate(fitslist):

        print('%sZ: working on %s (%s/%s)' %
              (datetime.utcnow().isoformat(),
               fits, i+1, len(fitslist)))

        fiphotresult = do_photometry(fits,
                                     fovcatalog,
                                     outdir=outdir,
                                     ccdextent=ccdextent,
                                     pixborders=pixborders,
                                     aperturelist=aperturelist,
                                     removesourcetemp=removesourcetemp,
                                     removesourcelist=removesourcelist,
                                     binaryoutput=binaryoutput)
        resultdict[i] = [fits, fiphotresult]

    return resultdict



def parallel_photometry_worker(task):
    '''
    This is the parallel photometry worker function for use with
    parallel_fitsdir_photometry below. Just calls do_photometry with expanded
    args and kwargs from the two element task list. task[0] is a tuple of args,
    and task[1] is a dictionary of kwargs.

    Returns a tuple of form: (fits, fiphot)

    '''

    # task[0] is args and first arg is the path to the fits file
    return (task[0][0], do_photometry(*task[0], **task[1]))



def parallel_fitsdir_photometry(
        fitsdir,
        outdir,
        fovcatalog,
        ccdextent=None,
        pixborders=0.0,
        aperturelist='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
        removesourcetemp=True,
        removesourcelist=False,
        binaryoutput=True,
        nworkers=16,
        maxtasksperworker=1000,
        resultstojson=None
        ):
    '''
    This does photometry for all FITS files in a directory using nworkers
    parallel workers.

    '''

    # get a list of all fits files in the directory
    fitslist = glob.glob(os.path.join(fitsdir,'*.fits'))

    print('%sZ: found %s FITS files in %s, starting photometry...' %
          (datetime.utcnow().isoformat(),
           len(fitslist), fitsdir))

    if outdir and not os.path.exists(outdir):

        print('%sZ: making new output directory %s' %
              (datetime.utcnow().isoformat(),
               outdir))
        os.mkdir(outdir)

    pool = mp.Pool(nworkers,maxtasksperchild=maxtasksperworker)

    tasks = [[(x, fovcatalog),
              {'outdir':outdir,
               'ccdextent':ccdextent,
               'pixborders':pixborders,
               'aperturelist':aperturelist,
               'removesourcetemp':removesourcetemp,
               'removesourcelist':removesourcelist,
               'binaryoutput':binaryoutput}] for x in fitslist]

    # fire up the pool of workers
    results = pool.map(parallel_photometry_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    # this is the return dictionary
    returndict = {x:y for (x,y) in results}

    if resultstojson:
        resultsfile = open(os.path.join(outdir,
                                        '%s-photometry-results.json' %
                                        os.path.basename(outdir)),'wb')
        resultsfile.write(json.dumps(returndict,ensure_ascii=True))
        resultsfile.close()
    else:
        return returndict



# remaining steps:
# 0) find a frame to use as the single photometric reference
# 1) run MagnitudeFitting.py with the single reference mode
# 2) using do_masterphotref.py to generate the master frames
# 3) run MagnitudeFitting.py with the master reference mode
# 4) run run_fitpsf.py to generate .psftrans files
# 5) using lc_gen.sh to generate the .rlc files
# 6) using the do_epd.py to generate .epd files
# 7) run tfa on those. this can be done by calling tfa from command line.


#################################
## MAGNITUDE FITTING FUNCTIONS ##
#################################

def get_magfit_frames(fitsdir,
                      fitsglob,
                      photdir,
                      workdir=None,
                      headerfilters=None,
                      selectreference=True,
                      framestats=False,
                      linkfiles=True,
                      outlistfile=None):
    '''
    fitsdir = directory where the FITS object frames are
    fitsglob = glob to select a subset of the FITS object frames
    photdir = directory where the TEXT fiphot files are (to select a
              singlephotref)
    workdir = directory where to create FITS and fiphot symlinks so
    MagnitudeFitting.py can work there

    This does the following:

    1. gets a list of frames in fitsdir that match fitsglob
    2. associates these files with their photometry files in photdir
    3. optionally filters them by the definitions in headerfilters (TBD)
    4. optionally creates symlinks to the chosen frame files in the workdir
       (photdir by default, directory will be created if non-existent)
    5. optionally selects a frame that might be a good reference for single
       reference magfit
    6. optionally writes the frame list to outlistfile.

    headerfilters is a list of string elements of the form:

    '<FITS header key> <operator> <FITS header value>'

    where <operator> is a standard binary operator. This string will be evaled
    to get the filter results.

    if selectreference is True, the following algorithm chooses a reference
    frame:

    1. from frames, get zenith distances (Z), moon distance (MOONDIST), moon
       elevation (MOONELEV).
    2. from (text) fiphot files, get nsources with G flags, MAD of magnitudes in
       largest aperture, median magnitude error in largest aperture.
    3. sort Z asc, MOONDIST desc, MOONELEV asc, nsources desc, good nsources
       desc, mag MAD asc, magerr asc, use these sort indices to sort framelists
    4. create sets of first 200 frames in each framelist above (or all frames if
       total nframes with phot < 200), and then intersect them all to find a set
       of frames that have the best properties. pick the first one out of the
       final intersection set as the reference frame.

    Returns a dict of the form:

    {'framelist':<list of frames passing fitsglob and headerfilters>,
     'photlist':<list of fitphot files associated with selected frames>,
     'framelistfile':<path to the outlistfile if created>,
     'referenceframe':<path to the chosen reference frame>,
     'referencestats':<stats dictionary for the reference frame>,
     'framestats':<a stats dictionary for all the frames>}

    '''

    # first, get the frames
    fitslist = glob.glob(os.path.join(fitsdir, fitsglob))

    print('%sZ: %s FITS files found in %s matching glob %s' %
          (datetime.utcnow().isoformat(),
           len(fitslist), fitsdir, fitsglob))

    # TODO: add filtering by header keywords using headerfilters

    photfits = []
    photlist = []

    # associate the frames with their fiphot files
    for fits in fitslist:

        searchpath = os.path.join(
            photdir,
            os.path.basename(fits).strip('.fits.fz') + '.fiphot'
            )

        if os.path.exists(searchpath):
            photfits.append(fits)
            photlist.append(searchpath)

    # remove all frames with no fiphot files
    fitslist = photfits

    # check if the workdir exists. if it doesn't, create it
    if workdir and not os.path.exists(workdir):
        os.mkdir(workdir)
        print('%sZ: made work directory %s' %
              (datetime.utcnow().isoformat(),
               workdir))

    # if the workdir is not provided, use the photdir as workdir
    if not workdir:
        workdir = photdir

    # make sure the workdir has a stats and MPHOTREF directory
    if not os.path.exists(os.path.join(workdir, 'stats')):
        os.mkdir(os.path.join(workdir,'stats'))
    if not os.path.exists(os.path.join(workdir,'MPHOTREF')):
        os.mkdir(os.path.join(workdir,'MPHOTREF'))

    workfitslist, workphotlist = [], []

    # link the frames in fitsdir to workdir
    if linkfiles:
        # temporarily change the directory so symlinking works
        cwd = os.getcwd()

        try:

            os.chdir(workdir)
            for fits in fitslist:
                os.symlink(fits, os.path.basename(fits))
                workfitslist.append(
                    os.path.join(workdir, os.path.basename(fits))
                    )
            # change back at the end
            os.chdir(cwd)
            print('%sZ: linked frames to workdir %s' %
                  (datetime.utcnow().isoformat(),
                   workdir))

        except Exception as e:
            print('%sZ: linking fiphot files to workdir %s failed!' %
                  (datetime.utcnow().isoformat(),
                   workdir))
            os.chdir(cwd)
            raise

    # if the workdir != photdir, then link the fiphot files too
    if workdir != photdir and linkfiles:
        # temporarily change the directory so symlinking works
        cwd = os.getcwd()

        try:
            os.chdir(workdir)
            for fiphot in photlist:
                os.symlink(fiphot, os.path.basename(fiphot))
                workphotlist.append(
                    os.path.join(workdir, os.path.basename(fiphot))
                    )
            # change back at the end
            os.chdir(cwd)
            print('%sZ: linked fiphot files to workdir %s' %
                  (datetime.utcnow().isoformat(),
                   workdir))

        except Exception as e:
            print('%sZ: linking fiphot files to workdir %s failed!' %
                  (datetime.utcnow().isoformat(),
                   workdir))
            os.chdir(cwd)
            raise

    if not workfitslist:
        workfitslist = fitslist
    if not workphotlist:
        workphotlist = photlist

    # collect the minimum returndict
    returndict = {'origframelist':fitslist,
                  'origphotlist':photlist,
                  'workframelist':workfitslist,
                  'workphotlist':workphotlist,
                  'framelistfile':outlistfile}

    # find a reference frame
    goodframes, goodphots = [], []
    if selectreference:

        print('%sZ: selecting a reference frame...' %
              (datetime.utcnow().isoformat(),))

        zenithdist, moondist, moonelev = [], [], []
        ngoodobjects, medmagerr, magerrmad = [], [], []

        # for each FITS and fiphot combo, collect stats
        for fits, fiphot in zip(workfitslist, workphotlist):

            headerdata = imageutils.get_header_keyword_list(
                fits,
                ['Z','MOONDIST','MOONELEV']
                )

            # decide if the fiphot file is binary or not. read the first 600
            # bytes and look for the '--binary-output' text
            with open(fiphot,'rb') as fiphotf:
                header = fiphotf.read(600)

            if '--binary-output' in header and HAVEBINPHOT:

                photdata_f = read_fiphot(fiphot)
                photdata = {
                    'mag':np.array(photdata_f['per aperture'][2]['mag']),
                    'err':np.array(photdata_f['per aperture'][2]['mag err']),
                    'flag':np.array(
                        photdata_f['per aperture'][2]['status flag']
                        )
                    }

            elif '--binary-output' in header and not HAVEBINPHOT:
                print('%sZ: %s is a binary fiphot file, '
                      'but no binary fiphot reader is present, skipping...' %
                      (datetime.utcnow().isoformat(), fiphot))
            else:

                # read in the phot file
                photdata = np.genfromtxt(
                    fiphot,
                    usecols=(12,13,14),
                    dtype='f8,f8,S5',
                    names=['mag','err','flag']
                    )

            # calculate stats for photometry

            # find good frames
            if '--binary-output' in header:
                goodind = np.where(photdata['flag'] == 0)
            else:
                goodind = np.where(photdata['flag'] == 'G')

            ngood = len(goodind[0])
            median_mag = np.median(photdata['mag'][goodind])
            median_magerr = np.median(photdata['err'][goodind])
            medabsdev_mag = np.median(
                np.abs(photdata['mag'][goodind] - median_mag)
                )

            # append to result lists
            goodframes.append(fits)
            goodphots.append(fiphot)
            zenithdist.append(headerdata['Z'])
            moondist.append(headerdata['MOONDIST'])
            moonelev.append(headerdata['MOONELEV'])
            ngoodobjects.append(ngood)
            medmagerr.append(median_magerr)
            magerrmad.append(medabsdev_mag)

            if DEBUG:
                print('frame = %s, phot = %s, Z = %s, MOONDIST = %s, '
                      'MOONELEV = %s, ngood = %s, medmagerr = %.5f, '
                      'magerrmad = %.5f' %
                      (fits, fiphot,
                       headerdata['Z'],
                       headerdata['MOONDIST'],
                       headerdata['MOONELEV'],
                       ngood, median_magerr, medabsdev_mag))

        # now that we're done collecting data, sort them in orders we want
        goodframes = np.array(goodframes)
        goodphots = np.array(goodphots)
        zenithdist_ind = np.argsort(zenithdist)
        moondist_ind = np.argsort(moondist)[::-1]
        moonelev_ind = np.argsort(moonelev)
        ngood_ind = np.argsort(ngoodobjects)[::-1]
        mederr_ind = np.argsort(medmagerr)
        magmad_ind = np.argsort(magerrmad)

        # get the first 200 of these or all 200 if n < 200
        if len(goodframes) > 200:
            zenithdist_ind = zenithdist_ind[:200]
            moondist_ind = moondist_ind[:200]
            moonelev_ind = moonelev_ind[:200]

            ngood_ind = ngood_ind[:200]
            mederr_ind = mederr_ind[:200]
            magmad_ind = magmad_ind[:200]

        # intersect all arrays to find a set of common indices that belong to
        # the likely reference frames
        photgood_ind = np.intersect1d(np.intersect1d(ngood_ind,
                                                     magmad_ind,
                                                     assume_unique=True),
                                      mederr_ind,assume_unique=True)

        headergood_ind =  np.intersect1d(np.intersect1d(moondist_ind,
                                                        moonelev_ind,
                                                        assume_unique=True),
                                         zenithdist_ind,assume_unique=True)

        allgood_ind = np.intersect1d(photgood_ind, headergood_ind,
                                     assume_unique=True)

        # if the headers and photometry produce a good reference frame, use
        # that. if they don't, use the photometry to choose a good reference
        # frame
        if len(allgood_ind) > 0:
            selectedreference = goodframes[allgood_ind[0]]
            selectedind = allgood_ind[0]
        elif len(photgood_ind) > 0:
            selectedreference = goodframes[photgood_ind[0]]
            selectedind = photgood_ind[0]
        elif len(headergood_ind) > 0:
            selectedreference = goodframes[headergood_ind[0]]
            selectedind = headergood_ind[0]
        else:
            selectedreference = None
            selectedind = None

        print('%sZ: selected reference frame = %s' %
              (datetime.utcnow().isoformat(), selectedreference))
        print('%sZ: selected reference phot = %s' %
              (datetime.utcnow().isoformat(), goodphots[selectedind]))

        # update the returndict with reference frame and stats
        returndict['referenceframe'] = selectedreference
        returndict['referencephot'] = goodphots[selectedind]
        if selectedreference:
            returndict['referencestats'] = {
                'zenithdist':zenithdist[selectedind],
                'moondist':moondist[selectedind],
                'moonelev':moonelev[selectedind],
                'ngood':ngoodobjects[selectedind],
                'magmad':magerrmad[selectedind],
                'mederr':medmagerr[selectedind],
                }
        else:
            returndict['referencestats'] = None

        # add all frame stats to the returndict
        if framestats:
            returndict['framestats'] = {'frame':goodframes,
                                        'phot':goodphots,
                                        'zenithdist':zenithdist,
                                        'moondist':moondist,
                                        'moonelev':moonelev,
                                        'ngood':ngood,
                                        'magmad':magerrmad,
                                        'mederr':medmagerr}
        # done with reference frame selection #


    # update the workfitslist and workphotlist if we did reference frame
    # selection
    if (len(goodframes) > 0 and len(goodphots) > 0 and
        len(goodframes) == len(goodphots)):
        returndict['workframelist'] = goodframes
        returndict['workphotlist'] = goodphots

    # write the framelistfile using the new frame locations
    if outlistfile:
        outf = open(outlistfile,'wb')
        for fitsline, photline in zip(workfitslist,workphotlist):
            outf.write('%s %s\n' % (photline, fitsline))
        outf.close()

        print('%sZ: wrote good frame list to %s' %
              (datetime.utcnow().isoformat(), os.path.abspath(outlistfile)))

    print('%sZ: done with fitsdir = %s' %
          (datetime.utcnow().isoformat(), fitsdir))

    # at the end, return returndict
    return returndict



def textphot_links_to_binphot_links(workdir,
                                    binphotdir):
    '''
    This is used to convert the fiphot links in workdir (which are text fiphot
    files) to the equivalent binary fiphot links in the same directory. Useful
    only for MagnitudeFitting.py.

    '''

    # get a list of text fiphot links
    text_fiphot_links = glob.glob(os.path.join(workdir,'*.fiphot'))

    print('%sZ: found %s text fiphot links in %s' %
          (datetime.utcnow().isoformat(),
           len(text_fiphot_links), workdir))

    # temporarily change working directory to the workdir
    cwd = os.getcwd()
    os.chdir(workdir)

    for textphot in text_fiphot_links:

        binphot = os.path.join(os.path.abspath(binphotdir),
                               os.path.basename(textphot))

        if os.path.exists(binphot):
            print('%sZ: converting textphot link %s -> binphot link for %s' %
                  (datetime.utcnow().isoformat(),
                   os.path.basename(textphot), binphot))
            os.remove(os.path.abspath(textphot))
            os.symlink(binphot, os.path.basename(textphot))
        else:
            print('%sZ: textphot link has no '
                  'equivalent binphot file, removing it!' %
                  (datetime.utcnow().isoformat(),
                   os.path.basename(textphot),))
            os.remove(os.path.abspath(textphot))



def make_magfit_config(configoutfile,
                       phot_apertures=3,
                       phot_band='r',
                       phot_brightmag=0.0,
                       phot_faintmag=16.0,
                       phot_fovcatalog='',
                       singlephot_statdir='',
                       masterphot_statdir='',
                       masterphot_outdir=''):
    '''
    This creates a magfit config file.

    '''

    # open the output file object
    outf = open(configoutfile,'wb')

    # write the top of the config file
    outf.write('num_apertures=%i\n' % phot_apertures)

    #
    # [mcat] section
    #
    outf.write('\n[mcat]\n')
    outf.write('rawphot_ver=0\n')
    # the mcat stat template entry
    templatestr = phot_fovcatalog
    outf.write("template=Template('%s')\n" % templatestr)
    outf.write('faint_mag=%.1f\n' %  phot_faintmag)
    outf.write('bright_mag=%.1f\n' %  phot_brightmag)
    outf.write('round_pointing=1\n')
    outf.write(("columns="
                "['id','ra','dec','xi','eta','J','K','qlt','I','%s']\n") %
               phot_band)
    outf.write('fov_alarm=20.0\n')
    outf.write('fov_safety_fac=1.1\n')

    #
    # [magfit.single] section
    #
    outf.write('\n[magfit.single]\n')
    # the magfit.single stat template entry
    templatestr = '/'.join(
        [singlephot_statdir,'SPHOTREFSTAT_object_${OBJECT}_${CMPOS}']
        )
    outf.write("stat_template=Template('%s')\n" % templatestr)
    outf.write("first_column=15\n")
    outf.write("reference_subpix=True\n")
    outf.write("file_extension='.sphotref'\n")
    outf.write("column_precision=1.0e-5\n")
    outf.write("column_name='mprmag'\n")
    outf.write("version=0\n")
    outf.write("count_weight=1.0\n")
    outf.write("error_avg='weightedmean'\n")
    outf.write("max_rej_iter=20\n")
    outf.write("rej_level=3.0\n")
    outf.write("max_mag_err=0.03\n")
    outf.write("noise_offset=0.0005\n")
    outf.write("ref_frame_stars=10000\n")
    outf.write("max_JmK=1.0\n")
    outf.write("min_JmK=0.0\n")
    outf.write("bright_mag_min=8.5\n")
    outf.write("AAAonly=True\n")
    outf.write("param_str='spatial:4;r:2,2;JmK:1,2;subpix:1,2'\n")
    outf.write("fntype=linear\n")

    #
    # [magfit.master] section
    #
    outf.write('\n[magfit.master]\n')
    # the magfit.master stat template entry
    templatestr = '/'.join(
        [masterphot_statdir,'MPHOTREFSTAT_object_${OBJECT}_${CMPOS}']
        )
    outf.write("stat_template=Template('%s')\n" % templatestr)
    outf.write("first_column=15\n")
    outf.write("reference_subpix=False\n")
    outf.write("file_extension='.mphotref'\n")
    outf.write("column_precision=1.0e-5\n")
    outf.write("column_name='mprmag'\n")
    outf.write("version=0\n")
    outf.write("count_weight=1.0\n")
    outf.write("error_avg='weightedmean'\n")
    outf.write("max_rej_iter=20\n")
    outf.write("rej_level=3.0\n")
    outf.write("max_mag_err=0.03\n")
    outf.write("noise_offset=0.0005\n")
    outf.write("ref_frame_stars=10000\n")
    outf.write("max_JmK=1.0\n")
    outf.write("min_JmK=0.0\n")
    outf.write("bright_mag_min=8.5\n")
    outf.write("AAAonly=True\n")
    outf.write("param_str='spatial:4;r:2,2;JmK:1,2;subpix:1,2'\n")
    outf.write("fntype='linear'\n")

    #
    # [mphotref] section
    #
    outf.write('\n[mphotref]\n')
    outf.write('version=0\n')
    # the mphotref template entry
    templatestr = '/'.join(
        [masterphot_outdir,'mphotref_${CMPOS}.AP${APERTURE}']
        )
    outf.write("template=Template('%s')" % templatestr)
    outf.write("rms_fit_bright_mag_min=9.0\n")
    outf.write("max_rms_quantile=0.1\n")
    outf.write("max_rms_above_fit=4\n")
    outf.write("rms_fit_rej_lvl=3.0\n")
    outf.write("rms_fit_err_avg='median'\n")
    outf.write("rms_fit_param='spatial:2;r:2,2'\n")
    outf.write("grcollect_tempdir='/dev/shm'\n")
    outf.write("min_meas_med=0.9\n")
    outf.write("min_measurements=0.01\n")
    outf.write("rej_iterations=20\n")
    outf.write("rej_outliers='iterations=20,median,meddev=8'\n")

    outf.close()

    print('%sZ: wrote magfit config to %s' %
          (datetime.utcnow().isoformat(),
           os.path.abspath(configoutfile),))

    return configoutfile



def make_fiphot_list(searchdirs,
                     searchglob,
                     listfile):
    '''This makes a list of fiphot files in the list of directories specified in
    searchdirs, using the searchglob to filter by filename. Returns a list of
    absolute paths to the fiphot files and writes this to the file specified in
    listfile.

    '''

    fiphotlist = []

    for fdir in searchdirs:
        fiphotfiles = glob.glob(fdir, searchglob)
        fiphotlist.extend([os.path.abspath(x) for x in fiphotfiles])

    # write to the output file
    outf = open(listfile,'wb')

    for fdir in fiphotlist:
        outf.write('%s\n' % fdir)

    outf.close()

    return listfile



def run_magfit(sphotref_frame,
               sphotref_phot,
               magfit_frame_list,
               magfit_config_file,
               magfit_type,
               nprocs=16,
               hatnetwork='HATSouth',
               magfitexec='MagnitudeFitting.py'):
    '''
    This runs magfit in single photometric reference mode.

    command to use is:

    first create the directories stats and MPHOTREF in the work directory

    python MagnitudeFitting.py HATSouth single /nfs/phs3/ar1/S/HP0/PHOT_WB/1-20141120-work-CCD5/1-377740d_5.fits /nfs/phs3/ar1/S/HP0/PHOT_WB/1-20141120-work-CCD5/1-377740d_5.fiphot -p 8 --config-file=/home/wbhatti/scripts/magfit_rG548.cfg --manual-frame-list=/nfs/phs3/ar1/S/HP0/PHOT_WB/1-20141120-work-CCD5/ccd5-20141120.framelist --stat

    '''

    if not (os.path.exists(magfit_frame_list) and
            os.path.exists(sphotref_frame) and
            os.path.exists(sphotref_phot) and
            os.path.exists(magfit_config_file)):

        print('%sZ: some required files are missing!' %
              (datetime.utcnow().isoformat(),))
        return None

    magfitcmd = MAGFITCMD.format(
        magfitexec=os.path.abspath(magfitexec),
        network=network,
        fit_type='single',
        sphotref_frame=sphotref_frame,
        sphotref_phot=sphotref_phot,
        nprocs=nprocs,
        magfit_config_file=configfile,
        magfit_frame_list=magfit_frame_list
        )

    if DEBUG:
        print(magfitcmd)

    print('%sZ: starting %s magfit with %s processes...' %
          (datetime.utcnow().isoformat(), fit_type, nprocs))

    # execute the magfit shell command
    magfitproc = subprocess.Popen(shlex.split(magfitcmd),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    # get results
    magfit_stdout, magfit_stderr = magfitproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if magfitproc.returncode == 0:

        print('%sZ: %s magfit completed.' %
              (datetime.utcnow().isoformat(), fit_type))

        return True

    else:

        print('%sZ: %s magfit failed!' %
              (datetime.utcnow().isoformat(), fit_type))

        return False



def get_master_photref(sphotref_frame,
                       fiphot_list,
                       magfit_config_file,
                       hatnetwork='HATSouth',
                       photrefexec='do_masterphotref.py'):
    '''
    Generates the master photometric reference from single ref photometry.

    python do_masterphotref.py HATSouth /nfs/phs3/ar1/S/HP0/PHOT_WB/1-20141120-work-CCD5/1-377740d_5.fits --manual-frame-list=/nfs/phs3/ar1/S/HP0/PHOT_WB/1-20141120-work-CCD5/1-20141120-CCD5-binphot.list --config-file=/home/wbhatti/scripts/magfit_rG548.cfg --nostat

    '''

    if not (os.path.exists(fiphot_list) and
            os.path.exists(sphotref_frame) and
            os.path.exists(magfit_config_file)):

        print('%sZ: some required files are missing!' %
              (datetime.utcnow().isoformat(),))
        return None

    photrefcmd = MPHOTREFCMD.format(
        mphotrefexec=os.path.abspath(photrefexec),
        network=network,
        sphotref_frame=sphotref_frame,
        fiphot_list=fiphot_list,
        photref_config_file=configfile
        )

    if DEBUG:
        print(photrefcmd)

    print('%sZ: starting photref...' %
          (datetime.utcnow().isoformat(),))

    # execute the photref shell command
    photrefproc = subprocess.Popen(shlex.split(photrefcmd),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    # get results
    photref_stdout, photref_stderr = photrefproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if photrefproc.returncode == 0:

        print('%sZ: photref completed.' %
              (datetime.utcnow().isoformat(),))

        return True

    else:

        print('%sZ: photref failed!' %
              (datetime.utcnow().isoformat(),))

        return False



def run_psf_trans():
    '''
    This runs the psf trans action. This probably gets the S, D, K values for
    the future EPD process. Not needed now since we get this from fistar in the
    current implementation.

    '''

#####################################
## LIGHTCURVE GENERATION FUNCTIONS ##
#####################################


def generate_lightcurves():
    '''
    This generates the raw lightcurves.

    Basically:

    0. dump the binary fiphot files to text fiphot files
    1. we need to get JDs for each frame
    2. get a master list of objects and the fiphot files they're in
    3. for each object, using linecache, get the lines for each measurement
       and concatenate into a output file along with JDs

    lc_gen.sh

    '''


###############################
## POST LIGHTCURVE FUNCTIONS ##
###############################

def run_epd():
    '''
    This runs EPD on the lightcurves from the pipeline.

    See do_epd.py for the coefficients and orders of the fit function.

    '''

def run_tfa():
    '''
    This runs TFA on the EPD lightcurves.

    use the commandline tfa program.

    '''

#########################
## REPORTING FUNCTIONS ##
#########################

# FIXME: break these out into their own module
# 1. RMS vs. MAG plot for EPD and TFA lightcurves
# 2. MAD vs. MAG plot for EPD and TFA lightcurves
# 3. ratios of RMS and MAD vs. MAG for CCD 6,7,8 to that of CCD 5
# 4. binned LC versions of these plots, using 10, 30, and 60 minute binning


##########################
## EXECUTABLE FUNCTIONS ##
##########################

# functions that roll up the parts of the pipeline into aggregate functions go
# here, along with a main function

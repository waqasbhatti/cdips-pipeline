#!/usr/bin/env python

'''
aperturephot.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Dec 2014

Contains aperture photometry routines for HATPI. Needs reduced frames.

The usual sequence is:

1. run parallel_extract_sources on all frames with threshold ~ 10000 to get
   bright stars for astrometry.

2. run parallel_anet to get WCS headers for all frames.

3. run make_fov_catalog to get a FOV source catalog for the field.

4. run reform_fov_catalog to cut this down to the columns needed for magfit
   only.

5. run parallel_fitsdir_photometry for photometry on all frames

6. run get_magfit_frames to select a single magfit photometry reference and set
   up per-CCD work directories, symlinks, etc. for the next steps.

7. run make_magfit_config to generate magfit config files for
   MagnitudeFitting.py

8. run make_fiphot_list to make lists of fiphot files for each CCD.

9. run MagnitudeFitting.py in single reference mode.

10. run do_masterphotref.py to get the master mag fit reference.

11. run MagnitudeFitting.py in master reference mode.

12. run parallel_collect_lightcurves to collect all lightcurves into .rlc files.

13. run serial_run_epd or parallel_run_epd to do EPD on all LCs.

14. run parallel_lc_statistics to collect stats on .epdlc files.

15. run choose_tfa_template to choose TFA template stars using the .epdlc stats.

16. run parallel_run_tfa for TFA to get .tfalc.TF{1,2,3} files (FIXME: still
    need to collect into single .tfalc files for all apertures)

17. run parallel_lc_statistics to collect stats on .tfalc files.

18. run parallel_bin_lightcurves to bin LCs to desired time-bins.

19. run parallel_binnedlc_statistics to collect stats for the binned LCs.

20. run plot_stats_file to make RMS vs. mag plots for all unbinned and binned
    LCs.

21. run plot_magrms_comparison to compare the mag-RMS relation for various CCDs.

22. run plot_ismphot_comparison to compare against ISM photometry statistics for
    the same field (requires common stars).

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
             "{sphotref_frame} {sphotref_phot} "
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
    if (anetproc.returncode == 0 and
        os.path.exists(os.path.abspath(wcsout)) and
        os.stat(os.path.abspath(wcsout)).st_size > 0):

        print('%sZ: anet WCS %s generated for frame sourcelist %s' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(wcsout), os.path.abspath(srclist)))

        return os.path.abspath(wcsout)

    # if astrometry did not succeed, complain and remove the zero-size wcs file
    # anet inexplicably makes anyway
    else:

        print('%sZ: anet WCS %s failed for frame sourcelist %s! Error was: %s' %
              (datetime.utcnow().isoformat(),
               os.path.abspath(wcsout),
               os.path.abspath(srclist),
               anet_stderr))

        # remove the broken wcs if astrometry failed
        if os.path.exists(os.path.abspath(wcsout)):
            os.remove(os.path.abspath(wcsout))

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
                  fistarglob='*_?.fistar',
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
    fistarlist = glob.glob(os.path.join(srclistdir, fistarglob))

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



##########################
## PHOTOMETRY FUNCTIONS ##
##########################

def make_fov_catalog(ra=None, dec=None, size=None,
                     brightrmag=6.0,
                     faintrmag=13.0,
                     fits=None,
                     outfile=None,
                     catalog='2MASS',
                     catalogpath=None,
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
                                                    catbox)

    print('%sZ: making FOV catalog for '
          'center RA, DEC = %.5f, %.5f with size = %.5f deg' %
          (datetime.utcnow().isoformat(),
           catra, catdec, catbox))


    catalogcmd = CATALOGS[catalog]['cmd'].format(
        ra=catra,
        dec=catdec,
        boxlen=catbox,
        catalogpath=catalogpath if catalogpath else CATALOGS[catalog]['path'],
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
                       columns='id,ra,dec,xi,eta,J,K,qlt,I,r,i,z'):
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
    fitslist = glob.glob(os.path.join(fitsdir,'*_?.fits'))

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

    # make sure the wcs file makes sense before trying the transform
    if (framewcsfile and
        os.path.exists(os.path.abspath(framewcsfile)) and
        os.stat(os.path.abspath(framewcsfile)).st_size > 0):

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

            # get indices for the lines to be kept optionally, remove all
            # sources within pixborders pixels of the edges of the image.
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

    # the wcs file doesn't make sense for this FITS image, complain and return
    else:
        print('%sZ: WCS transform file does not work for %s, '
              'skipping this frame...' %
              (datetime.utcnow().isoformat(), fits))
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

        else:

            outsourcelist = projcatfile


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

    try:

        # task[0] is args and first arg is the path to the fits file
        result = (task[0][0], do_photometry(*task[0], **task[1]))

    except Exception as e:
        print('photometry failed! reason: %s' % e)
        result = None

    return result


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
    fitslist = glob.glob(os.path.join(fitsdir,'*_?.fits'))

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
                del photdata_f

            elif '--binary-output' in header and not HAVEBINPHOT:

                print('%sZ: %s is a binary fiphot file, '
                      'but no binary fiphot reader is present, skipping...' %
                      (datetime.utcnow().isoformat(), fiphot))
                continue

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
            zenithdist_ind = zenithdist_ind[:500]
            moondist_ind = moondist_ind[:500]
            moonelev_ind = moonelev_ind[:500]

            ngood_ind = ngood_ind[:500]
            mederr_ind = mederr_ind[:500]
            magmad_ind = magmad_ind[:500]

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
                                        'ngood':ngoodobjects,
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
    outf.write("column_name='sprmag'\n")
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
    outf.write("template=Template('%s')\n" % templatestr)
    outf.write("rms_fit_bright_mag_min=9.0\n")
    outf.write("max_rms_quantile=0.1\n")
    outf.write("max_rms_above_fit=4\n")
    outf.write("rms_fit_rej_lvl=3.0\n")
    outf.write("rms_fit_err_avg=median\n")
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
        fiphotfiles = glob.glob(os.path.join(fdir, searchglob))
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
    This runs magfit in single/master photometric reference mode.

    lcohpsrv1 invocation:

    python /home/hatuser/wbhatti/src/MagnitudeFittingOrig.py HATSouth single /nfs/lcohpsrv1/ar1/scratch/PHOT_WB/projid8/ccd5-fits/1-404411d_5.fits /nfs/lcohpsrv1/ar1/scratch/PHOT_WB/projid8/ccd5-fits/1-404411d_5.fits -p 8 --log-config=/home/hatuser/wbhatti/src/logging.conf --config-file=/nfs/lcohpsrv1/ar1/scratch/PHOT_WB/projid8/photometry-ap/ccd5-magfit.cfg --manual-frame-list=/nfs/lcohpsrv1/ar1/scratch/PHOT_WB/projid8/photometry-ap/ccd5-magfit-frames.list --stat


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
        network=hatnetwork,
        fit_type=magfit_type,
        sphotref_frame=sphotref_frame,
        sphotref_phot=sphotref_phot,
        nprocs=nprocs,
        magfit_config_file=magfit_config_file,
        magfit_frame_list=magfit_frame_list
        )

    if DEBUG:
        print(magfitcmd)

    print('%sZ: starting %s magfit with %s processes...' %
          (datetime.utcnow().isoformat(), magfit_type, nprocs))

    # execute the magfit shell command
    magfitproc = subprocess.Popen(shlex.split(magfitcmd),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    # get results
    magfit_stdout, magfit_stderr = magfitproc.communicate()

    # get results if succeeded, log outcome, and return path of outfile
    if magfitproc.returncode == 0:

        print('%sZ: %s magfit completed.' %
              (datetime.utcnow().isoformat(), magfit_type))

        return True

    else:

        print('%sZ: %s magfit failed!' %
              (datetime.utcnow().isoformat(), magfit_type))
        print('%sZ: error returned was %s' % magfit_stderr)

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
        network=hatnetwork,
        sphotref_frame=sphotref_frame,
        fiphot_list=fiphot_list,
        magfit_config_file=magfit_config_file
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


#####################################
## LIGHTCURVE GENERATION FUNCTIONS ##
#####################################


def collect_lightcurve(hatid,
                       framefiles,
                       srcfiles,
                       photfiles,
                       jdlist,
                       outdir,
                       ignorecollected=True):
    '''
    This collects all photometric info into an LC for a given HATID. Returns
    path of collected LC.

    rjd    Reduced Julian Date (RJD = JD - 2400000.0)
    rstfc  Unique frame key ({STID}-{FRAMENUMBER}_{CCDNUM})
    xcc    X coordinate on CCD
    ycc    Y coordinate on CCD
    bgv    Background value
    bge    Background measurement error
    fsv    Measured S value
    fdv    Measured D value
    fkv    Measured K value
    im1    Instrumental magnitude in aperture 1
    ie1    Instrumental magnitude error for aperture 1
    iq1    Instrumental magnitude quality flag for aperture 1 (0 or G OK, X bad)
    im2    Instrumental magnitude in aperture 2
    ie2    Instrumental magnitude error for aperture 2
    iq2    Instrumental magnitude quality flag for aperture 2 (0 or G OK, X bad)
    im3    Instrumental magnitude in aperture 3
    ie3    Instrumental magnitude error for aperture 3
    iq3    Instrumental magnitude quality flag for aperture 3 (0 or G OK, X bad)
    rm1    Reduced fit magnitude in aperture 1 after magnitude fitting
    rm2    Reduced fit magnitude in aperture 2 after magnitude fitting
    rm3    Reduced fit magnitude in aperture 3 after magnitude fitting

    '''

    outfile = os.path.join(outdir, '%s.rlc' % hatid)

    # if this LC is already collected and ignorecollected is True, then don't
    # process it
    if ignorecollected and os.path.exists(outfile):
        outlines = open(outfile,'rb').readlines()
        if len(outlines) > 0:
            print('%sZ: LC for %s already collected with '
                  'ndet = %s, ignoring...' %
                  (datetime.utcnow().isoformat(), hatid, len(outlines)))
            return outfile

    outf = open(outfile, 'wb')

    outlineformat = ('%.8f %s %.3f %.3f %.3f %.3f %s %s %s '
                     '%.5f %.5f %s %.5f %.5f %s %.5f %.5f %s '
                     '%.5f %.5f %.5f\n')

    print('%sZ: collecting LC for %s...' %
          (datetime.utcnow().isoformat(), hatid))

    for frame, srcfile, photfile in zip(framefiles, srcfiles, photfiles):

        # get the JD and RSTFC key for this frame
        jd = jdlist[frame]
        rstfc = os.path.basename(frame).rstrip('.fits.fz')

        # get the sourcelist line for the hatid
        matchedsrcline = [x.split() for x in file(srcfile) if hatid in x]

        # if the hatid is found in the sourcelist, then try to find it in the
        # fiphot file
        if len(matchedsrcline) == 1:

            matchedsrcline = matchedsrcline[0]

            phot = read_fiphot(photfile)
            phot_hatids = ['HAT-%03i-%07i' % (x,y)
                           for x,y in zip(phot['field'],
                                          phot['source'])]

            # get the index of the HATID line in the fiphot file and get all the
            # photometric info
            try:

                phot_index = phot_hatids.index(hatid)

                objx, objy = phot['x'][phot_index], phot['y'][phot_index]
                objbgv, objbge = (phot['bg'][phot_index],
                                  phot['bg err'][phot_index])
                objs, objd, objk = (matchedsrcline[-3],
                                    matchedsrcline[-2],
                                    matchedsrcline[-1])

                objim1, objie1, objiq1 = (
                    phot['per aperture'][0]['mag'][phot_index],
                    phot['per aperture'][0]['mag err'][phot_index],
                    phot['per aperture'][0]['status flag'][phot_index]
                    )
                objim2, objie2, objiq2 = (
                    phot['per aperture'][1]['mag'][phot_index],
                    phot['per aperture'][1]['mag err'][phot_index],
                    phot['per aperture'][1]['status flag'][phot_index]
                    )
                objim3, objie3, objiq3 = (
                    phot['per aperture'][2]['mag'][phot_index],
                    phot['per aperture'][2]['mag err'][phot_index],
                    phot['per aperture'][2]['status flag'][phot_index]
                    )

                objrm1, objrm2, objrm3 = (phot['mprmag[0]'][phot_index],
                                          phot['mprmag[1]'][phot_index],
                                          phot['mprmag[2]'][phot_index])

                # format the line string and write it to the outfile
                outline = outlineformat % (
                    jd, rstfc, objx, objy, objbgv, objbge, objs, objd, objk,
                    objim1, objie1, objiq1,
                    objim2, objie2, objiq2,
                    objim3, objie3, objiq3,
                    objrm1, objrm2, objrm3
                    )
                outf.write(outline)

            # if the hatid isn't in the fiphot file, we can't do anything. skip
            # to the next hatid
            except ValueError as e:
                continue

    # at the end of collection
    print('%sZ: collected LC for %s' % (datetime.utcnow().isoformat(), hatid))

    outf.close()
    return outfile


def serial_collect_lightcurves(finalsourcesfile,
                               fitsdir,
                               fitsglob,
                               fiphotdir,
                               sourcelistdir,
                               outdir):
    '''
    This generates the raw lightcurves.

    finalsourcesfile = the file from which to take the list of final sources to
                       extract from the fiphot files. the sphotref/mphotref stat
                       files appear to be the right ones to use

    fitsdir = directory to search for FITS files

    fitsglob = glob to apply in search for FITS files

    fiphotdir = directory to search for .fiphot files matching FITS files. the
                xy coordinates, background measurements, instrumental magnitudes
                for each aperture, and the reduced mags for each aperture will
                be taken from here.

    sourcelistdir = directory to search for .sourcelist files matching FITS
                    files. the S, D, K values will be taken from here


    outdir = directory where the output .rlc LC files will go, one per HAT
             object.

    Basically:

    0. dump the binary fiphot files to text fiphot files
    1. we need to get JDs for each frame
    2. get a master list of objects and the fiphot files they're in
    3. for each object, using linecache, get the lines for each measurement
       and concatenate into a output file along with JDs

    lc_gen.sh

    The output LC format is:

    jd rstfc x y bgv bge s d k
    im1 ie1 iq1 im2 ie2 iq2 im3 ie3 iq3
    rm1 rm2 rm3 ep1 ep2 ep3 tf1 tf2 tf3


    '''

    # get all the FITS files
    fitsfiles = glob.glob(os.path.join(fitsdir, fitsglob))

    print('%sZ: found %s FITS files matching %s, getting JDs from headers...' %
          (datetime.utcnow().isoformat(), len(fitsfiles), fitsglob))

    jds = {x:get_header_keyword(x,'JD') for x in fitsfiles}

    # list the .fiphot files
    fiphotfiles = [
        os.path.join(fiphotdir,
                     '%s.fiphot' % os.path.basename(x).strip('.fits.fz'))
        for x in fitsfiles
        ]

    # list the .sourcelist files
    sourcelistfiles = [
        os.path.join(sourcelistdir,
                     '%s.sourcelist' % os.path.basename(x).strip('.fits.fz'))
        for x in fitsfiles
        ]

    # filter all the FITs, .sourcelist and fiphot files to make sure all three
    # files exist for a frame
    final_fits_files, final_fiphot_files, final_sourcelist_files = [], [], []

    for fits, fiphot, sourcelist in zip(fitsfiles, fiphotfiles, sourcelistfiles):

        if (os.path.exists(fits) and
            os.path.exists(fiphot) and
            os.path.exists(sourcelist)):

            final_fits_files.append(fits)
            final_fiphot_files.append(fiphot)
            final_sourcelist_files.append(sourcelist)

        else:

            print('%sZ: sourcelist/fiphot missing for %s, ignoring...' %
                  (datetime.utcnow().isoformat(), fits))


    # read the finalsourcesfile and get the HAT IDs to use
    hatid_list = np.loadtxt(finalsourcesfile,usecols=(0,),dtype='S17')

    print('%sZ: %s objects to process, starting collection...' %
          (datetime.utcnow().isoformat(), len(hatid_list)))

    # create a task list for the parallel collection processes
    collect_tasks = [(str(x), final_fits_files, final_sourcelist_files,
                      final_fiphot_files, jds, outdir)
                     for x in hatid_list]

    # make the output directory if not ready
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # launch the collection
    for task in collect_tasks:
        result = collect_lightcurve(*task)

    print('%sZ: done. %s objects collected into lightcurves.' %
          (datetime.utcnow().isoformat(), len(hatid_list)))



def lightcurve_collection_worker(task):
    '''
    Wrapper around collect_lightcurve for use in parallel_collect_lightcurves
    below.

    '''
    try:

        result =  collect_lightcurve(*task)

    except Exception as e:

        print('%sZ: failed to collect LC for %s, error was %s' %
              (datetime.utcnow().isoformat(), task[0], e))
        result = None

    return result


def parallel_collect_lightcurves(finalsourcesfile,
                                 fitsdir,
                                 fitsglob,
                                 fiphotdir,
                                 sourcelistdir,
                                 outdir,
                                 nworkers=16,
                                 workerntasks=500):
    '''
    This generates the raw lightcurves.

    finalsourcesfile = the file from which to take the list of final sources to
                       extract from the fiphot files. the sphotref/mphotref stat
                       files appear to be the right ones to use

    fitsdir = directory to search for FITS files

    fitsglob = glob to apply in search for FITS files

    fiphotdir = directory to search for .fiphot files matching FITS files. the
                xy coordinates, background measurements, instrumental magnitudes
                for each aperture, and the reduced mags for each aperture will
                be taken from here.

    sourcelistdir = directory to search for .sourcelist files matching FITS
                    files. the S, D, K values will be taken from here


    outdir = directory where the output .rlc LC files will go, one per HAT
             object.

    Basically:

    0. dump the binary fiphot files to text fiphot files
    1. we need to get JDs for each frame
    2. get a master list of objects and the fiphot files they're in
    3. for each object, using linecache, get the lines for each measurement
       and concatenate into a output file along with JDs

    lc_gen.sh

    The output LC format is:

    jd rstfc x y bgv bge s d k
    im1 ie1 iq1 im2 ie2 iq2 im3 ie3 iq3
    rm1 rm2 rm3 ep1 ep2 ep3 tf1 tf2 tf3

    '''

    # get all the FITS files
    fitsfiles = glob.glob(os.path.join(fitsdir, fitsglob))

    print('%sZ: found %s FITS files matching %s, getting JDs from headers...' %
          (datetime.utcnow().isoformat(), len(fitsfiles), fitsglob))

    jds = {x:get_header_keyword(x,'JD') for x in fitsfiles}

    # list the .fiphot files
    fiphotfiles = [
        os.path.join(fiphotdir,
                     '%s.fiphot' % os.path.basename(x).strip('.fits.fz'))
        for x in fitsfiles
        ]

    # list the .sourcelist files
    sourcelistfiles = [
        os.path.join(sourcelistdir,
                     '%s.sourcelist' % os.path.basename(x).strip('.fits.fz'))
        for x in fitsfiles
        ]

    # filter all the FITs, .sourcelist and fiphot files to make sure all three
    # files exist for a frame
    final_fits_files, final_fiphot_files, final_sourcelist_files = [], [], []

    for fits, fiphot, sourcelist in zip(fitsfiles, fiphotfiles, sourcelistfiles):

        if (os.path.exists(fits) and
            os.path.exists(fiphot) and
            os.path.exists(sourcelist)):

            final_fits_files.append(fits)
            final_fiphot_files.append(fiphot)
            final_sourcelist_files.append(sourcelist)

        else:

            print('%sZ: sourcelist/fiphot missing for %s, ignoring...' %
                  (datetime.utcnow().isoformat(), fits))

    # read the finalsourcesfile and get the HAT IDs to use
    hatid_list = np.loadtxt(finalsourcesfile,usecols=(0,),dtype='S17')

    print('%sZ: %s objects to process, starting collection...' %
          (datetime.utcnow().isoformat(), len(hatid_list)))

    # create a task list for the parallel collection processes
    collect_tasks = [(str(x), final_fits_files, final_sourcelist_files,
                      final_fiphot_files, jds, outdir)
                     for x in hatid_list]

    # make the output directory if not ready
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # make the output directory if not ready
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    results = pool.map(lightcurve_collection_worker, collect_tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s objects collected into lightcurves.' %
          (datetime.utcnow().isoformat(), len(hatid_list)))

    return {x:y for x,y in zip(hatid_list, results)}


###################
## EPD FUNCTIONS ##
###################

def epd_diffmags(coeff, fsv, fdv, fkv, xcc, ycc, bgv, bge, mag):
    '''
    This calculates the difference in mags after EPD coefficients are
    calculated.

    final EPD mags = median(magseries) + epd_diffmags()

    '''

    return -(coeff[0]*fsv**2. +
             coeff[1]*fsv +
             coeff[2]*fdv**2. +
             coeff[3]*fdv +
             coeff[4]*fkv**2. +
             coeff[5]*fkv +
             coeff[6] +
             coeff[7]*fsv*fdv +
             coeff[8]*fsv*fkv +
             coeff[9]*fdv*fkv +
             coeff[10]*np.sin(2*np.pi*xcc) +
             coeff[11]*np.cos(2*np.pi*xcc) +
             coeff[12]*np.sin(2*np.pi*ycc) +
             coeff[13]*np.cos(2*np.pi*ycc) +
             coeff[14]*np.sin(4*np.pi*xcc) +
             coeff[15]*np.cos(4*np.pi*xcc) +
             coeff[16]*np.sin(4*np.pi*ycc) +
             coeff[17]*np.cos(4*np.pi*ycc) +
             coeff[18]*bgv +
             coeff[19]*bge -
             mag)


def epd_magseries(mag, fsv, fdv, fkv, xcc, ycc, bgv, bge,
                  smooth=21, sigmaclip=3.0):
    '''
    Detrends a magnitude series given in mag using accompanying values of S in
    fsv, D in fdv, K in fkv, x coords in xcc, y coords in ycc, background in
    bgv, and background error in bge. smooth is used to set a smoothing
    parameter for the fit function. Does EPD voodoo.

    '''

    # find all the finite values of the magnitude
    finiteind = np.isfinite(mag)

    # calculate median and stdev
    mag_median = np.median(mag[finiteind])
    mag_stdev = np.nanstd(mag)

    # if we're supposed to sigma clip, do so
    if sigmaclip:
        excludeind = abs(mag - mag_median) < sigmaclip*mag_stdev
        finalind = finiteind & excludeind
    else:
        finalind = finiteind

    final_mag = mag[finalind]
    final_len = len(final_mag)

    if DEBUG:
        print('final epd fit mag len = %s' % final_len)

    # smooth the signal
    smoothedmag = medfilt(final_mag, smooth)

    # make the linear equation matrix
    epdmatrix = np.c_[fsv[finalind]**2.0,
                      fsv[finalind],
                      fdv[finalind]**2.0,
                      fdv[finalind],
                      fkv[finalind]**2.0,
                      fkv[finalind],
                      np.ones(final_len),
                      fsv[finalind]*fdv[finalind],
                      fsv[finalind]*fkv[finalind],
                      fdv[finalind]*fkv[finalind],
                      np.sin(2*np.pi*xcc[finalind]),
                      np.cos(2*np.pi*xcc[finalind]),
                      np.sin(2*np.pi*ycc[finalind]),
                      np.cos(2*np.pi*ycc[finalind]),
                      np.sin(4*np.pi*xcc[finalind]),
                      np.cos(4*np.pi*xcc[finalind]),
                      np.sin(4*np.pi*ycc[finalind]),
                      np.cos(4*np.pi*ycc[finalind]),
                      bgv[finalind],
                      bge[finalind]]

    # solve the equation epdmatrix * x = smoothedmag
    # return the EPD differential mags if the solution succeeds
    try:

        coeffs, residuals, rank, singulars = lstsq(epdmatrix, smoothedmag)

        if DEBUG:
            print('coeffs = %s, residuals = %s' % (coeffs, residuals))

        return epd_diffmags(coeffs, fsv, fdv, fkv, xcc, ycc, bgv, bge, mag)

    # if the solution fails, return nothing
    except Exception as e:

        print('%sZ: EPD solution did not converge! Error was: %s' %
              (datetime.utcnow().isoformat(), e))
        return None


def epd_lightcurve(rlcfile,
                   mags=[18,19,20],
                   sdk=[6,7,8],
                   xy=[2,3],
                   backgnd=[4,5],
                   smooth=21,
                   sigmaclip=3.0,
                   rlcext='rlc',
                   outfile=None):
    '''
    Runs the EPD process on rlcfile, using columns specified to get the required
    parameters. If outfile is None, the .epdlc will be placeed in the same
    directory as rlcfile.

    '''

    # read the lightcurve in
    rlc = np.genfromtxt(rlcfile,
                        usecols=tuple(xy + backgnd + sdk + mags),
                        dtype='f8,f8,f8,f8,f8,f8,f8,f8,f8,f8',
                        names=['xcc','ycc','bgv','bge','fsv','fdv','fkv',
                               'rm1','rm2','rm3'])

    # calculate the EPD differential mags
    epddiffmag1 = epd_magseries(rlc['rm1'],rlc['fsv'],rlc['fdv'],rlc['fkv'],
                                rlc['xcc'],rlc['ycc'],rlc['bgv'],rlc['bge'],
                                smooth=smooth, sigmaclip=sigmaclip)
    epddiffmag2 = epd_magseries(rlc['rm2'],rlc['fsv'],rlc['fdv'],rlc['fkv'],
                                rlc['xcc'],rlc['ycc'],rlc['bgv'],rlc['bge'],
                                smooth=smooth, sigmaclip=sigmaclip)
    epddiffmag3 = epd_magseries(rlc['rm3'],rlc['fsv'],rlc['fdv'],rlc['fkv'],
                                rlc['xcc'],rlc['ycc'],rlc['bgv'],rlc['bge'],
                                smooth=smooth, sigmaclip=sigmaclip)

    # add the EPD diff mags back to the median mag to get the EPD mags
    if epddiffmag1 is not None:
        mag_median = np.median(rlc['rm1'][np.isfinite(rlc['rm1'])])
        epdmag1 = epddiffmag1 + mag_median
    else:
        epdmag1 = np.array([np.nan for x in rlc['rm1']])
        print('%sZ: no EP1 mags available for %s!' %
              (datetime.utcnow().isoformat(), rlcfile))

    if epddiffmag2 is not None:
        mag_median = np.median(rlc['rm2'][np.isfinite(rlc['rm2'])])
        epdmag2 = epddiffmag2 + mag_median
    else:
        epdmag2 = np.array([np.nan for x in rlc['rm2']])
        print('%sZ: no EP2 mags available for %s!' %
              (datetime.utcnow().isoformat(), rlcfile))

    if epddiffmag3 is not None:
        mag_median = np.median(rlc['rm3'][np.isfinite(rlc['rm3'])])
        epdmag3 = epddiffmag3 + mag_median
    else:
        epdmag3 = np.array([np.nan for x in rlc['rm3']])
        print('%sZ: no EP3 mags available for %s!' %
              (datetime.utcnow().isoformat(), rlcfile))

    # now write the EPD LCs out to the outfile
    if not outfile:
        outfile = '%s.epdlc' % rlcfile.strip('.%s' % rlcext)

    inf = open(rlcfile,'rb')
    inflines = inf.readlines()
    inf.close()
    outf = open(outfile,'wb')

    for line, epd1, epd2, epd3 in zip(inflines, epdmag1, epdmag2, epdmag3):
        outline = '%s %.6f %.6f %.6f\n' % (line.rstrip('\n'), epd1, epd2, epd3)
        outf.write(outline)

    outf.close()
    return outfile


def serial_run_epd(rlcdir,
                   rlcglob='*.rlc',
                   outdir=None,
                   smooth=21,
                   sigmaclip=3.0):
    '''
    This runs EPD on the lightcurves from the pipeline.

    '''

    if not outdir:
        outdir = rlcdir

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    rlcfiles = glob.glob(os.path.join(rlcdir, rlcglob))

    for rlc in rlcfiles:

        outepd = os.path.join(outdir,
                              os.path.basename(rlc).replace('.rlc','.epdlc'))

        print('%sZ: doing EPD for %s...' %
              (datetime.utcnow().isoformat(), rlc))

        try:
            outfilename = epd_lightcurve(rlc,
                                         outfile=outepd,
                                         smooth=smooth,
                                         sigmaclip=sigmaclip,
                                         rlcext=os.path.splitext(rlcglob)[-1])
        except Exception as e:
            print('EPD failed for %s, error was: %s' % (rlc, e))



def parallel_epd_worker():
    '''
    Function to wrap the epd_lightcurve function for use with mp.Pool.

    '''

def parallel_run_epd():
    '''
    This runs EPD in parallel on the lightcurves from the pipeline.

    '''


###################
## TFA FUNCTIONS ##
###################

def choose_tfa_template(statsfile,
                        fovcatalog,
                        epdlcdir,
                        fovcat_idcol=0,
                        fovcat_xicol=3,
                        fovcat_etacol=4,
                        fovcat_magcol=9,
                        min_ndet=100,
                        min_nstars=50,
                        max_nstars=1000,
                        brightest_mag=8.5,
                        faintest_mag=12.0,
                        max_rms=0.1,
                        max_sigma_above_rmscurve=4.0,
                        outprefix=None,
                        tfastage1=True):
    '''This chooses suitable stars for TFA template purposes.

    statsfile = the file to use to get LC statistics from

    fovcatalog = the fovcatalog file, this must have xi and eta coordinates,
                 ras, decs, and magnitudes

    Returns a dict with lists of stars chosen, their stats, and filenames of
    where each star list was written.

    '''

    # read in the stats file
    stats = read_stats_file(statsfile)

    # read in the fovcatalog
    fovcat = np.genfromtxt(fovcatalog,
                           usecols=(fovcat_idcol,
                                    fovcat_xicol,
                                    fovcat_etacol,
                                    fovcat_magcol),
                           dtype='S17,f8,f8,f8',
                           names=['objid','xi','eta','mag'])

    # figure out the number of stars to use in the initial TFA template
    # number of stars = TFA_TEMPLATE_FRACTION * median ndet
    TFA_TEMPLATE_FRACTION = 0.1

    # 1. ndet >= median_ndet
    # 2. max rms <= 0.1
    # 3. brightest_mag < median_mag < faintest_mag
    # 4. fit rms-mag, then discard anything above max_sigma_above_rmscurve
    # find the objects in the fovcat and stats file that match these
    # conditions, then pick up to 1000 random stars

    outdict = {'statsfile':os.path.abspath(statsfile),
               'fovcat':os.path.abspath(fovcatalog),
               'lcdir':os.path.abspath(epdlcdir)}

    # do this per aperture
    for aperture in [1,2,3]:

        outdict[aperture] = {}

        # first, pick the stars that meet our stats requirements

        epdstr = 'ep%s' % aperture

        objid_col = 'lcobj'
        median_mag_col = 'med_sc_%s' % epdstr
        mad_mag_col = 'mad_sc_%s' % epdstr
        ndet_col = 'ndet_sc_%s' % epdstr

        objectid = stats[objid_col]
        mags_median = stats[median_mag_col]
        mags_mad = stats[mad_mag_col]
        obj_ndet = stats[ndet_col]

        goodind = np.isfinite(mags_median) & np.isfinite(mags_mad)
        objectid = objectid[goodind]
        mags_median = mags_median[goodind]
        mags_mad = mags_mad[goodind]
        obj_ndet = obj_ndet[goodind]

        print('\naperture %s: total good objects = %s' % (aperture,
                                                        len(objectid)))

        median_ndet = scipy.stats.nanmedian(obj_ndet)
        print('aperture %s: median ndet = %s' % (aperture, median_ndet))
        print('aperture %s: target TFA template size = %s' %
              (aperture, int(median_ndet*TFA_TEMPLATE_FRACTION)))
        outdict[aperture]['target_tfa_nstars'] = (
            median_ndet*TFA_TEMPLATE_FRACTION
        )

        stars_ndet_condition = obj_ndet >= median_ndet
        print('aperture %s: objects with ndet condition = %s' %
              (aperture, len(objectid[stars_ndet_condition])))

        stars_rms_condition = mags_mad < max_rms
        print('aperture %s: objects with rms condition = %s' %
              (aperture, len(objectid[stars_rms_condition])))

        rmsfit_condition = stars_ndet_condition
        print('aperture %s: objects with rmsfit condition = %s' %
              (aperture, len(objectid[rmsfit_condition])))

        # selection 1: fit a parabola to the median mag - mag MAD relation and
        #              reject all stars with RMS > max_sigma_above_rmscurve
        polyfit_coeffs = np.polyfit(mags_median[rmsfit_condition],
                             mags_mad[rmsfit_condition],
                             2)

        print('aperture %s: rms fit params = %s' % (aperture,polyfit_coeffs))

        # generate the model RMS curve with fit parameters
        model_rms = (polyfit_coeffs[0]*mags_median*mags_median +
                     polyfit_coeffs[1]*mags_median +
                     polyfit_coeffs[2])

        # find objects that lie below the requested rms threshold from this
        # curve
        threshold_condition = (mags_mad/model_rms) < max_sigma_above_rmscurve
        print('aperture %s: objects with threshold condition = %s' %
              (aperture, len(objectid[threshold_condition])))

        final_statistics_ind = (threshold_condition &
                                stars_rms_condition &
                                stars_ndet_condition)

        print('aperture %s: stars with good stats = %s' % (
            aperture,
            len(objectid[final_statistics_ind])
        ))

        good_stats_objects = objectid[final_statistics_ind]

        # selection 2: get the stars with good magnitudes
        mag_condition = ((fovcat['mag'] < faintest_mag) &
                         (fovcat['mag'] > brightest_mag))
        good_mag_objects = fovcat['objid'][mag_condition]
        print('aperture %s: stars with good mags = %s' %
              (aperture,len(good_mag_objects)))

        # finally, intersect these two arrays to find a set of good TFA objects
        good_tfa_objects = np.intersect1d(good_stats_objects,
                                          good_mag_objects)
        print('aperture %s: stars suitable for TFA = %s' %
              (aperture,len(good_tfa_objects)))

        # put this initial list into the outdict
        outdict[aperture]['tfa_suitable_objects'] = good_tfa_objects

        # selection 3: pick up to 1000 random objects for TFA. we need at least
        # TFA_TEMPLATE_FRACTION * median ndet number of objects
        if ((len(good_tfa_objects) > 1000) and
            (len(good_tfa_objects) > TFA_TEMPLATE_FRACTION*median_ndet)):
            tfa_stars = nprand.choice(good_tfa_objects,
                                      replace=False,
                                      size=500)
        elif (TFA_TEMPLATE_FRACTION*median_ndet <= len(good_tfa_objects) <= 1000):
            tfa_stars = good_tfa_objects
        else:
            print("aperture %s: not enough stars suitable for TFA!" %
                  aperture)
            tfa_stars = None

        # now get these stars IDs, LC fnames, xis, etas, and other things needed
        # for the first stage of TFA (this will choose exactly
        # TFA_TEMPLATE_FRACTION*median_ndet stars to use as the template for the
        # final stage of TFA)
        if tfa_stars is not None:

            print('aperture %s: %s objects chosen as TFA templates for stage 1' %
                  (aperture,len(tfa_stars)))
            outdict[aperture]['tfa_chosen_objects'] = tfa_stars

            tfa_stars_catmag = []
            tfa_stars_statrms = []
            tfa_stars_statndet = []
            tfa_stars_xi = []
            tfa_stars_eta = []
            tfa_stars_lcfile = []

            print('aperture %s: getting object info...' %
                  (aperture,))

            # get the stats for these objects
            for i, tfaobj in enumerate(tfa_stars):

                # search for the LC file for this object and make sure it exists
                lcfile_searchpath = os.path.join(epdlcdir, '%s.epdlc' % tfaobj)

                if os.path.exists(lcfile_searchpath):

                    print('getting LC for %s' % tfaobj)

                    tfa_stars_lcfile.append(
                        os.path.abspath(lcfile_searchpath)
                    )
                    tfa_stars_catmag.append(
                        fovcat['mag'][fovcat['objid'] == tfaobj]
                    )
                    tfa_stars_statrms.append(
                        mags_mad[objectid == tfaobj]
                    )
                    tfa_stars_statndet.append(
                        obj_ndet[objectid == tfaobj]
                    )
                    tfa_stars_xi.append(
                        fovcat['xi'][fovcat['objid'] == tfaobj]
                    )
                    tfa_stars_eta.append(
                        fovcat['eta'][fovcat['objid'] == tfaobj]
                    )

                # if it doesn't, then add nans to the file
                else:

                    print('ERR! couldn\'t find an LC for %s' % tfaobj)

                    tfa_stars_lcfile.append(None)
                    tfa_stars_catmag.append(np.nan)
                    tfa_stars_statrms.append(np.nan)
                    tfa_stars_statndet.append(np.nan)
                    tfa_stars_xi.append(np.nan)
                    tfa_stars_eta.append(np.nan)

            outdict[aperture]['tfa_chosen_lcfile'] = tfa_stars_lcfile
            outdict[aperture]['tfa_chosen_mag'] = np.ravel(tfa_stars_catmag)
            outdict[aperture]['tfa_chosen_rms'] = np.ravel(tfa_stars_statrms)
            outdict[aperture]['tfa_chosen_ndet'] = np.ravel(tfa_stars_statndet)
            outdict[aperture]['tfa_chosen_xi'] = np.ravel(tfa_stars_xi)
            outdict[aperture]['tfa_chosen_eta'] = np.ravel(tfa_stars_eta)

        # if no TFA stars could be chosen, return Nones for this aperture
        else:
            outdict[aperture]['tfa_chosen_lcfile'] = None
            outdict[aperture]['tfa_chosen_objects'] = None
            outdict[aperture]['tfa_chosen_mag'] = None
            outdict[aperture]['tfa_chosen_rms'] = None
            outdict[aperture]['tfa_chosen_ndet'] = None


        # make the input file for TFA stage 1 for this aperture
        if not outprefix:
            outfile = os.path.abspath(
                os.path.join(os.getcwd(),
                             'tfa-stage1-input-aperture-%s.txt' % aperture)
            )
        else:
            outfile = os.path.abspath(
                os.path.join(outprefix,
                             'tfa-stage1-input-aperture-%s.txt' % aperture)
            )


        outf = open(outfile,'wb')

        outline = '%s %s %.6f %.6f %i %.6f %.6f\n'

        for objid, lcf, mag, rms, ndet, xi, eta in zip(
                outdict[aperture]['tfa_chosen_objects'],
                outdict[aperture]['tfa_chosen_lcfile'],
                outdict[aperture]['tfa_chosen_mag'],
                outdict[aperture]['tfa_chosen_rms'],
                outdict[aperture]['tfa_chosen_ndet'],
                outdict[aperture]['tfa_chosen_xi'],
                outdict[aperture]['tfa_chosen_eta']
        ):

            outf.write(outline % (objid, lcf, mag, rms, ndet, xi, eta))
        outf.close()
        print('aperture %s: wrote object info to %s' %
              (aperture, outfile))
        outdict[aperture]['info_file'] = os.path.abspath(outfile)

        # END OF PER APERTURE STUFF


    # run TFAs stage 1 if we're supposed to do so
    if tfastage1:

        print('\nrunning TFA stage 1...')

        # run TFA stage 1 to pick the good objects
        tfa_stage1 = run_tfa_stage1(outdict)

        for aperture in tfa_stage1:

            outdict[aperture]['stage1_templatelist'] = tfa_stage1[aperture]

            # make the TFA template list file for this aperture
            if outdict[aperture]['stage1_templatelist']:
                templatelistfname = os.path.join(
                    outdict['lcdir'],
                    'aperture-%s-tfa-template.list' % aperture
                )
                outf = open(templatelistfname,'wb')
                for tfaobjid in outdict[aperture]['stage1_templatelist']:

                    templatelc = os.path.join(
                        outdict['lcdir'],
                        tfaobjid + '.epdlc'
                    )

                    if os.path.exists(templatelc):
                        outf.write('%s\n' % os.path.abspath(templatelc))

                outf.close()
                outdict[aperture]['stage1_tfa_templatefile'] = (
                    templatelistfname
                )
                print('aperture %s: wrote TFA template list to %s' %
                      (aperture, templatelistfname))


    return outdict


def run_tfa_stage1(tfainfo):
    '''This just runs the TFA in fake mode to generate a list of template
    stars. Uses the tfainfo dict created in choose_tfa_template above.

    Communicates with the tfa program over pipe and then writes its output to a
    tfa input file for the next stage.

       tfa -r <REFFILE> --col-ref-id <REFFILE_IDCOLNUM>
                        --col-ref-x <REFFILE_XCOLNUM>
                        --col-ref-y <REFFILE_YCOLNUM>
                        -n <NTEMPLATES_TO_USE>
                        -T - (for stdout)
                        -i /dev/null (no input file?)

    '''

    tfa_stage1_results = {}

    for aperture in [1,2,3]:

        tfacmdstr = ("tfa -r {inputfile} --col-ref-id 1 "
                     "--col-ref-x 6 --col-ref-y 7 "
                     "-n {ntemplates} -T - -i /dev/null")

        tfacmd = tfacmdstr.format(
            inputfile=tfainfo[aperture]['info_file'],
            ntemplates=int(tfainfo[aperture]['target_tfa_nstars'])
        )

        print('aperture %s: starting TFA stage 1...' % aperture)

        if DEBUG:
            print(tfacmd)

        tfaproc = subprocess.Popen(shlex.split(tfacmd),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        tfa_stdout, tfa_stderr = tfaproc.communicate()

        # get results if succeeded, log outcome, and return path of outfile
        if tfaproc.returncode == 0 or tfa_stdout:
            tfaobjects = tfa_stdout.split('\n')
            tfaobjects = [x for x in tfaobjects if x.startswith('HAT-')]
            print('aperture %s: TFA stage 1 completed, %s templates selected' %
                  (aperture, len(tfaobjects)))
            tfa_stage1_results[aperture] = tfaobjects
        else:
            print('aperture %s: TFA stage 1 failed, error was: %s' %
                  (aperture, tfa_stderr))
            tfa_stage1_results[aperture] = None

    return tfa_stage1_results



def run_tfa_singlelc(epdlc,
                     templatefiles,
                     outfile=None,
                     epdlc_jdcol=0,
                     epdlc_magcol=(21,22,23),
                     template_sigclip=5.0,
                     epdlc_sigclip=5.0):
    '''This runs TFA for all apertures defined in epdlc_magcol for the input
    epdlc file, given an existing TFA template list in templatefile. If outfile
    is None, the output TFA LC will be in the same directory as epdlc but with
    an extension of .tfalc.

    '''

    tfacmdstr = ("tfa -i {epdlc} -t {templatefile} "
                 "--col-jd {epdlc_jdcol} "
                 "--col-mag {epdlc_magcol} "
                 "--col-mag-out {tfalc_magcol} "
                 "--epsilon-time 1e-9 "
                 "--join-by-time "
                 "--templ-filter --templ-outl-limit {template_sigclip} "
                 "--lc-filter --lc-outl-limit {epdlc_sigclip} "
                 "--log-svn-version "
                 "-o {out_tfalc}")

    if not outfile:
        outfile = epdlc.replace('.epdlc','.tfalc')

    tfalc_output = []

    # run tfa for each aperture
    for templatef, magcol, magind in zip(templatefiles,
                                         epdlc_magcol,
                                         range(len(epdlc_magcol))):

        in_jdcol = epdlc_jdcol + 1
        in_magcol = magcol + 1
        out_magcol = in_magcol + 3

        aperture_outfile = outfile + ('.TF%s' % (magind+1))

        tfacmd = tfacmdstr.format(epdlc=epdlc,
                                  templatefile=templatef,
                                  epdlc_jdcol=in_jdcol,
                                  epdlc_magcol=in_magcol,
                                  tfalc_magcol=out_magcol,
                                  template_sigclip=template_sigclip,
                                  epdlc_sigclip=epdlc_sigclip,
                                  out_tfalc=aperture_outfile)

        if DEBUG:
            print(tfacmd)

        # execute the tfa shell command
        tfaproc = subprocess.Popen(shlex.split(tfacmd),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # get results
        tfa_stdout, tfa_stderr = tfaproc.communicate()

        # get results if succeeded, log outcome, and return path of outfile
        if tfaproc.returncode == 0:

            tfalc_output.append(aperture_outfile)

        else:

            print('%sZ: aperture %s TFA failed for %s! Error was: %s' %
                  (datetime.utcnow().isoformat(), magind+1, epdlc), tfa_stderr)

            tfalc_output.append(None)

    return tfalc_output


def parallel_tfa_worker(task):
    '''
    This wraps run_tfa_singlelc above.

    '''

    try:

        result = run_tfa_singlelc(task[0], task[1], **task[2])

    except Exception as e:

        print('%sZ: TFA failed for %s! Error was: %s' %
              (datetime.utcnow().isoformat(), task[0], e))
        result = None

    return result


def parallel_run_tfa(lcdir,
                     templatefiles,
                     epdlc_glob='*.epdlc',
                     epdlc_jdcol=0,
                     epdlc_magcol=(21,22,23),
                     template_sigclip=5.0,
                     epdlc_sigclip=5.0,
                     nworkers=16,
                     workerntasks=500):
    '''
    This runs TFA on the EPD lightcurves.

    '''

    epdlcfiles = glob.glob(os.path.join(lcdir, epdlc_glob))

    tasks = [(x, templatefiles, {'epdlc_jdcol':epdlc_jdcol,
                                 'epdlc_magcol':epdlc_magcol,
                                 'template_sigclip':template_sigclip,
                                 'epdlc_sigclip':epdlc_sigclip})
             for x in epdlcfiles]

    print('%sZ: %s objects to process, starting parallel TFA...' %
          (datetime.utcnow().isoformat(), len(epdlcfiles)))

    pool = mp.Pool(nworkers, maxtasksperchild=workerntasks)
    results = pool.map(parallel_tfa_worker, tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s LCs processed.' %
          (datetime.utcnow().isoformat(), len(epdlcfiles)))

    return {x:y for x,y in zip(epdlcfiles, results)}


#############################
## LC STATISTICS FUNCTIONS ##
#############################

# FIXME: break these out into their own module
# 1. RMS vs. MAG plot for EPD and TFA lightcurves
# 2. MAD vs. MAG plot for EPD and TFA lightcurves
# 3. ratios of RMS and MAD vs. MAG for CCD 6,7,8 to that of CCD 5
# 4. binned LC versions of these plots, using 10, 30, and 60 minute binning

def get_lc_statistics(lcfile,
                      rmcols=[18,19,20],
                      epcols=[21,22,23],
                      tfcols=[24,25,26],
                      sigclip=4.0):
    '''
    This calculates the following statistics for the magnitude columns in the
    given lcfile.

    mean
    median
    MAD
    stdev

    # TODO
    iterative stdev after sigma-clipping (RMS)

    IMPORTANT: the lcfile is always the .epdlc file (which contains the rlc, and
    is used to derive the filenames of the tfalcs)


    '''

    try:

        # get the reduced magnitude columns
        (rm1, rm2, rm3,
         ep1, ep2, ep3) = np.genfromtxt(lcfile,
                                        usecols=tuple(rmcols + epcols),
                                        unpack=True)

        tf1 = np.genfromtxt(lcfile.replace('.epdlc','.tfalc.TF1'),
                            usecols=(24,), unpack=True)
        tf2 = np.genfromtxt(lcfile.replace('.epdlc','.tfalc.TF2'),
                            usecols=(24,), unpack=True)
        tf3 = np.genfromtxt(lcfile.replace('.epdlc','.tfalc.TF3'),
                            usecols=(24,), unpack=True)

    # if we don't have TF columns, cut down to RM and EP only
    except Exception as e:

        print('%sZ: no TFA mags available for %s!' %
              (datetime.utcnow().isoformat(), lcfile))

        try:

            (rm1, rm2, rm3,
             ep1, ep2, ep3) = np.genfromtxt(lcfile,
                                            usecols=tuple(rmcols + epcols),
                                            unpack=True)

            tf1, tf2, tf3 = [], [], []

        except Exception as e:

            print('%sZ: no EPD mags available for %s!' %
                  (datetime.utcnow().isoformat(), lcfile))

            rm1, rm2, rm3 = np.genfromtxt(lcfile,
                                          usecols=tuple(rmcols),
                                          unpack=True)

            ep1, ep2, ep3, tf1, tf2, tf3 = [], [], [], [], [], []



    # get statistics for each column

    # RM1
    if len(rm1) > 4:

        finiteind = np.isfinite(rm1)
        rm1 = rm1[finiteind]
        median_rm1 = np.median(rm1)
        mad_rm1 = np.median(np.fabs(rm1 - median_rm1))
        mean_rm1 = np.mean(rm1)
        stdev_rm1 = np.std(rm1)
        ndet_rm1 = len(rm1)

        if sigclip:
            sigclip_rm1, lo, hi = stats_sigmaclip(rm1,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_rm1 = np.median(sigclip_rm1)
            mad_sigclip_rm1 = np.median(np.fabs(sigclip_rm1 -
                                                median_sigclip_rm1))
            mean_sigclip_rm1 = np.mean(sigclip_rm1)
            stdev_sigclip_rm1 = np.std(sigclip_rm1)
            ndet_sigclip_rm1 = len(sigclip_rm1)

        else:
            median_sigclip_rm1 = np.nan
            mad_sigclip_rm1 = np.nan
            mean_sigclip_rm1 = np.nan
            stdev_sigclip_rm1 = np.nan
            ndet_sigclip_rm1 = np.nan

    else:

        median_rm1, mad_rm1, mean_rm1, stdev_rm1 = np.nan, np.nan, np.nan, np.nan
        ndet_rm1 = np.nan
        median_sigclip_rm1, mad_sigclip_rm1 = np.nan, np.nan
        mean_sigclip_rm1, stdev_sigclip_rm1 = np.nan, np.nan
        ndet_sigclip_rm1 = np.nan

    # RM2
    if len(rm2) > 4:

        finiteind = np.isfinite(rm2)
        rm2 = rm2[finiteind]
        median_rm2 = np.median(rm2)
        mad_rm2 = np.median(np.fabs(rm2 - median_rm2))
        mean_rm2 = np.mean(rm2)
        stdev_rm2 = np.std(rm2)
        ndet_rm2 = len(rm2)

        if sigclip:
            sigclip_rm2, lo, hi = stats_sigmaclip(rm2,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_rm2 = np.median(sigclip_rm2)
            mad_sigclip_rm2 = np.median(np.fabs(sigclip_rm2 -
                                                median_sigclip_rm2))
            mean_sigclip_rm2 = np.mean(sigclip_rm2)
            stdev_sigclip_rm2 = np.std(sigclip_rm2)
            ndet_sigclip_rm2 = len(sigclip_rm2)

        else:
            median_sigclip_rm2 = np.nan
            mad_sigclip_rm2 = np.nan
            mean_sigclip_rm2 = np.nan
            stdev_sigclip_rm2 = np.nan
            ndet_sigclip_rm2 = np.nan

    else:

        median_rm2, mad_rm2, mean_rm2, stdev_rm2 = np.nan, np.nan, np.nan, np.nan
        ndet_rm2 = np.nan
        median_sigclip_rm2, mad_sigclip_rm2 = np.nan, np.nan
        mean_sigclip_rm2, stdev_sigclip_rm2 = np.nan, np.nan
        ndet_sigclip_rm2 = np.nan

    # RM3
    if len(rm3) > 4:

        finiteind = np.isfinite(rm3)
        rm3 = rm3[finiteind]
        median_rm3 = np.median(rm3)
        mad_rm3 = np.median(np.fabs(rm3 - median_rm3))
        mean_rm3 = np.mean(rm3)
        stdev_rm3 = np.std(rm3)
        ndet_rm3 = len(rm3)

        if sigclip:
            sigclip_rm3, lo, hi = stats_sigmaclip(rm3,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_rm3 = np.median(sigclip_rm3)
            mad_sigclip_rm3 = np.median(np.fabs(sigclip_rm3 -
                                                median_sigclip_rm3))
            mean_sigclip_rm3 = np.mean(sigclip_rm3)
            stdev_sigclip_rm3 = np.std(sigclip_rm3)
            ndet_sigclip_rm3 = len(sigclip_rm3)

        else:
            median_sigclip_rm3 = np.nan
            mad_sigclip_rm3 = np.nan
            mean_sigclip_rm3 = np.nan
            stdev_sigclip_rm3 = np.nan
            ndet_sigclip_rm3 = np.nan

    else:

        median_rm3, mad_rm3, mean_rm3, stdev_rm3 = (np.nan, np.nan,
                                                    np.nan, np.nan)
        ndet_rm3 = np.nan
        median_sigclip_rm3, mad_sigclip_rm3 = np.nan, np.nan
        mean_sigclip_rm3, stdev_sigclip_rm3 = np.nan, np.nan
        ndet_sigclip_rm3 = np.nan

    # EP1
    if len(ep1) > 4:

        finiteind = np.isfinite(ep1)
        ep1 = ep1[finiteind]
        median_ep1 = np.median(ep1)
        mad_ep1 = np.median(np.fabs(ep1 - median_ep1))
        mean_ep1 = np.mean(ep1)
        stdev_ep1 = np.std(ep1)
        ndet_ep1 = len(ep1)

        if sigclip:
            sigclip_ep1, lo, hi = stats_sigmaclip(ep1,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep1 = np.median(sigclip_ep1)
            mad_sigclip_ep1 = np.median(np.fabs(sigclip_ep1 -
                                                median_sigclip_ep1))
            mean_sigclip_ep1 = np.mean(sigclip_ep1)
            stdev_sigclip_ep1 = np.std(sigclip_ep1)
            ndet_sigclip_ep1 = len(sigclip_ep1)

        else:
            median_sigclip_ep1 = np.nan
            mad_sigclip_ep1 = np.nan
            mean_sigclip_ep1 = np.nan
            stdev_sigclip_ep1 = np.nan
            ndet_sigclip_ep1 = np.nan

    else:

        median_ep1, mad_ep1, mean_ep1, stdev_ep1 = np.nan, np.nan, np.nan, np.nan
        ndet_ep1 = np.nan
        median_sigclip_ep1, mad_sigclip_ep1 = np.nan, np.nan
        mean_sigclip_ep1, stdev_sigclip_ep1 = np.nan, np.nan
        ndet_sigclip_ep1 = np.nan

    # EP2
    if len(ep2) > 4:

        finiteind = np.isfinite(ep2)
        ep2 = ep2[finiteind]
        median_ep2 = np.median(ep2)
        mad_ep2 = np.median(np.fabs(ep2 - median_ep2))
        mean_ep2 = np.mean(ep2)
        stdev_ep2 = np.std(ep2)
        ndet_ep2 = len(ep2)

        if sigclip:
            sigclip_ep2, lo, hi = stats_sigmaclip(ep2,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep2 = np.median(sigclip_ep2)
            mad_sigclip_ep2 = np.median(np.fabs(sigclip_ep2 -
                                                median_sigclip_ep2))
            mean_sigclip_ep2 = np.mean(sigclip_ep2)
            stdev_sigclip_ep2 = np.std(sigclip_ep2)
            ndet_sigclip_ep2 = len(sigclip_ep2)

        else:
            median_sigclip_ep2 = np.nan
            mad_sigclip_ep2 = np.nan
            mean_sigclip_ep2 = np.nan
            stdev_sigclip_ep2 = np.nan
            ndet_sigclip_ep2 = np.nan

    else:

        median_ep2, mad_ep2, mean_ep2, stdev_ep2 = np.nan, np.nan, np.nan, np.nan
        ndet_ep2 = np.nan
        median_sigclip_ep2, mad_sigclip_ep2 = np.nan, np.nan
        mean_sigclip_ep2, stdev_sigclip_ep2 = np.nan, np.nan
        ndet_sigclip_ep2 = np.nan

    # EP3
    if len(ep3) > 4:

        finiteind = np.isfinite(ep3)
        ep3 = ep3[finiteind]
        median_ep3 = np.median(ep3)
        mad_ep3 = np.median(np.fabs(ep3 - median_ep3))
        mean_ep3 = np.mean(ep3)
        stdev_ep3 = np.std(ep3)
        ndet_ep3 = len(ep3)

        if sigclip:
            sigclip_ep3, lo, hi = stats_sigmaclip(ep3,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep3 = np.median(sigclip_ep3)
            mad_sigclip_ep3 = np.median(np.fabs(sigclip_ep3 -
                                                median_sigclip_ep3))
            mean_sigclip_ep3 = np.mean(sigclip_ep3)
            stdev_sigclip_ep3 = np.std(sigclip_ep3)
            ndet_sigclip_ep3 = len(sigclip_ep3)

        else:
            median_sigclip_ep3 = np.nan
            mad_sigclip_ep3 = np.nan
            mean_sigclip_ep3 = np.nan
            stdev_sigclip_ep3 = np.nan
            ndet_sigclip_ep3 = np.nan

    else:

        median_ep3, mad_ep3, mean_ep3, stdev_ep3 = np.nan, np.nan, np.nan, np.nan
        ndet_ep3 = np.nan
        median_sigclip_ep3, mad_sigclip_ep3 = np.nan, np.nan
        mean_sigclip_ep3, stdev_sigclip_ep3 = np.nan, np.nan
        ndet_sigclip_ep3 = np.nan

    # TF1
    if len(tf1) > 4:

        finiteind = np.isfinite(tf1)
        tf1 = tf1[finiteind]
        median_tf1 = np.median(tf1)
        mad_tf1 = np.median(np.fabs(tf1 - median_tf1))
        mean_tf1 = np.mean(tf1)
        stdev_tf1 = np.std(tf1)
        ndet_tf1 = len(tf1)

        if sigclip:
            sigclip_tf1, lo, hi = stats_sigmaclip(tf1,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf1 = np.median(sigclip_tf1)
            mad_sigclip_tf1 = np.median(np.fabs(sigclip_tf1 -
                                                median_sigclip_tf1))
            mean_sigclip_tf1 = np.mean(sigclip_tf1)
            stdev_sigclip_tf1 = np.std(sigclip_tf1)
            ndet_sigclip_tf1 = len(sigclip_tf1)

        else:
            median_sigclip_tf1 = np.nan
            mad_sigclip_tf1 = np.nan
            mean_sigclip_tf1 = np.nan
            stdev_sigclip_tf1 = np.nan
            ndet_sigclip_tf1 = np.nan

    else:

        median_tf1, mad_tf1, mean_tf1, stdev_tf1 = np.nan, np.nan, np.nan, np.nan
        ndet_tf1 = np.nan
        median_sigclip_tf1, mad_sigclip_tf1 = np.nan, np.nan
        mean_sigclip_tf1, stdev_sigclip_tf1 = np.nan, np.nan
        ndet_sigclip_tf1 = np.nan

    # TF2
    if len(tf2) > 4:

        finiteind = np.isfinite(tf2)
        tf2 = tf2[finiteind]
        median_tf2 = np.median(tf2)
        mad_tf2 = np.median(np.fabs(tf2 - median_tf2))
        mean_tf2 = np.mean(tf2)
        stdev_tf2 = np.std(tf2)
        ndet_tf2 = len(tf2)

        if sigclip:
            sigclip_tf2, lo, hi = stats_sigmaclip(tf2,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf2 = np.median(sigclip_tf2)
            mad_sigclip_tf2 = np.median(np.fabs(sigclip_tf2 -
                                                median_sigclip_tf2))
            mean_sigclip_tf2 = np.mean(sigclip_tf2)
            stdev_sigclip_tf2 = np.std(sigclip_tf2)
            ndet_sigclip_tf2 = len(sigclip_tf2)

        else:
            median_sigclip_tf2 = np.nan
            mad_sigclip_tf2 = np.nan
            mean_sigclip_tf2 = np.nan
            stdev_sigclip_tf2 = np.nan
            ndet_sigclip_tf2 = np.nan

    else:

        median_tf2, mad_tf2, mean_tf2, stdev_tf2 = np.nan, np.nan, np.nan, np.nan
        ndet_tf2 = np.nan
        median_sigclip_tf2, mad_sigclip_tf2 = np.nan, np.nan
        mean_sigclip_tf2, stdev_sigclip_tf2 = np.nan, np.nan
        ndet_sigclip_tf2 = np.nan

    # TF3
    if len(tf3) > 4:

        finiteind = np.isfinite(tf3)
        tf3 = tf3[finiteind]
        median_tf3 = np.median(tf3)
        mad_tf3 = np.median(np.fabs(tf3 - median_tf3))
        mean_tf3 = np.mean(tf3)
        stdev_tf3 = np.std(tf3)
        ndet_tf3 = len(tf3)

        if sigclip:
            sigclip_tf3, lo, hi = stats_sigmaclip(tf3,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf3 = np.median(sigclip_tf3)
            mad_sigclip_tf3 = np.median(np.fabs(sigclip_tf3 -
                                                median_sigclip_tf3))
            mean_sigclip_tf3 = np.mean(sigclip_tf3)
            stdev_sigclip_tf3 = np.std(sigclip_tf3)
            ndet_sigclip_tf3 = len(sigclip_tf3)

        else:
            median_sigclip_tf3 = np.nan
            mad_sigclip_tf3 = np.nan
            mean_sigclip_tf3 = np.nan
            stdev_sigclip_tf3 = np.nan
            ndet_sigclip_tf3 = np.nan

    else:

        median_tf3, mad_tf3, mean_tf3, stdev_tf3 = np.nan, np.nan, np.nan, np.nan
        ndet_tf3 = np.nan
        median_sigclip_tf3, mad_sigclip_tf3 = np.nan, np.nan
        mean_sigclip_tf3, stdev_sigclip_tf3 = np.nan, np.nan
        ndet_sigclip_tf3 = np.nan


    ## COLLECT STATS

    print('%sZ: done with statistics for %s' %
          (datetime.utcnow().isoformat(), lcfile))

    return {'lcfile':lcfile,
            'lcobj':os.path.splitext(os.path.basename(lcfile))[0],
            # reduced mags aperture 1
            'median_rm1':median_rm1,
            'mad_rm1':mad_rm1,
            'mean_rm1':mean_rm1,
            'stdev_rm1':stdev_rm1,
            'ndet_rm1':ndet_rm1,
            'median_sigclip_rm1':median_sigclip_rm1,
            'mad_sigclip_rm1':mad_sigclip_rm1,
            'mean_sigclip_rm1':mean_sigclip_rm1,
            'stdev_sigclip_rm1':stdev_sigclip_rm1,
            'ndet_sigclip_rm1':ndet_sigclip_rm1,
            # reduced mags aperture 2
            'median_rm2':median_rm2,
            'mad_rm2':mad_rm2,
            'mean_rm2':mean_rm2,
            'stdev_rm2':stdev_rm2,
            'ndet_rm2':ndet_rm2,
            'median_sigclip_rm2':median_sigclip_rm2,
            'mad_sigclip_rm2':mad_sigclip_rm2,
            'mean_sigclip_rm2':mean_sigclip_rm2,
            'stdev_sigclip_rm2':stdev_sigclip_rm2,
            'ndet_sigclip_rm2':ndet_sigclip_rm2,
            # reduced mags aperture 3
            'median_rm3':median_rm3,
            'mad_rm3':mad_rm3,
            'mean_rm3':mean_rm3,
            'stdev_rm3':stdev_rm3,
            'ndet_rm3':ndet_rm3,
            'median_sigclip_rm3':median_sigclip_rm3,
            'mad_sigclip_rm3':mad_sigclip_rm3,
            'mean_sigclip_rm3':mean_sigclip_rm3,
            'stdev_sigclip_rm3':stdev_sigclip_rm3,
            'ndet_sigclip_rm3':ndet_sigclip_rm3,
            # EPD mags aperture 1
            'median_ep1':median_ep1,
            'mad_ep1':mad_ep1,
            'mean_ep1':mean_ep1,
            'stdev_ep1':stdev_ep1,
            'ndet_ep1':ndet_ep1,
            'median_sigclip_ep1':median_sigclip_ep1,
            'mad_sigclip_ep1':mad_sigclip_ep1,
            'mean_sigclip_ep1':mean_sigclip_ep1,
            'stdev_sigclip_ep1':stdev_sigclip_ep1,
            'ndet_sigclip_ep1':ndet_sigclip_ep1,
            # EPD mags aperture 2
            'median_ep2':median_ep2,
            'mad_ep2':mad_ep2,
            'mean_ep2':mean_ep2,
            'stdev_ep2':stdev_ep2,
            'ndet_ep2':ndet_ep2,
            'median_sigclip_ep2':median_sigclip_ep2,
            'mad_sigclip_ep2':mad_sigclip_ep2,
            'mean_sigclip_ep2':mean_sigclip_ep2,
            'stdev_sigclip_ep2':stdev_sigclip_ep2,
            'ndet_sigclip_ep2':ndet_sigclip_ep2,
            # EPD mags aperture 3
            'median_ep3':median_ep3,
            'mad_ep3':mad_ep3,
            'mean_ep3':mean_ep3,
            'stdev_ep3':stdev_ep3,
            'ndet_ep3':ndet_ep3,
            'median_sigclip_ep3':median_sigclip_ep3,
            'mad_sigclip_ep3':mad_sigclip_ep3,
            'mean_sigclip_ep3':mean_sigclip_ep3,
            'stdev_sigclip_ep3':stdev_sigclip_ep3,
            'ndet_sigclip_ep3':ndet_sigclip_ep3,
            # TFA mags aperture 1
            'median_tf1':median_tf1,
            'mad_tf1':mad_tf1,
            'mean_tf1':mean_tf1,
            'stdev_tf1':stdev_tf1,
            'ndet_tf1':ndet_tf1,
            'median_sigclip_tf1':median_sigclip_tf1,
            'mad_sigclip_tf1':mad_sigclip_tf1,
            'mean_sigclip_tf1':mean_sigclip_tf1,
            'stdev_sigclip_tf1':stdev_sigclip_tf1,
            'ndet_sigclip_tf1':ndet_sigclip_tf1,
            # TFA mags aperture 2
            'median_tf2':median_tf2,
            'mad_tf2':mad_tf2,
            'mean_tf2':mean_tf2,
            'stdev_tf2':stdev_tf2,
            'ndet_tf2':ndet_tf2,
            'median_sigclip_tf2':median_sigclip_tf2,
            'mad_sigclip_tf2':mad_sigclip_tf2,
            'mean_sigclip_tf2':mean_sigclip_tf2,
            'stdev_sigclip_tf2':stdev_sigclip_tf2,
            'ndet_sigclip_tf2':ndet_sigclip_tf2,
            # TFA mags aperture 3
            'median_tf3':median_tf3,
            'mad_tf3':mad_tf3,
            'mean_tf3':mean_tf3,
            'stdev_tf3':stdev_tf3,
            'ndet_tf3':ndet_tf3,
            'median_sigclip_tf3':median_sigclip_tf3,
            'mad_sigclip_tf3':mad_sigclip_tf3,
            'mean_sigclip_tf3':mean_sigclip_tf3,
            'stdev_sigclip_tf3':stdev_sigclip_tf3,
            'ndet_sigclip_tf3':ndet_sigclip_tf3}



def lc_statistics_worker(task):
    '''
    This is a worker that runs the function above in a parallel worker pool.

    '''

    try:
        return get_lc_statistics(task[0], **task[1])
    except Exception as e:
        print('SOMETHING WENT WRONG! task was %s' % task)
        return None



def parallel_lc_statistics(lcdir,
                           lcglob,
                           fovcatalog,
                           fovcatcols=(0,9), # objectid, magcol to use
                           fovcatmaglabel='r',
                           outfile=None,
                           nworkers=16,
                           workerntasks=500,
                           rmcols=[18,19,20],
                           epcols=[21,22,23],
                           tfcols=[24,25,26],
                           sigclip=4.0):
    '''
    This calculates statistics on all lc files in lcdir using lcglob to find the
    files. Puts the results in text file outfile. Needs the fovcatalog to get
    the catalog magnitude to use as the canonical magnitude.

    outfile contains the following columns:

    object, ndet,
    median RM[1-3], MAD RM[1-3], mean RM[1-3], stdev RM[1-3],
    median EP[1-3], MAD EP[1-3], mean EP[1-3], stdev EP[1-3],
    median TF[1-3], MAD TF[1-3], mean TF[1-3], stdev TF[1-3]

    if a value is missing, it will be np.nan.

    '''

    lcfiles = glob.glob(os.path.join(lcdir, lcglob))

    tasks = [[x, {'rmcols':rmcols,
                  'epcols':epcols,
                  'tfcols':tfcols,
                  'sigclip':sigclip}] for x in lcfiles]

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    results = pool.map(lc_statistics_worker, tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    if not outfile:
        outfile = os.path.join(lcdir, 'lightcurve-statistics.txt')

    outf = open(outfile,'wb')

    outlineformat = (
        '%s %.3f  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s\n'
        )

    outheader = '# total objects: %s, sigmaclip used: %s\n' % (len(lcfiles),
                                                               sigclip)
    outf.write(outheader)

    outcolumnkey = (
        '# columns are:\n'
        '# 0,1: object, catalog mag %s\n'
        '# 2,3,4,5,6: median RM1, MAD RM1, mean RM1, stdev RM1, ndet RM1\n'
        '# 7,8,9,10,11: sigma-clipped median RM1, MAD RM1, mean RM1, stdev RM1, '
        'ndet RM1\n'
        '# 12,13,14,15,16: median RM2, MAD RM2, mean RM2, stdev RM2, ndet RM2\n'
        '# 17,18,19,20,21: sigma-clipped median RM2, MAD RM2, mean RM2, '
        'stdev RM2, ndet RM2\n'
        '# 22,23,24,25,26: median RM3, MAD RM3, mean RM3, stdev RM3, ndet RM3\n'
        '# 27,28,29,30,31: sigma-clipped median RM3, MAD RM3, mean RM3, '
        'stdev RM3, ndet RM3\n'
        '# 32,33,34,35,36: median EP1, MAD EP1, mean EP1, stdev EP1, ndet EP1\n'
        '# 37,38,39,40,41: sigma-clipped median EP1, MAD EP1, mean EP1, '
        'stdev EP1, ndet EP1\n'
        '# 42,43,44,45,46: median EP2, MAD EP2, mean EP2, stdev EP2, ndet EP2\n'
        '# 47,48,49,50,51: sigma-clipped median EP2, MAD EP2, mean EP2, '
        'stdev EP2, ndet EP2\n'
        '# 52,53,54,55,56: median EP3, MAD EP3, mean EP3, stdev EP3, ndet EP3\n'
        '# 57,58,59,60,61: sigma-clipped median EP3, MAD EP3, mean EP3, '
        'stdev EP3, ndet EP3\n'
        '# 62,63,64,65,66: median TF1, MAD TF1, mean TF1, stdev TF1, ndet TF1\n'
        '# 67,68,69,70,71: sigma-clipped median TF1, MAD TF1, mean TF1, '
        'stdev TF1, ndet TF1\n'
        '# 72,73,74,75,76: median TF2, MAD TF2, mean TF2, stdev TF2, ndet TF2\n'
        '# 77,78,79,80,81: sigma-clipped median TF2, MAD TF2, mean TF2, '
        'stdev TF2, ndet TF2\n'
        '# 82,83,84,85,86: median TF3, MAD TF3, mean TF3, stdev TF3, ndet TF3\n'
        '# 87,88,89,90,91: sigma-clipped median TF3, MAD TF3, mean TF3, '
        'stdev TF3, ndet TF3\n'
        ) % fovcatmaglabel
    outf.write(outcolumnkey)

    # open the fovcatalog and read in the column magnitudes and hatids
    fovcat = np.genfromtxt(fovcatalog,
                           usecols=fovcatcols,
                           dtype='S17,f8',
                           names=['objid','mag'])


    for stat in results:

        if stat is not None:

            # find the catalog mag for this object
            try:
                catmag = fovcat['mag'][np.where(fovcat['objid'] == stat['lcobj'])]
                if not catmag:
                    print('no catalog mag for %s, using median TF3 mag' %
                          stat['lcobj'])
                    catmag = stat['median_tf3']
            except Exception as e:
                print('no catalog mag for %s, using median TF3 mag' %
                      stat['lcobj'])
                catmag = stat['median_tf3']

            outline = outlineformat % (
                stat['lcobj'],
                catmag,

                stat['median_rm1'],
                stat['mad_rm1'],
                stat['mean_rm1'],
                stat['stdev_rm1'],
                stat['ndet_rm1'],
                stat['median_sigclip_rm1'],
                stat['mad_sigclip_rm1'],
                stat['mean_sigclip_rm1'],
                stat['stdev_sigclip_rm1'],
                stat['ndet_sigclip_rm1'],

                stat['median_rm2'],
                stat['mad_rm2'],
                stat['mean_rm2'],
                stat['stdev_rm2'],
                stat['ndet_rm2'],
                stat['median_sigclip_rm2'],
                stat['mad_sigclip_rm2'],
                stat['mean_sigclip_rm2'],
                stat['stdev_sigclip_rm2'],
                stat['ndet_sigclip_rm2'],

                stat['median_rm3'],
                stat['mad_rm3'],
                stat['mean_rm3'],
                stat['stdev_rm3'],
                stat['ndet_rm3'],
                stat['median_sigclip_rm3'],
                stat['mad_sigclip_rm3'],
                stat['mean_sigclip_rm3'],
                stat['stdev_sigclip_rm3'],
                stat['ndet_sigclip_rm3'],

                stat['median_ep1'],
                stat['mad_ep1'],
                stat['mean_ep1'],
                stat['stdev_ep1'],
                stat['ndet_ep1'],
                stat['median_sigclip_ep1'],
                stat['mad_sigclip_ep1'],
                stat['mean_sigclip_ep1'],
                stat['stdev_sigclip_ep1'],
                stat['ndet_sigclip_ep1'],

                stat['median_ep2'],
                stat['mad_ep2'],
                stat['mean_ep2'],
                stat['stdev_ep2'],
                stat['ndet_ep2'],
                stat['median_sigclip_ep2'],
                stat['mad_sigclip_ep2'],
                stat['mean_sigclip_ep2'],
                stat['stdev_sigclip_ep2'],
                stat['ndet_sigclip_ep2'],

                stat['median_ep3'],
                stat['mad_ep3'],
                stat['mean_ep3'],
                stat['stdev_ep3'],
                stat['ndet_ep3'],
                stat['median_sigclip_ep3'],
                stat['mad_sigclip_ep3'],
                stat['mean_sigclip_ep3'],
                stat['stdev_sigclip_ep3'],
                stat['ndet_sigclip_ep3'],

                stat['median_tf1'],
                stat['mad_tf1'],
                stat['mean_tf1'],
                stat['stdev_tf1'],
                stat['ndet_tf1'],
                stat['median_sigclip_tf1'],
                stat['mad_sigclip_tf1'],
                stat['mean_sigclip_tf1'],
                stat['stdev_sigclip_tf1'],
                stat['ndet_sigclip_tf1'],

                stat['median_tf2'],
                stat['mad_tf2'],
                stat['mean_tf2'],
                stat['stdev_tf2'],
                stat['ndet_tf2'],
                stat['median_sigclip_tf2'],
                stat['mad_sigclip_tf2'],
                stat['mean_sigclip_tf2'],
                stat['stdev_sigclip_tf2'],
                stat['ndet_sigclip_tf2'],

                stat['median_tf3'],
                stat['mad_tf3'],
                stat['mean_tf3'],
                stat['stdev_tf3'],
                stat['ndet_tf3'],
                stat['median_sigclip_tf3'],
                stat['mad_sigclip_tf3'],
                stat['mean_sigclip_tf3'],
                stat['stdev_sigclip_tf3'],
                stat['ndet_sigclip_tf3']
                )
            outf.write(outline)

    outf.close()

    print('%sZ: wrote statistics to file %s' %
          (datetime.utcnow().isoformat(), outfile))

    return results


def read_stats_file(statsfile):
    '''
    Reads the stats file into a numpy recarray.

    '''

    # open the statfile and read all the columns
    stats = np.genfromtxt(
        statsfile,
        dtype=('S17,f8,'
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM1
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM2
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM3
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP1
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP2
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP3
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF1
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF2
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8'), # TF3
        names=['lcobj','cat_mag',
               'med_rm1','mad_rm1','mean_rm1','stdev_rm1','ndet_rm1',
               'med_sc_rm1','mad_sc_rm1','mean_sc_rm1','stdev_sc_rm1',
               'ndet_sc_rm1',
               'med_rm2','mad_rm2','mean_rm2','stdev_rm2','ndet_rm2',
               'med_sc_rm2','mad_sc_rm2','mean_sc_rm2','stdev_sc_rm2',
               'ndet_sc_rm2',
               'med_rm3','mad_rm3','mean_rm3','stdev_rm3','ndet_rm3',
               'med_sc_rm3','mad_sc_rm3','mean_sc_rm3','stdev_sc_rm3',
               'ndet_sc_rm3',
               'med_ep1','mad_ep1','mean_ep1','stdev_ep1','ndet_ep1',
               'med_sc_ep1','mad_sc_ep1','mean_sc_ep1','stdev_sc_ep1',
               'ndet_sc_ep1',
               'med_ep2','mad_ep2','mean_ep2','stdev_ep2','ndet_ep2',
               'med_sc_ep2','mad_sc_ep2','mean_sc_ep2','stdev_sc_ep2',
               'ndet_sc_ep2',
               'med_ep3','mad_ep3','mean_ep3','stdev_ep3','ndet_ep3',
               'med_sc_ep3','mad_sc_ep3','mean_sc_ep3','stdev_sc_ep3',
               'ndet_sc_ep3',
               'med_tf1','mad_tf1','mean_tf1','stdev_tf1','ndet_tf1',
               'med_sc_tf1','mad_sc_tf1','mean_sc_tf1','stdev_sc_tf1',
               'ndet_sc_tf1',
               'med_tf2','mad_tf2','mean_tf2','stdev_tf2','ndet_tf2',
               'med_sc_tf2','mad_sc_tf2','mean_sc_tf2','stdev_sc_tf2',
               'ndet_sc_tf2',
               'med_tf3','mad_tf3','mean_tf3','stdev_tf3','ndet_tf3',
               'med_sc_tf3','mad_sc_tf3','mean_sc_tf3','stdev_sc_tf3',
               'ndet_sc_tf3']
        )

    return stats


##########################
## LC BINNING FUNCTIONS ##
##########################

def time_bin_lightcurve(lcprefix,
                        lcexts=('epdlc',
                                'tfalc.TF1','tfalc.TF2','tfalc.TF3'),
                        jdcol=0,
                        lcmagcols=([21,22,23],[24,],[24,],[24,]),
                        binsize=540,
                        minperbin=10,
                        outfile=None):
    '''
    This bins a lightcurve in time using the binsize given. binsize is in
    seconds.

    Needs only the lcprefix; figures out the fnames and cols using the lcexts
    and lcmagcols kwargs.

    '''

    collected_binned_mags = {}

    for ext, magcolspec in zip(lcexts, lcmagcols):

        lcfname = '%s.%s' % (lcprefix, ext)

        if not os.path.exists(lcfname):
            print('LC: %s does not exist! skipping...' % lcfname)
            continue

        lcmagcols = ['AP%s' % x for x in range(len(magcolspec))]

        collected_binned_mags[ext] = {x:{} for x in lcmagcols}

        # extract the JD and magnitude columns
        lcdata = np.genfromtxt(lcfname,
                               usecols=tuple([jdcol] + magcolspec),
                               names=['rjd'] + lcmagcols)

        if lcdata['rjd'].shape and len(lcdata['rjd']) > 10:

            # convert binsize in seconds to JD units
            binsizejd = binsize/(86400.0)
            nbins = int(np.ceil((np.max(lcdata['rjd']) -
                             np.min(lcdata['rjd']))/binsizejd) + 1)

            minjd = np.min(lcdata['rjd'])
            jdbins = [(minjd + x*binsizejd) for x in range(nbins)]

            # make a KD-tree on the JDs so we can do fast distance calculations
            rjd_coords = np.array([[x,1.0] for x in lcdata['rjd']])
            jdtree = kdtree(rjd_coords)
            binned_timeseries_indices = []

            for jd in jdbins:
                # find all bin indices close to within binsizejd of this point
                # using the kdtree query. we use the p-norm = 1 (I think this
                # means straight-up pairwise distance? FIXME: check this)
                bin_indices = jdtree.query_ball_point(np.array([jd,1.0]),
                                                      binsizejd/2.0, p=1.0)

                # if the bin_indices have already been collected, then we're
                # done with this bin, move to the next one. if they haven't,
                # then this is the start of a new bin.
                if (bin_indices not in binned_timeseries_indices and
                    len(bin_indices)) > 0:
                    binned_timeseries_indices.append(bin_indices)


            # convert to ndarrays
            binned_timeseries_indices = [np.array(x) for x in
                                         binned_timeseries_indices]

            collected_binned_mags[ext]['jdbins'] = binned_timeseries_indices
            collected_binned_mags[ext]['nbins'] = len(binned_timeseries_indices)

            # collect the JDs
            binned_jd = np.array([np.median(lcdata['rjd'][x])
                                  for x in binned_timeseries_indices])
            collected_binned_mags[ext]['RJD'] = binned_jd

            # now collect the mags
            for magcol in lcmagcols:

                print('%s.%s: binning aperture %s mags...' %
                      (lcprefix,ext,magcol))

                mags = lcdata[magcol]
                collected_binned_mags[ext][magcol] = (
                    np.array([np.median(lcdata[magcol][x])
                              for x in binned_timeseries_indices])
                    )


    collected_binned_mags['binsize'] = binsize

    # write everything to a file. This is a pickled file because we might have
    # different row numbers for each column and I don't feel like handling this
    # right now.
    if not outfile:

        outfile = lcprefix + '.binned-%ssec-lc.pkl' % binsize
        outf = open(outfile,'wb')
        pickle.dump(collected_binned_mags, outf, pickle.HIGHEST_PROTOCOL)
        outf.close()

    return outfile


def serial_bin_lightcurves(lcdir,
                           epdlc_glob='*.epdlc',
                           binsizes=[540,3600],
                           lcexts=('epdlc',
                                   'tfalc.TF1','tfalc.TF2','tfalc.TF3'),
                           jdcol=0,
                           lcmagcols=([21,22,23],[24,],[24,],[24,]),
                           nworkers=16,
                           workerntasks=500):

    epdlcfiles = glob.glob(os.path.join(lcdir, epdlc_glob))

    tasks = [(os.path.splitext(x)[0],
              binsizes,
              {'lcexts':lcexts,
               'jdcol':jdcol,
               'lcmagcols':lcmagcols})
             for x in epdlcfiles]

    print('%sZ: %s objects to process, starting serial binning...' %
          (datetime.utcnow().isoformat(), len(epdlcfiles)))

    for epdlc in epdlcfiles:
        lcprefix = os.path.splitext(epdlc)[0]
        for binsize in binsizes:
            binlc = time_bin_lightcurve(lcprefix,binsize=binsize)


def parallel_lcbinning_worker(task):
    '''
    This calls time_bin_lightcurve above with all binsizes specified in
    task[1]. task[0] contains the lcprefix, and task[3] is a dict containing
    kwargs, which are expanded.

    '''

    try:
        results = [time_bin_lightcurve(task[0],binsize=x,**task[2])
                   for x in task[1]]
    except Exception as e:
        print('something went wrong with %s!' % task)
        results = None

    return results


def parallel_bin_lightcurves(lcdir,
                             epdlc_glob='*.epdlc',
                             binsizes=[540,3600],
                             lcexts=('epdlc',
                                     'tfalc.TF1','tfalc.TF2','tfalc.TF3'),
                             jdcol=0,
                             lcmagcols=([21,22,23],[24,],[24,],[24,]),
                             nworkers=16,
                             workerntasks=500):


    epdlcfiles = glob.glob(os.path.join(lcdir, epdlc_glob))

    tasks = [(os.path.splitext(x)[0],
              binsizes,
              {'lcexts':lcexts,
               'jdcol':jdcol,
               'lcmagcols':lcmagcols})
             for x in epdlcfiles]

    print('%sZ: %s objects to process, starting parallel binning...' %
          (datetime.utcnow().isoformat(), len(epdlcfiles)))

    pool = mp.Pool(nworkers, maxtasksperchild=workerntasks)
    results = pool.map(parallel_lcbinning_worker, tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s LCs processed.' %
          (datetime.utcnow().isoformat(), len(epdlcfiles)))


def read_binned_lc(binnedlc):
    '''
    This reads back the binnedlc pkl file to a dictionary.
    '''
    with open(binnedlc,'rb') as lcf:
        lcdict = pickle.load(lcf)
    return lcdict


def get_binnedlc_statistics(lcfile,
                            sigclip=4.0):
    '''
    This collects stats for the binned LCs in the lcfile.

    '''

    # read in the lcfile
    binnedlc = read_binned_lc(lcfile)

    # read each lckey and get its statistics
    stats_dict = {}

    # get the EPDLC first
    ep1 = binnedlc['epdlc']['AP0']
    ep2 = binnedlc['epdlc']['AP1']
    ep3 = binnedlc['epdlc']['AP2']

    # then get the TFALC
    tf1 = binnedlc['tfalc.TF1']['AP0']
    tf2 = binnedlc['tfalc.TF2']['AP0']
    tf3 = binnedlc['tfalc.TF3']['AP0']

    # get stats for each column
    # EP1
    if len(ep1) > 4:

        finiteind = np.isfinite(ep1)
        ep1 = ep1[finiteind]
        median_ep1 = np.median(ep1)
        mad_ep1 = np.median(np.fabs(ep1 - median_ep1))
        mean_ep1 = np.mean(ep1)
        stdev_ep1 = np.std(ep1)
        ndet_ep1 = len(ep1)

        if sigclip:
            sigclip_ep1, lo, hi = stats_sigmaclip(ep1,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep1 = np.median(sigclip_ep1)
            mad_sigclip_ep1 = np.median(np.fabs(sigclip_ep1 -
                                                median_sigclip_ep1))
            mean_sigclip_ep1 = np.mean(sigclip_ep1)
            stdev_sigclip_ep1 = np.std(sigclip_ep1)
            ndet_sigclip_ep1 = len(sigclip_ep1)

        else:
            median_sigclip_ep1 = np.nan
            mad_sigclip_ep1 = np.nan
            mean_sigclip_ep1 = np.nan
            stdev_sigclip_ep1 = np.nan
            ndet_sigclip_ep1 = np.nan

    else:

        median_ep1, mad_ep1, mean_ep1, stdev_ep1 = np.nan, np.nan, np.nan, np.nan
        ndet_ep1 = np.nan
        median_sigclip_ep1, mad_sigclip_ep1 = np.nan, np.nan
        mean_sigclip_ep1, stdev_sigclip_ep1 = np.nan, np.nan
        ndet_sigclip_ep1 = np.nan

    # EP2
    if len(ep2) > 4:

        finiteind = np.isfinite(ep2)
        ep2 = ep2[finiteind]
        median_ep2 = np.median(ep2)
        mad_ep2 = np.median(np.fabs(ep2 - median_ep2))
        mean_ep2 = np.mean(ep2)
        stdev_ep2 = np.std(ep2)
        ndet_ep2 = len(ep2)

        if sigclip:
            sigclip_ep2, lo, hi = stats_sigmaclip(ep2,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep2 = np.median(sigclip_ep2)
            mad_sigclip_ep2 = np.median(np.fabs(sigclip_ep2 -
                                                median_sigclip_ep2))
            mean_sigclip_ep2 = np.mean(sigclip_ep2)
            stdev_sigclip_ep2 = np.std(sigclip_ep2)
            ndet_sigclip_ep2 = len(sigclip_ep2)

        else:
            median_sigclip_ep2 = np.nan
            mad_sigclip_ep2 = np.nan
            mean_sigclip_ep2 = np.nan
            stdev_sigclip_ep2 = np.nan
            ndet_sigclip_ep2 = np.nan

    else:

        median_ep2, mad_ep2, mean_ep2, stdev_ep2 = np.nan, np.nan, np.nan, np.nan
        ndet_ep2 = np.nan
        median_sigclip_ep2, mad_sigclip_ep2 = np.nan, np.nan
        mean_sigclip_ep2, stdev_sigclip_ep2 = np.nan, np.nan
        ndet_sigclip_ep2 = np.nan

    # EP3
    if len(ep3) > 4:

        finiteind = np.isfinite(ep3)
        ep3 = ep3[finiteind]
        median_ep3 = np.median(ep3)
        mad_ep3 = np.median(np.fabs(ep3 - median_ep3))
        mean_ep3 = np.mean(ep3)
        stdev_ep3 = np.std(ep3)
        ndet_ep3 = len(ep3)

        if sigclip:
            sigclip_ep3, lo, hi = stats_sigmaclip(ep3,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_ep3 = np.median(sigclip_ep3)
            mad_sigclip_ep3 = np.median(np.fabs(sigclip_ep3 -
                                                median_sigclip_ep3))
            mean_sigclip_ep3 = np.mean(sigclip_ep3)
            stdev_sigclip_ep3 = np.std(sigclip_ep3)
            ndet_sigclip_ep3 = len(sigclip_ep3)

        else:
            median_sigclip_ep3 = np.nan
            mad_sigclip_ep3 = np.nan
            mean_sigclip_ep3 = np.nan
            stdev_sigclip_ep3 = np.nan
            ndet_sigclip_ep3 = np.nan

    else:

        median_ep3, mad_ep3, mean_ep3, stdev_ep3 = np.nan, np.nan, np.nan, np.nan
        ndet_ep3 = np.nan
        median_sigclip_ep3, mad_sigclip_ep3 = np.nan, np.nan
        mean_sigclip_ep3, stdev_sigclip_ep3 = np.nan, np.nan
        ndet_sigclip_ep3 = np.nan

    # TF1
    if len(tf1) > 4:

        finiteind = np.isfinite(tf1)
        tf1 = tf1[finiteind]
        median_tf1 = np.median(tf1)
        mad_tf1 = np.median(np.fabs(tf1 - median_tf1))
        mean_tf1 = np.mean(tf1)
        stdev_tf1 = np.std(tf1)
        ndet_tf1 = len(tf1)

        if sigclip:
            sigclip_tf1, lo, hi = stats_sigmaclip(tf1,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf1 = np.median(sigclip_tf1)
            mad_sigclip_tf1 = np.median(np.fabs(sigclip_tf1 -
                                                median_sigclip_tf1))
            mean_sigclip_tf1 = np.mean(sigclip_tf1)
            stdev_sigclip_tf1 = np.std(sigclip_tf1)
            ndet_sigclip_tf1 = len(sigclip_tf1)

        else:
            median_sigclip_tf1 = np.nan
            mad_sigclip_tf1 = np.nan
            mean_sigclip_tf1 = np.nan
            stdev_sigclip_tf1 = np.nan
            ndet_sigclip_tf1 = np.nan

    else:

        median_tf1, mad_tf1, mean_tf1, stdev_tf1 = np.nan, np.nan, np.nan, np.nan
        ndet_tf1 = np.nan
        median_sigclip_tf1, mad_sigclip_tf1 = np.nan, np.nan
        mean_sigclip_tf1, stdev_sigclip_tf1 = np.nan, np.nan
        ndet_sigclip_tf1 = np.nan

    # TF2
    if len(tf2) > 4:

        finiteind = np.isfinite(tf2)
        tf2 = tf2[finiteind]
        median_tf2 = np.median(tf2)
        mad_tf2 = np.median(np.fabs(tf2 - median_tf2))
        mean_tf2 = np.mean(tf2)
        stdev_tf2 = np.std(tf2)
        ndet_tf2 = len(tf2)

        if sigclip:
            sigclip_tf2, lo, hi = stats_sigmaclip(tf2,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf2 = np.median(sigclip_tf2)
            mad_sigclip_tf2 = np.median(np.fabs(sigclip_tf2 -
                                                median_sigclip_tf2))
            mean_sigclip_tf2 = np.mean(sigclip_tf2)
            stdev_sigclip_tf2 = np.std(sigclip_tf2)
            ndet_sigclip_tf2 = len(sigclip_tf2)

        else:
            median_sigclip_tf2 = np.nan
            mad_sigclip_tf2 = np.nan
            mean_sigclip_tf2 = np.nan
            stdev_sigclip_tf2 = np.nan
            ndet_sigclip_tf2 = np.nan

    else:

        median_tf2, mad_tf2, mean_tf2, stdev_tf2 = np.nan, np.nan, np.nan, np.nan
        ndet_tf2 = np.nan
        median_sigclip_tf2, mad_sigclip_tf2 = np.nan, np.nan
        mean_sigclip_tf2, stdev_sigclip_tf2 = np.nan, np.nan
        ndet_sigclip_tf2 = np.nan

    # TF3
    if len(tf3) > 4:

        finiteind = np.isfinite(tf3)
        tf3 = tf3[finiteind]
        median_tf3 = np.median(tf3)
        mad_tf3 = np.median(np.fabs(tf3 - median_tf3))
        mean_tf3 = np.mean(tf3)
        stdev_tf3 = np.std(tf3)
        ndet_tf3 = len(tf3)

        if sigclip:
            sigclip_tf3, lo, hi = stats_sigmaclip(tf3,
                                                  low=sigclip,
                                                  high=sigclip)
            median_sigclip_tf3 = np.median(sigclip_tf3)
            mad_sigclip_tf3 = np.median(np.fabs(sigclip_tf3 -
                                                median_sigclip_tf3))
            mean_sigclip_tf3 = np.mean(sigclip_tf3)
            stdev_sigclip_tf3 = np.std(sigclip_tf3)
            ndet_sigclip_tf3 = len(sigclip_tf3)

        else:
            median_sigclip_tf3 = np.nan
            mad_sigclip_tf3 = np.nan
            mean_sigclip_tf3 = np.nan
            stdev_sigclip_tf3 = np.nan
            ndet_sigclip_tf3 = np.nan

    else:

        median_tf3, mad_tf3, mean_tf3, stdev_tf3 = np.nan, np.nan, np.nan, np.nan
        ndet_tf3 = np.nan
        median_sigclip_tf3, mad_sigclip_tf3 = np.nan, np.nan
        mean_sigclip_tf3, stdev_sigclip_tf3 = np.nan, np.nan
        ndet_sigclip_tf3 = np.nan

    ## COLLECT STATS

    print('%sZ: done with statistics for %s' %
          (datetime.utcnow().isoformat(), lcfile))

    return {'lcfile':lcfile,
            'lcobj':os.path.splitext(os.path.basename(lcfile))[0],
            # EPD mags aperture 1
            'median_ep1':median_ep1,
            'mad_ep1':mad_ep1,
            'mean_ep1':mean_ep1,
            'stdev_ep1':stdev_ep1,
            'ndet_ep1':ndet_ep1,
            'median_sigclip_ep1':median_sigclip_ep1,
            'mad_sigclip_ep1':mad_sigclip_ep1,
            'mean_sigclip_ep1':mean_sigclip_ep1,
            'stdev_sigclip_ep1':stdev_sigclip_ep1,
            'ndet_sigclip_ep1':ndet_sigclip_ep1,
            # EPD mags aperture 2
            'median_ep2':median_ep2,
            'mad_ep2':mad_ep2,
            'mean_ep2':mean_ep2,
            'stdev_ep2':stdev_ep2,
            'ndet_ep2':ndet_ep2,
            'median_sigclip_ep2':median_sigclip_ep2,
            'mad_sigclip_ep2':mad_sigclip_ep2,
            'mean_sigclip_ep2':mean_sigclip_ep2,
            'stdev_sigclip_ep2':stdev_sigclip_ep2,
            'ndet_sigclip_ep2':ndet_sigclip_ep2,
            # EPD mags aperture 3
            'median_ep3':median_ep3,
            'mad_ep3':mad_ep3,
            'mean_ep3':mean_ep3,
            'stdev_ep3':stdev_ep3,
            'ndet_ep3':ndet_ep3,
            'median_sigclip_ep3':median_sigclip_ep3,
            'mad_sigclip_ep3':mad_sigclip_ep3,
            'mean_sigclip_ep3':mean_sigclip_ep3,
            'stdev_sigclip_ep3':stdev_sigclip_ep3,
            'ndet_sigclip_ep3':ndet_sigclip_ep3,
            # TFA mags aperture 1
            'median_tf1':median_tf1,
            'mad_tf1':mad_tf1,
            'mean_tf1':mean_tf1,
            'stdev_tf1':stdev_tf1,
            'ndet_tf1':ndet_tf1,
            'median_sigclip_tf1':median_sigclip_tf1,
            'mad_sigclip_tf1':mad_sigclip_tf1,
            'mean_sigclip_tf1':mean_sigclip_tf1,
            'stdev_sigclip_tf1':stdev_sigclip_tf1,
            'ndet_sigclip_tf1':ndet_sigclip_tf1,
            # TFA mags aperture 2
            'median_tf2':median_tf2,
            'mad_tf2':mad_tf2,
            'mean_tf2':mean_tf2,
            'stdev_tf2':stdev_tf2,
            'ndet_tf2':ndet_tf2,
            'median_sigclip_tf2':median_sigclip_tf2,
            'mad_sigclip_tf2':mad_sigclip_tf2,
            'mean_sigclip_tf2':mean_sigclip_tf2,
            'stdev_sigclip_tf2':stdev_sigclip_tf2,
            'ndet_sigclip_tf2':ndet_sigclip_tf2,
            # TFA mags aperture 3
            'median_tf3':median_tf3,
            'mad_tf3':mad_tf3,
            'mean_tf3':mean_tf3,
            'stdev_tf3':stdev_tf3,
            'ndet_tf3':ndet_tf3,
            'median_sigclip_tf3':median_sigclip_tf3,
            'mad_sigclip_tf3':mad_sigclip_tf3,
            'mean_sigclip_tf3':mean_sigclip_tf3,
            'stdev_sigclip_tf3':stdev_sigclip_tf3,
            'ndet_sigclip_tf3':ndet_sigclip_tf3}


def binnedlc_statistics_worker(task):
    '''
    This is a worker that runs the function above in a parallel worker pool.

    '''

    try:
        return get_binnedlc_statistics(task[0], **task[1])
    except Exception as e:
        print('SOMETHING WENT WRONG! task was %s' % task)
        return None



def parallel_binnedlc_statistics(lcdir,
                                 lcglob,
                                 fovcatalog,
                                 fovcatcols=(0,9), # objectid, magcol to use
                                 fovcatmaglabel='r',
                                 outfile=None,
                                 nworkers=16,
                                 workerntasks=500,
                                 sigclip=4.0):
    '''
    This calculates statistics on all binned lc files in lcdir using lcglob to
    find the files. Puts the results in text file outfile.

    outfile contains the following columns:

    object, ndet,
    median EP[1-3], MAD EP[1-3], mean EP[1-3], stdev EP[1-3],
    median TF[1-3], MAD TF[1-3], mean TF[1-3], stdev TF[1-3]

    if a value is missing, it will be np.nan.

    '''

    lcfiles = glob.glob(os.path.join(lcdir, lcglob))

    tasks = [[x, {'sigclip':sigclip}] for x in lcfiles]

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    results = pool.map(binnedlc_statistics_worker, tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    if not outfile:
        outfile = os.path.join(lcdir, 'binned-lightcurve-statistics.txt')

    outf = open(outfile,'wb')

    outlineformat = (
        '%s %.3f  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s  '
        '%.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f %s\n'
        )

    outheader = '# total objects: %s, sigmaclip used: %s\n' % (len(lcfiles),
                                                               sigclip)
    outf.write(outheader)

    outcolumnkey = (
        '# columns are:\n'
        '# 0,1: object, catalog mag %s\n'
        '# 2,3,4,5,6: median EP1, MAD EP1, mean EP1, stdev EP1, ndet EP1\n'
        '# 7,8,9,10,11: sigma-clipped median EP1, MAD EP1, mean EP1, '
        'stdev EP1, ndet EP1\n'
        '# 12,13,14,15,16: median EP2, MAD EP2, mean EP2, stdev EP2, ndet EP2\n'
        '# 17,18,19,20,21: sigma-clipped median EP2, MAD EP2, mean EP2, '
        'stdev EP2, ndet EP2\n'
        '# 22,23,24,25,26: median EP3, MAD EP3, mean EP3, stdev EP3, ndet EP3\n'
        '# 27,28,29,30,31: sigma-clipped median EP3, MAD EP3, mean EP3, '
        'stdev EP3, ndet EP3\n'
        '# 32,33,34,35,36: median TF1, MAD TF1, mean TF1, stdev TF1, ndet TF1\n'
        '# 37,38,39,40,41: sigma-clipped median TF1, MAD TF1, mean TF1, '
        'stdev TF1, ndet TF1\n'
        '# 42,43,44,45,46: median TF2, MAD TF2, mean TF2, stdev TF2, ndet TF2\n'
        '# 47,48,49,50,51: sigma-clipped median TF2, MAD TF2, mean TF2, '
        'stdev TF2, ndet TF2\n'
        '# 52,53,54,55,56: median TF3, MAD TF3, mean TF3, stdev TF3, ndet TF3\n'
        '# 57,58,59,60,61: sigma-clipped median TF3, MAD TF3, mean TF3, '
        'stdev TF3, ndet TF3\n'
        ) % fovcatmaglabel
    outf.write(outcolumnkey)

    # open the fovcatalog and read in the column magnitudes and hatids
    fovcat = np.genfromtxt(fovcatalog,
                           usecols=fovcatcols,
                           dtype='S17,f8',
                           names=['objid','mag'])

    for stat in results:

        if stat is not None:

            # find the catalog mag for this object
            try:
                catmag = fovcat['mag'][
                    np.where(fovcat['objid'] == stat['lcobj'][:15])
                    ]
                if not catmag:
                    print('no catalog mag for %s, using median TF3 mag' %
                          stat['lcobj'])
                    catmag = stat['median_tf3']
            except Exception as e:
                print('no catalog mag for %s, using median TF3 mag' %
                      stat['lcobj'])
                catmag = stat['median_tf3']

            outline = outlineformat % (
                stat['lcobj'],
                catmag,

                stat['median_ep1'],
                stat['mad_ep1'],
                stat['mean_ep1'],
                stat['stdev_ep1'],
                stat['ndet_ep1'],
                stat['median_sigclip_ep1'],
                stat['mad_sigclip_ep1'],
                stat['mean_sigclip_ep1'],
                stat['stdev_sigclip_ep1'],
                stat['ndet_sigclip_ep1'],

                stat['median_ep2'],
                stat['mad_ep2'],
                stat['mean_ep2'],
                stat['stdev_ep2'],
                stat['ndet_ep2'],
                stat['median_sigclip_ep2'],
                stat['mad_sigclip_ep2'],
                stat['mean_sigclip_ep2'],
                stat['stdev_sigclip_ep2'],
                stat['ndet_sigclip_ep2'],

                stat['median_ep3'],
                stat['mad_ep3'],
                stat['mean_ep3'],
                stat['stdev_ep3'],
                stat['ndet_ep3'],
                stat['median_sigclip_ep3'],
                stat['mad_sigclip_ep3'],
                stat['mean_sigclip_ep3'],
                stat['stdev_sigclip_ep3'],
                stat['ndet_sigclip_ep3'],

                stat['median_tf1'],
                stat['mad_tf1'],
                stat['mean_tf1'],
                stat['stdev_tf1'],
                stat['ndet_tf1'],
                stat['median_sigclip_tf1'],
                stat['mad_sigclip_tf1'],
                stat['mean_sigclip_tf1'],
                stat['stdev_sigclip_tf1'],
                stat['ndet_sigclip_tf1'],

                stat['median_tf2'],
                stat['mad_tf2'],
                stat['mean_tf2'],
                stat['stdev_tf2'],
                stat['ndet_tf2'],
                stat['median_sigclip_tf2'],
                stat['mad_sigclip_tf2'],
                stat['mean_sigclip_tf2'],
                stat['stdev_sigclip_tf2'],
                stat['ndet_sigclip_tf2'],

                stat['median_tf3'],
                stat['mad_tf3'],
                stat['mean_tf3'],
                stat['stdev_tf3'],
                stat['ndet_tf3'],
                stat['median_sigclip_tf3'],
                stat['mad_sigclip_tf3'],
                stat['mean_sigclip_tf3'],
                stat['stdev_sigclip_tf3'],
                stat['ndet_sigclip_tf3']
                )
            outf.write(outline)

    outf.close()

    print('%sZ: wrote statistics to file %s' %
          (datetime.utcnow().isoformat(), outfile))

    return results


def read_binnedlc_stats_file(statsfile):
    '''
    Reads the stats file into a numpy recarray.

    '''

    # open the statfile and read all the columns
    stats = np.genfromtxt(
        statsfile,
        dtype=('S17,f8,'
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP1
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP2
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP3
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF1
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF2
               'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8'), # TF3
        names=['lcobj','cat_mag',
               'med_ep1','mad_ep1','mean_ep1','stdev_ep1','ndet_ep1',
               'med_sc_ep1','mad_sc_ep1','mean_sc_ep1','stdev_sc_ep1',
               'ndet_sc_ep1',
               'med_ep2','mad_ep2','mean_ep2','stdev_ep2','ndet_ep2',
               'med_sc_ep2','mad_sc_ep2','mean_sc_ep2','stdev_sc_ep2',
               'ndet_sc_ep2',
               'med_ep3','mad_ep3','mean_ep3','stdev_ep3','ndet_ep3',
               'med_sc_ep3','mad_sc_ep3','mean_sc_ep3','stdev_sc_ep3',
               'ndet_sc_ep3',
               'med_tf1','mad_tf1','mean_tf1','stdev_tf1','ndet_tf1',
               'med_sc_tf1','mad_sc_tf1','mean_sc_tf1','stdev_sc_tf1',
               'ndet_sc_tf1',
               'med_tf2','mad_tf2','mean_tf2','stdev_tf2','ndet_tf2',
               'med_sc_tf2','mad_sc_tf2','mean_sc_tf2','stdev_sc_tf2',
               'ndet_sc_tf2',
               'med_tf3','mad_tf3','mean_tf3','stdev_tf3','ndet_tf3',
               'med_sc_tf3','mad_sc_tf3','mean_sc_tf3','stdev_sc_tf3',
               'ndet_sc_tf3']
        )

    return stats



########################
## PLOTTING FUNCTIONS ##
########################

# lists all plots to make plus columns to use and metadata
STATS_PLOTS = {
    'median-rm1-vs-mad-rm1':{
        'xcol':'med_rm1',
        'ycol':'mad_rm1',
        'title':'RM1 median mag vs. RM1 median abs. dev.',
        'xlabel':'RM1 median magnitude',
        'ylabel':'RM1 median abs. dev.',
        'binned':False
        },
    'median-rm1-vs-mad-rm1-sigclipped':{
        'xcol':'med_sc_rm1',
        'ycol':'mad_sc_rm1',
        'title':'RM1 median mag vs. RM1 median abs. dev. (sigclip LCs)',
        'xlabel':'RM1 median magnitude',
        'ylabel':'RM1 median abs. dev.',
        'binned':False
        },
    'median-rm2-vs-mad-rm2':{
        'xcol':'med_rm2',
        'ycol':'mad_rm2',
        'title':'RM2 median mag vs. RM2 median abs. dev.',
        'xlabel':'RM2 median magnitude',
        'ylabel':'RM2 median abs. dev.',
        'binned':False
        },
    'median-rm2-vs-mad-rm2-sigclipped':{
        'xcol':'med_sc_rm2',
        'ycol':'mad_sc_rm2',
        'title':'RM2 median mag vs. RM2 median abs. dev. (sigclip LCs)',
        'xlabel':'RM2 median magnitude',
        'ylabel':'RM2 median abs. dev.',
        'binned':False
        },
    'median-rm3-vs-mad-rm3':{
        'xcol':'med_rm3',
        'ycol':'mad_rm3',
        'title':'RM3 median mag vs. RM3 median abs. dev.',
        'xlabel':'RM3 median magnitude',
        'ylabel':'RM3 median abs. dev.',
        'binned':False
        },
    'median-rm3-vs-mad-rm3-sigclipped':{
        'xcol':'med_sc_rm3',
        'ycol':'mad_sc_rm3',
        'title':'RM3 median mag vs. RM3 median abs. dev. (sigclip LCs)',
        'xlabel':'RM3 median magnitude',
        'ylabel':'RM3 median abs. dev.',
        'binned':False
        },
    'median-EP1-vs-mad-EP1':{
        'xcol':'med_ep1',
        'ycol':'mad_ep1',
        'title':'EP1 median mag vs. EP1 median abs. dev.',
        'xlabel':'EP1 median magnitude',
        'ylabel':'EP1 median abs. dev.',
        'binned':True
        },
    'median-EP1-vs-mad-EP1-sigclipped':{
        'xcol':'med_sc_ep1',
        'ycol':'mad_sc_ep1',
        'title':'EP1 median mag vs. EP1 median abs. dev. (sigclip LCs)',
        'xlabel':'EP1 median magnitude',
        'ylabel':'EP1 median abs. dev.',
        'binned':True
        },
    'median-EP2-vs-mad-EP2':{
        'xcol':'med_ep2',
        'ycol':'mad_ep2',
        'title':'EP2 median mag vs. EP2 median abs. dev.',
        'xlabel':'EP2 median magnitude',
        'ylabel':'EP2 median abs. dev.',
        'binned':True
        },
    'median-EP2-vs-mad-EP2-sigclipped':{
        'xcol':'med_sc_ep2',
        'ycol':'mad_sc_ep2',
        'title':'EP2 median mag vs. EP2 median abs. dev. (sigclip LCs)',
        'xlabel':'EP2 median magnitude',
        'ylabel':'EP2 median abs. dev.',
        'binned':True
        },
    'median-EP3-vs-mad-EP3':{
        'xcol':'med_ep3',
        'ycol':'mad_ep3',
        'title':'EP3 median mag vs. EP3 median abs. dev.',
        'xlabel':'EP3 median magnitude',
        'ylabel':'EP3 median abs. dev.',
        'binned':True
        },
    'median-EP3-vs-mad-EP3-sigclipped':{
        'xcol':'med_sc_ep3',
        'ycol':'mad_sc_ep3',
        'title':'EP3 median mag vs. EP3 median abs. dev. (sigclip LCs)',
        'xlabel':'EP3 median magnitude',
        'ylabel':'EP3 median abs. dev.',
        'binned':True
        },
    'median-TF1-vs-mad-TF1':{
        'xcol':'med_tf1',
        'ycol':'mad_tf1',
        'title':'TF1 median mag vs. TF1 median abs. dev.',
        'xlabel':'TF1 median magnitude',
        'ylabel':'TF1 median abs. dev.',
        'binned':True
        },
    'median-TF1-vs-mad-TF1-sigclipped':{
        'xcol':'med_sc_tf1',
        'ycol':'mad_sc_tf1',
        'title':'TF1 median mag vs. TF1 median abs. dev. (sigclip LCs)',
        'xlabel':'TF1 median magnitude',
        'ylabel':'TF1 median abs. dev.',
        'binned':True
        },
    'median-TF2-vs-mad-TF2':{
        'xcol':'med_tf2',
        'ycol':'mad_tf2',
        'title':'TF2 median mag vs. TF2 median abs. dev.',
        'xlabel':'TF2 median magnitude',
        'ylabel':'TF2 median abs. dev.',
        'binned':True
        },
    'median-TF2-vs-mad-TF2-sigclipped':{
        'xcol':'med_sc_tf2',
        'ycol':'mad_sc_tf2',
        'title':'TF2 median mag vs. TF2 median abs. dev. (sigclip LCs)',
        'xlabel':'TF2 median magnitude',
        'ylabel':'TF2 median abs. dev.',
        'binned':True
        },
    'median-TF3-vs-mad-TF3':{
        'xcol':'med_tf3',
        'ycol':'mad_tf3',
        'title':'TF3 median mag vs. TF3 median abs. dev.',
        'xlabel':'TF3 median magnitude',
        'ylabel':'TF3 median abs. dev.',
        'binned':True
        },
    'median-TF3-vs-mad-TF3-sigclipped':{
        'xcol':'med_sc_tf3',
        'ycol':'mad_sc_tf3',
        'title':'TF3 median mag vs. TF3 median abs. dev. (sigclip LCs)',
        'xlabel':'TF3 median magnitude',
        'ylabel':'TF3 median abs. dev.',
        'binned':True
        },
    'catalog-r-mag-vs-mad-TF3':{
        'xcol':'cat_mag',
        'ycol':'mad_tf3',
        'title':'FOV catalog SDSS r mag vs. TF3 median abs. dev.',
        'xlabel':'FOV catalog SDSS r mag',
        'ylabel':'TF3 median abs. dev.',
        'binned':True
        },
    'catalog-r-mag-vs-mad-TF3-sigclipped':{
        'xcol':'cat_mag',
        'ycol':'mad_sc_tf3',
        'title':'FOV catalog SDSS r mag vs. TF3 median abs. dev. (sigclip LCs)',
        'xlabel':'FOV catalog SDSS r mag',
        'ylabel':'TF3 median abs. dev.',
        'binned':True
        },
    }


def plot_stats_file(statsfile, outdir, outprefix,
                    binned=False,
                    logy=False,
                    logx=False,
                    rangex=(5.9,14.1)):
    '''
    This makes all the plots for a stats file.

    '''

    # read the stats file
    if binned:
        stats = read_binnedlc_stats_file(statsfile)
    else:
        stats = read_stats_file(statsfile)

    for plot in STATS_PLOTS:

        if binned and not STATS_PLOTS[plot]['binned']:
            print('not making %s for binned LCs' % plot)
            continue

        if logy:
            logypostfix = '-logy'
        else:
            logypostfix = ''
        if logx:
            logxpostfix = '-logx'
        else:
            logxpostfix = ''

        outfile = os.path.join(outdir, '%s-%s%s%s.png' % (outprefix,
                                                          plot,
                                                          logxpostfix,
                                                          logypostfix))

        xcol, ycol = (stats[STATS_PLOTS[plot]['xcol']],
                      stats[STATS_PLOTS[plot]['ycol']])
        xlabel, ylabel = (STATS_PLOTS[plot]['xlabel'],
                          STATS_PLOTS[plot]['ylabel'])
        title = '%s - %s - %s' % (outprefix,
                                  len(xcol),
                                  STATS_PLOTS[plot]['title'])

        # make the plot
        if logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.yscale('log',basex=10.0)
        elif logx:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log',basex=10.0)
        elif logx and logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log',basex=10.0)
            plt.yscale('log',basex=10.0)
        else:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')


        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xlim(rangex)

        if binned:
            plt.ylim((0.0001,1.0))
            plt.hlines([0.0005,0.001,0.002,0.003,0.01],
                       xmin=5.0,xmax=15.0,colors='b')

        else:
            # make the horizontal lines for 10, 5, 1 mmag
            plt.ylim((0.001,1.0))
            plt.hlines([0.001, 0.002, 0.003, 0.004, 0.005, 0.01],
                       xmin=5.0,xmax=15.0,colors='b')

        plt.savefig(outfile)
        plt.close()

        print('%sZ: made %s plot: %s' %
              (datetime.utcnow().isoformat(), title, outfile))



def plot_magrms_comparison(reference_stats_file,
                           comparison_stats_file,
                           ref_name, comp_name,
                           outfile,
                           ref_col='mad_tf3',
                           comp_col='mad_tf3',
                           logy=False, logx=False,
                           rangex=(5.9,14.1)):
    '''
    This makes magrms comparison plots using two lightcurve stats files.

    '''

    ref_stats = read_stats_file(reference_stats_file)
    comp_stats = read_stats_file(comparison_stats_file)

    ref_objects = ref_stats['lcobj']
    comp_objects = comp_stats['lcobj']

    common_objects = np.intersect1d(ref_objects, comp_objects)

    print('common objects = %s' % len(common_objects))

    if len(common_objects) > 0:

        # put together the data for the common objects
        ref_tf3_mag = [ref_stats['cat_mag'][ref_stats['lcobj'] == x]
                       for x in common_objects]
        ref_tf3_compcol = [ref_stats[ref_col][ref_stats['lcobj'] == x]
                           for x in common_objects]

        comp_tf3_mag = [comp_stats['cat_mag'][comp_stats['lcobj'] == x]
                        for x in common_objects]
        comp_tf3_compcol = [comp_stats[comp_col][comp_stats['lcobj'] == x]
                            for x in common_objects]

        tf3_compcol_ratios = (
            np.array(ref_tf3_compcol)/np.array(comp_tf3_compcol)
            )

        xcol, ycol = ref_tf3_mag, tf3_compcol_ratios

        xlabel, ylabel = ('FOV catalog SDSS r mag',
                          'TF3 median abs. dev. %s/%s' % (ref_name, comp_name))
        title = 'comparison of TF3 median abs. dev. - %s/%s' % (ref_name,
                                                                comp_name)

        # make the plot
        if logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.yscale('log')
        elif logx:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log')
        elif logx and logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log')
            plt.yscale('log')
        else:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')


        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xlim(rangex)
        plt.ylim(-0.5,4)

        # make the horizontal lines for 10, 5, 1 mmag
        plt.hlines([-2.0,-1.5,0.5,1,1.5,2.0],
                   xmin=5.0,xmax=15.0,colors='r')

        plt.savefig(outfile)
        plt.close()

    else:

        print('no common objects to use for comparison!')


def plot_ismphot_comparison(apphot_stats_file,
                            ismphot_stats_file,
                            ref_name, comp_name,
                            outfile,
                            comp_cols=(0,1,13),
                            logy=False, logx=False,
                            rangex=(5.9,14.1)):

    ref_stats = read_stats_file(apphot_stats_file)
    comp_stats = np.genfromtxt(ismphot_stats_file,
                               usecols=comp_cols,
                               dtype='S17,f8,f8',
                               names=['lcobj','cat_mag','mad_tf3'])

    ref_objects = ref_stats['lcobj']
    comp_objects = comp_stats['lcobj']

    common_objects = np.intersect1d(ref_objects, comp_objects)

    print('common objects = %s' % len(common_objects))

    if len(common_objects) > 0:

        # put together the data for the common objects
        ref_tf3_mag = [ref_stats['cat_mag'][ref_stats['lcobj'] == x]
                       for x in common_objects]
        ref_tf3_compcol = [ref_stats['mad_tf3'][ref_stats['lcobj'] == x]
                           for x in common_objects]

        comp_tf3_mag = [comp_stats['cat_mag'][comp_stats['lcobj'] == x]
                        for x in common_objects]
        comp_tf3_compcol = [comp_stats['mad_tf3'][comp_stats['lcobj'] == x]
                            for x in common_objects]

        tf3_compcol_ratios = np.array(comp_tf3_compcol)/np.array(ref_tf3_compcol)

        xcol, ycol = ref_tf3_mag, tf3_compcol_ratios

        xlabel, ylabel = ('FOV catalog SDSS r mag',
                          'TF3 median abs. dev. %s/%s' % (comp_name, ref_name))
        title = 'comparison of TF3 median abs. dev. - %s/%s' % (comp_name,
                                                                ref_name)

        # make the plot
        if logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.yscale('log')
        elif logx:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log')
        elif logx and logy:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')
            plt.xscale('log')
            plt.yscale('log')
        else:
            plt.scatter(xcol, ycol,
                        s=1,
                        marker='.')


        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xlim(rangex)
        plt.ylim(-0.5,4)

        # make the horizontal lines for 10, 5, 1 mmag
        plt.hlines([-2.0,-1.5,0.5,1,1.5,2.0],
                   xmin=5.0,xmax=15.0,colors='r')

        plt.savefig(outfile)
        plt.close()

    else:

        print('no common objects to use for comparison!')

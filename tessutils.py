'''
functions specific to wrangling TESS data
'''

##########################################

import aperturephot as ap
import imageutils as iu
import shared_variables as sv

import os
import numpy as np
from glob import glob

from astropy.io import fits
from astroquery.mast import Catalogs

from numpy import array as nparr, all as npall

from datetime import datetime
import multiprocessing as mp

##########################################

DEBUG = False

SATURATIONMASKCMD = ('fiign {fitsfile} -o {outfitsfile} '
                     '-s {saturationlevel} --an')

DQUALITYMASKCMD = ('fiign {fitsfile} -o {outfitsfile} '
                   '-s {minimumreadout} --ignore-nonpositive')

def mask_saturated_stars_worker(task):
    """
    Mask the saturated stars using the calibrated images. Overwrites image with
    fitsh-aware mask in the header. Note the mask is not actually applied to
    the pixel values until we need to visualize; the mask is metadata.

    At first thought, it would be better to do this with raw images because
    once the flat-field is applied, saturated stars do not have a constant
    value.  However inspection of the raw images shows that the saturated stars
    do not have constant values there either. So it's easier to first trim, and
    then apply this.
    """

    fitsname, saturationlevel = task

    cmdtorun = SATURATIONMASKCMD.format(fitsfile=fitsname, outfitsfile=fitsname,
                                        saturationlevel=saturationlevel)

    returncode = os.system(cmdtorun)

    if DEBUG:
        print(cmdtorun)

    if returncode == 0:
        print('%sZ: appended saturated stars mask %s' %
              (datetime.utcnow().isoformat(), fitsname))
        return fitsname
    else:
        print('ERR! %sZ: saturated star mask construction failed for %s' %
              (datetime.utcnow().isoformat(), fitsname))
        return None


def parallel_mask_saturated_stars(fitslist, saturationlevel=65535, nworkers=16,
                                  maxworkertasks=1000):
    '''
    Append mask of saturated stars to fits header using the calibrated TESS
    images, by applying a constant saturation level.

    Kwargs:
        saturationlevel (int): 2**16 - 1 = 65535 is a common choice. This is
        the number that is used to mask out cores of saturated stars.
    '''

    print('%sZ: %s files to mask saturated stars in' %
          (datetime.utcnow().isoformat(), len(fitslist)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, saturationlevel) for x in fitslist]

    # fire up the pool of workers
    results = pool.map(mask_saturated_stars_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {result for result in results}


def mask_dquality_flag_frame(task):
    """
    Mask the entire frame if dquality flag from SPOC is raised. Overwrites
    image with fitsh-aware mask in the header. Note the mask is not
    actually applied to the pixel values until we need to visualize; the
    mask is metadata.
    """

    fitsname, flagvalue = task

    # open the header, check if the quality flag matches the passed value.
    data, hdr = iu.read_fits(fitsname, ext=0)

    # if it matches, mask out everything. otherwise, do nothing.
    if hdr['DQUALITY'] == flagvalue:

        # if you pass `fiign` "0" as the saturation value, you don't get a
        # correct mask. so we mask out non-positive values, and set the
        # saturation value below the lowest positive value.
        minabsdata = np.min(np.abs(data))
        minimumreadout = minabsdata - minabsdata/2

        cmdtorun = DQUALITYMASKCMD.format(fitsfile=fitsname, outfitsfile=fitsname,
                                          minimumreadout=minimumreadout)

        returncode = os.system(cmdtorun)

        if returncode == 0:
            print('%sZ: appended DQUALITY=%d mask to %s' %
                  (datetime.utcnow().isoformat(), flagvalue, fitsname))
            return fitsname
        else:
            print('ERR! %sZ: DQUALITY mask construction failed for %s' %
                  (datetime.utcnow().isoformat(), fitsname))
            return None

    else:
        return 1



def parallel_mask_dquality_flag_frames(fitslist, flagvalue=32,
                                       nworkers=16, maxworkertasks=1000):
    '''
    Append mask to fits header of an ENTIRE CALIBRATED FRAME if it matches the
    passed data quality flag value.

    See e.g., EXP-TESS-ARC-ICD-TM-0014 for a description of the flag values.

    Kwargs:
        flagvalue (int): integer that tells you something about the data
        quality. E.g., "32" is a "reaction wheel desaturation event", AKA a
        "momentum dump".
    '''

    print('%sZ: %s files to check for DQUALITY=%d in' %
          (datetime.utcnow().isoformat(), len(fitslist), flagvalue))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, flagvalue) for x in fitslist]

    # fire up the pool of workers
    results = pool.map(mask_dquality_flag_frame, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    print('%sZ: finished checking for DQUALITY=%d' %
          (datetime.utcnow().isoformat(), flagvalue))

    return {result for result in results}


def parallel_trim_get_single_extension(fitslist, outdir, projid, nworkers=16,
                                       maxworkertasks=1000):
    '''
    see docstring for from_CAL_to_fitsh_compatible
    '''

    print('%sZ: %s files to trim and get single extension' %
          (datetime.utcnow().isoformat(), len(fitslist)))

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    path_exists = []
    for ix, fitsname in enumerate(fitslist):
        outname = (
            outdir + os.path.basename(
                fitsname.replace('-s_ffic.fits', '_cal_img.fits'))
        )
        if os.path.exists(outname):
            path_exists.append(1)
        else:
            path_exists.append(0)

    path_exists = nparr(path_exists)

    if npall(path_exists):
        print(
            'found all {:d} fitsh-compatible files, continuing'.
            format(len(path_exists))
        )
        return 0

    else:
        tasks = [(x, outdir, projid) for x in fitslist]

        # fire up the pool of workers
        results = pool.map(from_CAL_to_fitsh_compatible, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()

        return 1


def from_CAL_to_fitsh_compatible(task):
    '''
    This function takes a calibrated image from MAST and turns it into a
    fitsh-compatible image. This is done BEFORE constructing a mask for the
    image using tessutils.mask_saturated_stars

    This means:
        * get single extension (omit the uncertainty map, for now).
        * trim to remove virtual columns
        * append "PROJID" header keywork.

    Arg:
        task: tuple of (fitspath, outdir), where `fitspath` is a calibrated
        full frame image from MAST. This means its filename should end with
        "-s_ffic.fits".  (Each image is a single CCD, with 4 subarrays).
        outdir (str): directory to write the files to.

    Returns:
        nothing.
    '''
    fitsname, outdir, projid = task

    if fitsname.split('-')[-1] != 's_ffic.fits':
        raise AssertionError('expected calibrated FFI from MAST.')

    outname = (
        outdir + os.path.basename(
            fitsname.replace('-s_ffic.fits', '_cal_img.fits'))
    )

    data, hdr = iu.read_fits(fitsname, ext=1)

    # FITS header is 1-based counting, but python is 0-based. To convert
    # from FITS->python from slicing you need to subtract 1.
    rows = hdr['SCIROWS']-1 # row start.
    rowe = hdr['SCIROWE']-1+1 # row end (inclusive)
    cols = hdr['SCCSA']-1 # col start
    cole = hdr['SCCED']-1+1 # col end (inclusive)

    trim = data[slice(rows,rowe),slice(cols,cole)]

    hdr['PROJID'] = projid

    assert trim.shape == (2048, 2048)

    fits.writeto(outname, trim, header=hdr)

    print('Wrote {:s} to {:s}'.format(fitsname, outname))


def from_ete6_to_fitsh_compatible(fitslist, outdir, projid=42):
    '''
    This function takes a list of ETE6 reduced images and turns them into
    fitsh-compatible images.  Using it is a bad idea, because you shouldn't be
    using ETE6 any more anyway.
    '''

    for fitsname in fitslist:
        from_CAL_to_fitsh_compatible((fitsname, outdir, projid))

##############################
# UNDER CONSTRUCTION
##############################

def read_tess_lightcurve(
    lcfile,
    catfile='/home/lbouma/proj/ete6/data/RED_1-1/2MASS-RA250.546600342-DEC32.4748344421-SIZE24.reformed_catalog',
    infokeys=['id', 'ra', 'dec', 'xi', 'eta', '2massJ', '2massK', '2massqlt',
              '2massI', '2massr', '2massi', '2massz']
    ):
    # 00 rjd    Reduced Julian Date (RJD = JD - 2400000.0)
    # 01 hat    HAT ID of the object
    # 02 rstfc  Unique frame key ({STID}-{FRAMENUMBER}_{CCDNUM})
    # 03 xcc    original X coordinate on CCD
    # 04 ycc    original y coordinate on CCD
    # 05 bgv    Background value
    # 06 bge    Background measurement error
    # 07 fsv    Measured S value
    # 08 fdv    Measured D value
    # 09 fkv    Measured K value
    # 10 im1    Instrumental magnitude in aperture 1
    # 11 ie1    Instrumental magnitude error for aperture 1
    # 12 iq1    Instrumental magnitude quality flag for aperture 1 (0/G OK, X bad)
    # 13 im2    Instrumental magnitude in aperture 2
    # 14 ie2    Instrumental magnitude error for aperture 2
    # 15 iq2    Instrumental magnitude quality flag for aperture 2 (0/G OK, X bad)
    # 16 im3    Instrumental magnitude in aperture 3
    # 17 ie3    Instrumental magnitude error for aperture 3
    # 18 iq3    Instrumental magnitude quality flag for aperture 3 (0/G OK, X bad)
    # 19 rm1    Reduced Mags from magfit in aperture 1
    # 20 rm2    Reduced Mags from magfit in aperture 2
    # 21 rm3    Reduced Mags from magfit in aperture 3

    # read the LC into a numpy recarray
    recarr = np.genfromtxt(lcfile,
                           usecols=(0,3,4,10,11,13,14,16,17),
                           names=['rjd','xcc','ycc','im1','ie1','im2','ie2','im3','ie3'],
                           dtype='f8,f8,f8,f8,f8,f8,f8,f8,f8'
                          )

    # generate the objectid
    # here, we're basically stripping the filename to form the objectid
    objectid = os.path.basename(lcfile).rstrip('.rlc')

    catalog_recarray = read_object_catalog(catfile)

    # look up this object in the catalog
    objind = catalog_recarray['id'] == objectid

    if objind.size > 0:

        objectinfo = {x:np.asscalar(catalog_recarray[x][objind]) for x in infokeys}

    else:

        objectinfo = {'objectid': objectid}

    # this is the lcdict we need
    lcdict = {'objectid': objectid,
              'objectinfo':objectinfo,
              'rjd':recarr['rjd'],
              'xcc':recarr['xcc'],
              'ycc':recarr['ycc'],
              'im1':recarr['im1'],
              'ie1':recarr['ie1'],
              'im2':recarr['im2'],
              'ie2':recarr['ie2'],
              'im3':recarr['im3'],
              'ie3':recarr['ie3']
             }

    return lcdict


def read_object_catalog(catalogfile):
    # you often need an objectid, ra, dec, and a magnitude.
    # get these from the 2mass catalog.

    columns='id,ra,dec,xi,eta,2massJ,2massK,2massqlt,2massI,2massr,2massi,2massz'
    columns = columns.split(',')

    catarr = np.genfromtxt(catalogfile,
                           comments='#',
                           usecols=list(range(len(columns))),
                           names=columns,
                           dtype='U15,f8,f8,f8,f8,f8,f8,U3,f8,f8,f8,f8',
                           delimiter=' '
                           )

    return catarr




'''
functions specific to wrangling TESS data.

contents:

mask_saturated_stars_worker
parallel_mask_saturated_stars
mask_dquality_flag_frame
parallel_mask_dquality_flag_frames
parallel_trim_get_single_extension
from_CAL_to_fitsh_compatible
    (deprecated) from_ete6_to_fitsh_compatible

read_tess_lightcurve
read_object_catalog
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

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii

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


def are_known_HJs_in_field(ra_center, dec_center, outname):
    """
    Given a field center, find which known hot Jupiters are on chip.
    Dependencies: tessmaps, astropy, astroquery

    Args:
        ra_center, dec_center (floats): coordinate of field center in degrees

        outname (str): path to which csv file with details of HJs on chip will
        be written

    Returns:
        True if any HJs on chip, else False.
    """

    from tessmaps import get_time_on_silicon as gts
    from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

    eatab = NasaExoplanetArchive.get_confirmed_planets_table()

    cam_dirn = SkyCoord(ra_center*u.deg, dec_center*u.deg, frame='icrs')
    cam_tuple = (cam_dirn.barycentrictrueecliptic.lat.value,
                 cam_dirn.barycentrictrueecliptic.lon.value)

    is_jup = (eatab['pl_radj'] > 0.4*u.Rjup)
    is_hot = (eatab['pl_orbper'] < 10*u.day)

    is_hj = is_jup & is_hot

    hj_coords = eatab[is_hj]['sky_coord']

    hj_onchip = gts.given_one_camera_get_stars_on_silicon(
        hj_coords, cam_tuple, withgaps=False)

    if np.any(hj_onchip):
        desiredcols = ['pl_hostname', 'pl_letter', 'pl_name', 'pl_orbper',
                       'pl_radj', 'st_optmag', 'gaia_gmag','sky_coord']

        outtab = eatab[is_hj][desiredcols]
        outtab = outtab[hj_onchip.astype(bool)]

        ascii.write(outtab, output=outname, format='ecsv',
                    overwrite=True)
        print('wrote {}'.format(outname))

        return True

    else:
        return False


def measure_known_HJ_SNR(hjonchippath, projcatalogpath, lcdirectory,
                         minxmatchsep=3):
    """
    Args:

        hjonchippath (str): path to csv that tells you which HJs are expected
        to be on (or near) chip.

        projcatalogpath (str): path to projected catalog (from e.g., 2MASS or
        GAIA).

        lcdirectory (str): path to directory with lightcurves. Globs are
        constructed to search for LCs of the HJ matches.

        minxmatchsep (float): minimum separation, in arcseconds, for
        crossmatching of hot Jupiter exoarchive catalog positions to 2mass
        projected plate positions.
    """

    minxmatchsep = minxmatchsep*u.arcsec

    # read HJ information table from `tessutils.are_known_HJs_in_field`
    tab = ascii.read(hjonchippath)

    # if it's not HAT/WASP/KELT, it's probably not a HJ that this tuning test
    # cares about. if so, scrap it.
    wantstrs = ['HAT','WASP','KELT']
    is_wanted = []
    for hn in nparr(tab['pl_hostname']):
        want = False
        for wantstr in wantstrs:
            if wantstr in hn:
                want = True
                break
            else:
                pass
        is_wanted.append(want)
    is_wanted = nparr(is_wanted)

    tab = tab[is_wanted]

    # get the HATIDs that correspond to the HJs on-chip
    projcat = read_object_catalog(projcatalogpath)

    proj_ra, proj_dec = projcat['ra']*u.deg, projcat['dec']*u.deg

    proj_coords = SkyCoord(proj_ra, proj_dec, frame='icrs')
    hj_coords = tab['sky_coord']

    for colname in [
        'match_proj_ra','match_proj_dec','match_sep_arcsec'
    ]:
        tab[colname] = np.nan

    hatids = []
    for hj_coord, name in zip(hj_coords, tab['pl_name']):

        seps = hj_coord.separation(proj_coords)
        projsorted = proj_coords[np.argsort(seps)]
        sepssorted = proj_coords[np.argsort(seps)]

        # now get the HATIDs of the match
        if np.any(seps < minxmatchsep):
            sel = (seps < minxmatchsep)

            if len(sel[sel]) > 1:
                raise NotImplementedError('expected a single HJ HATID match')

            projcatmatch = projcat[np.argmin(seps)]

            thispl = (tab['pl_name']==name)
            hatids.append(projcatmatch[0].encode('ascii'))
            tab[thispl]['match_proj_ra'] = projcatmatch[1]
            tab[thispl]['match_proj_dec'] = projcatmatch[2]
            match_sep = np.min(seps).to(u.arcsec).value
            tab[thispl]['match_sep_arcsec'] = match_sep

            print('{}: got projcatalog match with separation {} arcsec'.
                  format(name, match_sep))

        else:
            hatids.append(np.nan)
            print('{}: did not get projcatalog match with separation < {}'.
                  format(name, minxmatchsep))
            continue

    tab['hatids'] = hatids

    # check if the HJs on chip with known HATIDs have lightcurves
    rle, ele, tle, rlcs, elcs, tlcs = [],[],[],[],[],[]
    for hatid in tab['hatids']:

        if hatid=='nan':
            rawlc = str(np.nan)
            epdlc = str(np.nan)
            tfalc = str(np.nan)
        else:
            rawlc = glob(os.path.join(lcdirectory, hatid+'.grcollectilc'))
            epdlc = glob(os.path.join(lcdirectory, hatid+'.epdlc'))
            tfalc = glob(os.path.join(lcdirectory, hatid+'.tfalc'))

            for lc in [rawlc, epdlc, tfalc]:
                if len(lc) > 1:
                    raise ValueError('expected at most 1 lightcurve match')

            rawlc = rawlc[0] if len(rawlc)==1 else str(np.nan)
            epdlc = epdlc[0] if len(epdlc)==1 else str(np.nan)
            tfalc = tfalc[0] if len(tfalc)==1 else str(np.nan)

        rawlcexists = os.path.exists(rawlc)
        epdlcexists = os.path.exists(epdlc)
        tfalcexists = os.path.exists(tfalc)

        rle.append(rawlcexists)
        ele.append(epdlcexists)
        tle.append(tfalcexists)
        rlcs.append(rawlc)
        elcs.append(epdlc)
        tlcs.append(tfalc)

    tab['rawlcexists'] = nparr(rle)
    tab['epdlcexists'] = nparr(ele)
    tab['tfalcexists'] = nparr(tle)
    tab['rawlc'] = nparr(rlcs)
    tab['epdlc'] = nparr(elcs)
    tab['tfalc'] = nparr(tlcs)

    print(tab['pl_name', 'hatids', 'rawlcexists', 'epdlcexists',
              'tfalcexists'])
    print(tab['rawlc'])

    if not np.any(tab['tfalcexists']):

        print('WRN! did not find any known hot jupiters that had TFA '
              'lightcurves')

        return 0

    elif not np.any(tab['tfalcexists']) and np.any(tab['epdlcexists']):

        print('WRN! but found EPD lcs for some HJs. Perhaps TFA is failing.')

        return 0

    elif not np.any(tab['tfalcexists']) and not np.any(tab['epdlcexists']):
        print('WRN! and did not find any EPD lightcurves either.')

        return 0

    # otherwise, continue to by measuring the HJ's SNR
    for tfalc in tab['tfalc'][tab['tfalcexists']]:
        _measure_hj_snr(tfalc)

def _measure_hj_snr(tfalc):

    #FIXME: implement
    raise NotImplementedError
    pass
    # read the lc

    # BLS

    # fit trapezoidal model

    # use standard model fitting to get the SNR

    # probably make a plot



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




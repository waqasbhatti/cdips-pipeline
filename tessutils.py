# -*- coding: utf-8 -*-
"""
functions for reducing TESS data.  contents are as follows, where "sub-tasks"
in a conceptual sense are listed below the main task:
------------------------------------------

parallel_mask_saturated_stars: mask saturated stars given saturation level

    mask_saturated_stars_worker

parallel_mask_dquality_flag_frames: mask entire frames based on DQUALITY flag

    mask_dquality_flag_frame

parallel_trim_get_single_extension: get single extension image, trim to remove
virtual columns, append "PROJID" header keywork

    from_CAL_to_fitsh_compatible
    (deprecated) from_ete6_to_fitsh_compatible

are_known_planets_in_field: precursor to measuring known planet properties is
knowing if there are any

measure_known_planet_SNR: if there are known planets in this frame, estimate
their SNR through a max-likelihood trapezoidal transit model

    _measure_planet_snr

read_object_reformed_catalog: get objectid, ra, dec from reformed catalog

read_object_catalog: read the Gaia .catalog file into a pandas DataFrame

read_tess_txt_lightcurve: Read grcollect dumped lightcurve file into a
dictionary (which can then be put into a FITS table, or some other data
structure better than text).
"""
from __future__ import division, print_function

##########################################

import aperturephot as ap
import imageutils as iu
import shared_variables as sv

import os, pickle, re
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from glob import glob
from parse import search

from astropy.io import fits
from astroquery.mast import Catalogs
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

from numpy import array as nparr, all as npall

from datetime import datetime
import multiprocessing as mp

from astropy.coordinates import SkyCoord
from astropy import units as u, constants as const
from astropy.io import ascii

from lcstatistics import read_tfa_lc, plot_raw_epd_tfa

from astrobase.periodbase import kbls
from astrobase.varbase import lcfit
from astrobase.varbase.transits import get_transit_times

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
    hdr.comments['PROJID'] = 'trex identifier for parameters in this run'

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


def are_known_planets_in_field(ra_center, dec_center, outname, use_NEA=False,
                               use_alerts=True):
    """
    Given a field center, find which known planets are on chip.
    Dependencies: tessmaps, astroquery

    Args:
        ra_center, dec_center (floats): coordinate of field center in degrees

        outname (str): path to which csv file with details of HJs on chip will
        be written

        useNEA (bool): whether to use NASA exoplanet archive to grab the known
        planet properties (in particular, known HJs), or not.

        use_alerts (bool): whether to use the "alert" file, from
        https://archive.stsci.edu/prepds/tess-data-alerts/index.html, to get
        candidate planet host (TOI) properties.

    Returns:
        True if any HJs on chip, else False.
    """

    from tessmaps import get_time_on_silicon as gts

    cam_dirn = SkyCoord(ra_center*u.deg, dec_center*u.deg, frame='icrs')
    cam_tuple = (cam_dirn.barycentrictrueecliptic.lat.value,
                 cam_dirn.barycentrictrueecliptic.lon.value)

    if use_NEA:
        eatab = NasaExoplanetArchive.get_confirmed_planets_table()

        is_jup = (eatab['pl_radj'] > 0.4*u.Rjup)
        is_hot = (eatab['pl_orbper'] < 10*u.day)
        is_hj = is_jup & is_hot
        hj_coords = eatab[is_hj]['sky_coord']

        pl_onchip = gts.given_one_camera_get_stars_on_silicon(
            hj_coords, cam_tuple, withgaps=False)

    elif use_alerts:
        df = pd.read_csv(
            '/home/lbouma/proj/pipe-trex/data/'
            'hlsp_tess-data-alerts_tess_phot_alert-summary-s01+s02+s03+s04_tess_v9_spoc.csv'
        )
        kp_coords = SkyCoord(nparr(df['RA'])*u.deg,
                             nparr(df['Dec'])*u.deg,
                             frame='icrs')

        pl_onchip = gts.given_one_camera_get_stars_on_silicon(
            kp_coords, cam_tuple, withgaps=False)

    else:
        raise NotImplementedError('use_alerts or use_NEA must be true.')


    if np.any(pl_onchip) and use_NEA:
        desiredcols = ['pl_hostname', 'pl_letter', 'pl_name', 'pl_orbper',
                       'pl_radj', 'st_optmag', 'gaia_gmag','sky_coord']

        outtab = eatab[is_hj][desiredcols]
        outtab = outtab[pl_onchip.astype(bool)]

        ascii.write(outtab, output=outname, format='ecsv',
                    overwrite=True)
        print('wrote {}'.format(outname))

        return True

    elif np.any(pl_onchip) and use_alerts:

        outdf = df[pl_onchip.astype(bool)]
        outdf.to_csv(outname, index=False)
        print('wrote {}'.format(outname))

        return True

    else:
        return False


def are_clusters_in_field(ra_center, dec_center, outname, dist_cut=2000,
                          sectornum=1):
    """
    Given a field center, find which clusters from Kharchenko+ 2013 are on chip.
    Limit to clusters within some distance, say 2kpc.
    Dependencies: tessmaps

    Args:
        ra_center, dec_center (floats): coordinate of field center in degrees

        outname (str): path to which csv file with details of HJs on chip will
        be written

    Kwargs:

        dist_cut (float): number of parsecs beyond which we are not interested

        sectornum (int): 1-based "science sector" count.

    Returns:
        True if any HJs on chip, else False.
    """

    from tessmaps import get_time_on_silicon as gts

    # get the coordinates and properties of the clusters
    k13path = ('/home/lbouma/proj/pipe-trex/data/'
               'Kharchenko_2013_MWSC_tess_sectors.csv')
    k13 = pd.read_csv(k13path)

    sel_clusters = (k13['sector_{:d}'.format(sectornum-1)]==1)
    sel_clusters &= (k13['d'] < dist_cut)

    df = k13[sel_clusters]
    desiredcols = ['MWSC', 'Name', 'Type', 'n_Type', 'RAJ2000', 'DEJ2000',
                   'r0', 'r1', 'r2', 'N1sr0', 'N1sr1', 'N1sr2', 'd', 'logt',
                   'ra', 'dec']
    df = df[desiredcols]

    # get coordinates
    cluster_coords = SkyCoord(nparr(df['ra'])*u.deg, nparr(df['dec'])*u.deg,
                              frame='icrs')

    ##########################################
    cam_dirn = SkyCoord(ra_center*u.deg, dec_center*u.deg, frame='icrs')
    cam_tuple = (cam_dirn.barycentrictrueecliptic.lat.value,
                 cam_dirn.barycentrictrueecliptic.lon.value)

    ##########################################
    cluster_onchip = gts.given_one_camera_get_stars_on_silicon(
        cluster_coords, cam_tuple, withgaps=False)

    if np.any(cluster_onchip):

        outdf = df
        outdf = outdf[cluster_onchip.astype(bool)]

        outdf.to_csv(outname, index=False)
        print('wrote {}'.format(outname))

        return True

    else:
        return False


def cluster_cutout_jpg_worker(task):

    (calfitsimage, subfitsimage, wcsfile, fitsdir, cname, ctype, nstar, dist,
     logt, ra, dec, angrad, clusterdistancecutoff ) = task

    outcalpath = os.path.join(
        fitsdir, 'CUT-{}_{}_CAL.jpg'.
        format(cname,
               os.path.basename(calfitsimage.replace('.fits','')))
    )
    outsubpath = os.path.join(
        fitsdir, 'CUT-{}_{}_SUB.jpg'.
        format(cname,
               os.path.basename(subfitsimage.replace('.fits','')))
    )

    if dist > clusterdistancecutoff:
        print('skipping {}, d>{}pc'.
              format(cname, clusterdistancecutoff))
        return 1

    if os.path.exists(outcalpath) and os.path.exists(outsubpath):
        return 1

    thisctype, thisnstar, thislogt = str(ctype), int(nstar), float(dist)

    ra, dec = float(ra), float(dec)
    clusterangradius = (float(angrad)*u.deg).value
    boxwidth, boxheight = 5*clusterangradius, 5*clusterangradius

    radeccenter = [ra, dec, boxwidth, boxheight]

    annotatestr = '{:s}-{:s}, {:s}*, logt={:s}, d={:.1f}'.format(
        cname, ctype, repr(int(nstar)),
        '{:.1f}'.format(logt), dist/1000
    )

    try:
        iu.frame_radecbox_to_jpeg(calfitsimage, wcsfrom=wcsfile,
                                  radeccenter=radeccenter,
                                  out_fname=outcalpath,
                                  annotatejd=False,
                                  annotate=annotatestr,
                                  forcesquare=True,
                                  overplotscalebar=True,
                                  rescaleimage=True,
                                  scale_func=iu.clipped_logscale_img,
                                  scale_func_params={
                                      'cap':255.0, 'lomult':2.0,
                                      'himult':2.0, 'coeff':1000.0},
                                  verbose=False
                                 )
        iu.frame_radecbox_to_jpeg(subfitsimage, wcsfrom=wcsfile,
                                  radeccenter=radeccenter,
                                  out_fname=outsubpath.replace('.jpg','_grayscale.jpg'),
                                  annotatejd=False,
                                  annotate=annotatestr,
                                  forcesquare=True,
                                  overplotscalebar=True,
                                  rescaleimage=True,
                                  verbose=False)
        iu.frame_radecbox_to_jpeg(subfitsimage, wcsfrom=wcsfile,
                                  radeccenter=radeccenter,
                                  out_fname=outsubpath.replace('.jpg','_bwr.jpg'),
                                  annotatejd=False,
                                  annotate=annotatestr,
                                  forcesquare=True,
                                  overplotscalebar=True,
                                  rescaleimage=True,
                                  colorscheme='bwr',
                                  verbose=False)
    except Exception as e:
        print(e)
        print('{}, {}'.format(cname, repr(radeccenter)))
        print('failed to cut for {}, {}'.format(outcalpath, outsubpath))
        return -1

    return 1


def parallel_cluster_cutout_jpgs(calfitsimages, subfitsimages, wcsfiles,
                                 fitsdir, cname, ctype, nstar, dist, logt,
                                 ra, dec, angrad, clusterdistancecutoff,
                                 nworkers=16, maxworkertasks=1000):

    print('%sZ: %s files to make cluster cutouts for %s' %
          (datetime.utcnow().isoformat(), len(wcsfiles), cname))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, y, z, fitsdir, cname, ctype, nstar, dist, logt, ra, dec,
              angrad, clusterdistancecutoff)
             for x,y,z in zip(calfitsimages, subfitsimages, wcsfiles)
            ]

    # fire up the pool of workers
    results = pool.map(cluster_cutout_jpg_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {result for result in results}




def make_cluster_cutout_jpgs(sectornum, fitsdir, racenter, deccenter, field,
                             camera, ccd, statsdir,
                             clusterdistancecutoff=2000, nworkers=16):
    """
    (kw)args:
        clusterdistancecutoff (float): distance beyond which you are not making
        cutouts for clusters.
    """

    wcsfiles = glob(os.path.join(fitsdir,'*.wcs'))
    calfitsimages = [w.replace('.wcs','.fits') for w in wcsfiles]
    subfitsimages = glob(os.path.join(fitsdir,'rsub-*-tess*-xtrns.fits'))

    if not len(wcsfiles)==len(subfitsimages)==len(calfitsimages):
        raise AssertionError

    outcsv = os.path.join(
        statsdir, '{:s}_cam{:d}_ccd{:d}_kharchenko13_clusters.csv'.
        format(field, camera, ccd)
    )

    if are_clusters_in_field(racenter, deccenter, outcsv, sectornum=sectornum):
        pass
    else:
        return -1

    # columns: Index([u'MWSC', u'Name', u'Type', u'n_Type', u'RAJ2000', u'DEJ2000',
    # u'r0', u'r1', u'r2', u'N1sr0', u'N1sr1', u'N1sr2', u'd', u'logt', u'ra',
    # u'dec'], dtype='object')
    df = pd.read_csv(outcsv)

    names = np.array(df['Name'])
    ctype = np.array(df['Type'].fillna(''))
    nstar = np.array(df['N1sr2'])
    dist = np.array(df['d'])
    logt = np.array(df['logt'])
    ras = np.array(df['ra'])
    decs = np.array(df['dec'])
    angrads = np.array(df['r2'])

    for name, ct, ns, d, age, ra, dec, r2 in zip(names, ctype, nstar, dist,
                                                 logt, ras, decs, angrads):

        outcalmatches = glob(os.path.join(
            fitsdir, 'CUT-{}*_CAL.jpg'.format(name) ))
        if len(outcalmatches) > 10:
            print('found cuts for {}, continue'.format(name))
            continue

        parallel_cluster_cutout_jpgs(calfitsimages, subfitsimages, wcsfiles,
                                     fitsdir, name, ct, ns, d, age, ra, dec,
                                     r2, clusterdistancecutoff,
                                     nworkers=nworkers)



def measure_known_planet_SNR(kponchippath, projcatalogpath, lcdirectory,
                             statsdir, sectornum, minxmatchsep=3, nworkers=20,
                             use_NEA=False, use_alerts=True):
    """
    Args:

        kponchippath (str): path to csv that tells you which HJs are expected
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

    if use_NEA:
        # read information table from `tessutils.are_known_planets_in_field`
        tab = ascii.read(kponchippath)

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
    elif use_alerts:
        tab = pd.read_csv(kponchippath)
        tab['pl_name'] = tab['toi_id']
        tab['pl_hostname'] = tab['tic_id']
    else:
        raise NotImplementedError

    # get the starids that correspond to the HJs on-chip
    projcat = read_object_reformed_catalog(projcatalogpath, isgaiaid=True)

    proj_ra, proj_dec = projcat['ra']*u.deg, projcat['dec']*u.deg

    proj_coords = SkyCoord(proj_ra, proj_dec, frame='icrs')
    if use_NEA:
        kp_coords = tab['sky_coord']
    elif use_alerts:
        kp_coords = SkyCoord(nparr(tab['RA'])*u.deg,
                             nparr(tab['Dec'])*u.deg,
                             frame='icrs')

    for colname in [
        'match_proj_ra','match_proj_dec','match_sep_arcsec'
    ]:
        tab[colname] = np.nan

    starids = []
    for kp_coord, name in zip(kp_coords, tab['pl_name']):

        seps = kp_coord.separation(proj_coords)
        projsorted = proj_coords[np.argsort(seps)]
        sepssorted = proj_coords[np.argsort(seps)]

        # now get the starids of the match
        if np.any(seps < minxmatchsep):
            sel = (seps < minxmatchsep)

            if len(sel[sel]) > 1:
                print('got {} matches'.format(len(sel[sel])))
                print('separations are {}'.format(
                    repr(seps[sel].to(u.arcsec)))
                )
                print('taking minimum separation...')

            projcatmatch = projcat[np.argmin(seps)]

            thispl = (tab['pl_name']==name)
            starids.append(projcatmatch[0].encode('ascii'))
            tab[thispl]['match_proj_ra'] = projcatmatch[1]
            tab[thispl]['match_proj_dec'] = projcatmatch[2]
            match_sep = np.min(seps).to(u.arcsec).value
            tab[thispl]['match_sep_arcsec'] = match_sep

            print('{}: got projcatalog match with separation {} arcsec'.
                  format(name, match_sep))

        else:
            starids.append(np.nan)
            print('{}: did not get projcatalog match with separation < {}'.
                  format(name, minxmatchsep))
            continue

    tab['starids'] = starids

    # check if the HJs on chip with known starids have lightcurves
    rle, ele, tle, rlcs, elcs, tlcs = [],[],[],[],[],[]
    for starid in tab['starids']:

        if starid=='nan' or pd.isnull(starid):
            rawlc = str(np.nan)
            epdlc = str(np.nan)
            tfalc = str(np.nan)
        else:
            if type(starid) == bytes:
                starid = starid.decode('utf-8')
            rawlc = glob(os.path.join(lcdirectory, starid+'_llc.fits'))
            epdlc = glob(os.path.join(lcdirectory, starid+'_llc.fits'))
            tfalc = glob(os.path.join(lcdirectory, starid+'_llc.fits'))

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

    print(tab[['pl_name', 'starids', 'rawlcexists', 'epdlcexists',
               'tfalcexists']])
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
    for tfalc, plname in zip(
        tab['tfalc'][tab['tfalcexists']], tab['pl_name'][tab['tfalcexists']]
    ):
        _measure_planet_snr(plname, tfalc, statsdir, sectornum,
                            nworkers=nworkers, use_NEA=use_NEA,
                            use_alerts=use_alerts)


def _measure_planet_snr(plname, tfalc, statsdir, sectornum,
                        timename='TMID_BJD', nworkers=20,
                        n_transit_durations=4, use_NEA=False, use_alerts=True):

    plname = str(plname).replace(' ','')
    plname = str(plname).replace('.','-')

    # we only care about TFA lightcurves for planet SNR measurements.
    dtrstages = ['TFA']
    aps = [str(r) for r in range(1,4)]
    magnames = [dtrstage+ap for dtrstage in dtrstages for ap in aps]
    errnames = ['IRE'+ap for dtrstage in dtrstages for ap in aps]

    # read the lc
    hdulist = fits.open(tfalc)
    lc = hdulist[1].data

    time = lc[timename]

    plotpath = os.path.join(statsdir, str(plname)+'_AP1_lightcurve.png')
    plot_raw_epd_tfa(time, lc['IRM1'], lc['EP1'], lc['TFA1'], 1,
                     savpath=plotpath, xlabel='BTJD = BJD - 2457000')

    outdir = os.path.join(statsdir, plname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for magname, errname in zip(magnames, errnames):

        blsfit_savfile = os.path.join(
            outdir, str(plname)+'_'+str(magname)+'_BLS_fit.png')
        trapfit_savfile = os.path.join(
            outdir, str(plname)+'_'+str(magname)+'_trapezoidal_fit.png')
        snrfit_savfile = os.path.join(
            outdir, str(plname)+'_'+str(magname)+'_snr_of_fit.csv')

        mag = lc[magname]
        err_mag = lc[errname]

        mag_0 = 10  # these numbers do not matter; we care about relative flux
        flux_0 = 1000
        flux = flux_0 * 10**(-0.4 * (mag - mag_0))
        # sigma_flux = dg/d(mag) * sigma_mag, for g=f0 * 10**(-0.4*(mag-mag0)).
        err_flux = np.abs(
            -0.4 * np.log(10) * flux_0 * 10**(-0.4*(mag-mag_0)) * err_mag
        )

        fluxmedian = np.nanmedian(flux)
        flux /= fluxmedian
        err_flux /= fluxmedian

        # fit BLS model; plot resulting phased LC.
        endp = 1.05*(np.nanmax(time) - np.nanmin(time))/2
        blsdict = kbls.bls_parallel_pfind(time, flux, err_flux,
                                          magsarefluxes=True, startp=0.1,
                                          endp=endp, maxtransitduration=0.15,
                                          nworkers=nworkers, sigclip=None)
        fitd = kbls.bls_stats_singleperiod(time, flux, err_flux,
                                           blsdict['bestperiod'],
                                           maxtransitduration=0.15,
                                           magsarefluxes=True, sigclip=None,
                                           perioddeltapercent=5)
        bls_period = fitd['period']
        lcfit._make_fit_plot(fitd['phases'], fitd['phasedmags'], None,
                             fitd['blsmodel'], fitd['period'], fitd['epoch'],
                             fitd['epoch'], blsfit_savfile, magsarefluxes=True)
        ingduration_guess = fitd['transitduration']*0.2
        transitparams = [fitd['period'], fitd['epoch'], fitd['transitdepth'],
                         fitd['transitduration'], ingduration_guess ]

        # fit initial trapezoidal model; plot resulting phased LC.
        trapfit = lcfit.traptransit_fit_magseries(time, flux, err_flux,
                                                  transitparams,
                                                  magsarefluxes=True,
                                                  sigclip=None,
                                                  plotfit=trapfit_savfile)

        # isolate each transit to within +/- n_transit_durations
        tmids, t_starts, t_ends = (
            get_transit_times(fitd, time, n_transit_durations, trapd=trapfit)
        )

        fitmags = trapfit['fitinfo']['fitmags']
        fitflux = fitmags

        # estimate snr direct from the trapezoidal fit
        t_full_durn_phase = trapfit['fitinfo']['finalparams'][3]
        t_ing_durn_phase = trapfit['fitinfo']['finalparams'][4]
        period = trapfit['fitinfo']['finalparams'][0]

        t_full_durn = t_full_durn_phase * period
        t_ing_dur = t_ing_durn_phase * period

        per_point_cadence=30.*u.min

        npoints_in_transit = (
            float(((t_full_durn*u.day)/per_point_cadence).cgs.value)
        )

        transitdepth = np.abs(trapfit['fitinfo']['finalparams'][2])

        # get the indices in transit, for RMS measurement
        (in_fluxs, time_list,
         intra_inds_list, windowinds ) = [], [], [], []
        for t_start,t_end in zip(t_starts, t_ends):
            this_window_inds = (time > t_start) & (time < t_end)
            tmid = t_start + (t_end-t_start)/2
            # flag out slightly more than expected "in transit" points
            prefactor = 1.05
            transit_start = tmid - prefactor*t_full_durn/2
            transit_end = tmid + prefactor*t_full_durn/2

            this_window_intra = (
                (time[this_window_inds] > transit_start) &
                (time[this_window_inds] < transit_end)
            )
            this_window_oot = ~this_window_intra

            windowinds.append(this_window_inds)
            time_list.append( time[this_window_inds] )
            in_fluxs.append( flux[this_window_inds] )
            intra_inds_list.append( (time[this_window_inds]>transit_start) &
                                    (time[this_window_inds]<transit_end) )

        intra = np.hstack(np.array(intra_inds_list))
        sel_time = np.array(time_list)
        sel_flux = np.array(in_fluxs)
        windowinds = np.bitwise_or.reduce( np.array(windowinds), axis=0 )

        try:
            subtractedrms = np.std(flux[windowinds][~intra] -
                                   fitflux[windowinds][~intra] )

            ntransits = len(t_starts)
            trapz_snr = (
                np.sqrt(npoints_in_transit*ntransits) *
                transitdepth/subtractedrms
            )

            outdf = pd.DataFrame({'plname':plname, 'trapz_snr':trapz_snr,
                                  'tfalc':tfalc}, index=[0])
            outdf.to_csv(snrfit_savfile, index=False)
            print('made {}'.format(snrfit_savfile))

        except IndexError:

            outdf = pd.DataFrame({'plname':plname, 'trapz_snr':np.nan,
                                  'tfalc':tfalc}, index=[0])
            outdf.to_csv(snrfit_savfile, index=False)
            print('WRN! SNR calculation failed, but made {} anyway with nan'.
                  format(snrfit_savfile))
            print('made {}'.format(snrfit_savfile))

    return 1


def explore_ccd_temperature_timeseries():
    """
    This is a one-off exploration to understand the engineering files below
    better. Some take-aways:

        * 'tess2018331094053_sector01-eng.fits' and
        'tess2018323111417_sector01-eng.fits' have identical temperature time
        series. It's unclear why there are two of these data files, but it
        doesn't matter for our purposes. This was confirmed by email from Scott
        Fleming (2019/01/04).

        * We can extract the time-series of on-chip CCD temperature using the
        S_CAM{:d}_ALCU_sensor_CCD{:d} pattern.

        * The cadence of the above time-series is different from FFIs, and will
        need to be binned to match the FFI cadence and assign an average
        temperature to each FFI.

        * There is one other type of CCD-specific temperature available:
            * Board temperature S_CAM[1-4]_CCD[1-4]_board_temp

        * There are also a variety of camera-specific temperatures (see TESS
        Instrument Handbook, Appendix A.5 Engineering Data).
    """

    engdatadir = '/nfs/phtess1/ar1/TESS/FFI/ENGINEERING/'
    s01_first = 'tess2018331094053_sector01-eng.fits'
    s01_second = 'tess2018323111417_sector01-eng.fits'
    s02 = 'tess2018330083923_sector02-eng.fits'

    d = {}
    for dfile in [s01_first, s01_second, s02]:

        d[dfile] = {}

        hdulist = fits.open(os.path.join(engdatadir, dfile))

        ccds, cams = list(range(1,5)), list(range(1,5))

        # NOTE: unclear if "COOKED" is actually the temperature... plausible
        # but unclear.
        for cam in cams:
            for ccd in ccds:
                TEMPERATURE_HDR_NAME = 'S_CAM{:d}_ALCU_sensor_CCD{:d}'.format(cam, ccd)
                this_time = hdulist[TEMPERATURE_HDR_NAME].data['TIME']
                this_temperature = hdulist[TEMPERATURE_HDR_NAME].data['COOKED']

                d[dfile][TEMPERATURE_HDR_NAME] = {}
                d[dfile][TEMPERATURE_HDR_NAME]['time'] = this_time
                d[dfile][TEMPERATURE_HDR_NAME]['temperature'] = this_temperature


    temperature_hdr_names = ['S_CAM{:d}_ALCU_sensor_CCD{:d}'.format(cam, ccd)
                             for cam in cams for ccd in ccds]

    times = [d[k][thn]['time'] for k in d.keys() for thn in
             temperature_hdr_names]

    temperatures = [d[k][thn]['temperature'] for k in d.keys() for thn in
                    temperature_hdr_names]

    f,ax = plt.subplots()
    for k in d.keys():
        for temperature_hdr_name in temperature_hdr_names:
            if 'S_CAM1_ALCU_sensor_CCD1' not in temperature_hdr_name:
                continue
            if 'tess2018331094053_sector01-eng.fits' in k:
                ax.plot(d[k][temperature_hdr_name]['time'],
                        d[k][temperature_hdr_name]['temperature']+3,
                        label='{:s}_{:s}'.format(k, temperature_hdr_name))
            else:
                ax.plot(d[k][temperature_hdr_name]['time'],
                        d[k][temperature_hdr_name]['temperature'],
                        label='{:s}_{:s}'.format(k, temperature_hdr_name))

    ax.set_xlabel('BTJD', fontsize='small')
    ax.set_ylabel('On-chip CCD temperature sensor ("COOKED") [deg C]',
                  fontsize='small')

    ax.legend(loc='best',fontsize='xx-small')

    f.savefig('temperature_sanity_check.png', dpi=300)

    # ok, now check if they're actually THE SAME, for every single ccd. If so,
    # then having these two engineering files for sector 1 is pointless.
    for temperature_hdr_name in temperature_hdr_names:
        temperature_1 = d[s01_first][temperature_hdr_name]['temperature']
        temperature_2 = d[s01_second][temperature_hdr_name]['temperature']

        time_1 = d[s01_first][temperature_hdr_name]['time']
        time_2 = d[s01_second][temperature_hdr_name]['time']

        np.testing.assert_array_almost_equal(temperature_1, temperature_2,
                                             decimal=6)
        np.testing.assert_array_almost_equal(time_1, time_2,
                                             decimal=6)

    print('Test passed: {:s} and {:s} have identical temperature timeseries'.
          format(s01_first, s01_second))



def append_ccd_temperature_to_hdr_worker(task):

    fitsname, d = task

    framekey = os.path.splitext(os.path.basename(fitsname))[0]

    # get the temperature time series appropriate for this CAM/CCD pair using
    # the file name.  match: tess2018206192942-s0001-4-3-0120_cal_img.fits
    sr = search('{}/tess2{}-{}-{}-{}-{}_cal_img.fits', fitsname)
    cam = sr[3]
    ccd = sr[4]

    thiskey = 'S_CAM{:d}_ALCU_sensor_CCD{:d}'.format(int(cam), int(ccd))

    time = d[thiskey]['time'] # engineering data time in "TJD = BJD-2457000".
    temperature = d[thiskey]['temperature'] # ditto, for CCD temperature

    # take the mean of the temperature values within this time window. append
    # it to the header.
    data, hdr = iu.read_fits(fitsname, ext=0)
    tstart = hdr['TSTART'] # start as BTJD = BJD-2457000.
    tstop = hdr['TSTOP'] # stop in BTJD, ditto.

    sel_times = (time>=tstart) & (time<tstop)

    mean_temp = np.mean( temperature[sel_times] )
    n_temps = len(temperature[sel_times])

    if n_temps >= 1:
        hdr['CCDTEMP'] = mean_temp
        hdr['NTEMPS'] = n_temps
    else:
        hdr['CCDTEMP'] = np.nan
        hdr['NTEMPS'] = 0

    hdr.comments['CCDTEMP'] = (
        'mean temperature (S_CAM{:d}_ALCU_sensor_CCD{:d})'.
        format(int(cam), int(ccd))
    )
    hdr.comments['NTEMPS'] = (
        'number of temperatures avgd for CCDTEMP'
    )

    fits.writeto(fitsname, data, header=hdr, output_verify='silentfix+ignore',
                 overwrite=True)

    if n_temps>=1:
        print('Wrote CCDTEMP to {:s}'.format(fitsname))
        return framekey, mean_temp, n_temps
    else:
        print('WRN! Wrote NAN CCDTEMP to {:s}'.format(fitsname))
        return framekey, np.nan, 0


def parallel_append_ccd_temperature_to_hdr(fitslist, temperaturepklpath,
                                           nworkers=16, maxworkertasks=1000):
    # pass a list of trim, cut cal images:
    # tess2018206195942-s0001-4-3-0120_cal_img.fits

    d = pickle.load(open(temperaturepklpath, 'rb'))

    print('%sZ: %s files to append CCD temperature to header' %
          (datetime.utcnow().isoformat(), len(fitslist)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, d) for x in fitslist]

    # fire up the pool of workers
    results = pool.map(append_ccd_temperature_to_hdr_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    print('%sZ: finished appending CCD temperatures to headers' %
          (datetime.utcnow().isoformat()))

    # write CCDTEMP, NTEMPS, and FRAMEKEY csv file for later quick-read during
    # FITS LC creation.
    outdf = pd.DataFrame(results, columns=['framekey','ccdtemp','ntemps'])

    framekey = outdf['framekey'].iloc[0]
    res = search('tess{}-{}_cal_img', framekey)
    keystring = res[-1]
    outname = '{}_key_temperature_count.csv'.format(keystring)

    # output csv path example:
    # /nfs/phtess1/ar1/TESS/FFI/ENGINEERING/s0002-3-3-0121_key_temperature_count.csv
    outdir = os.path.dirname(temperaturepklpath)
    outpath = os.path.join(outdir, outname)

    outdf.to_csv(outpath, index=False)
    print('made {}'.format(outpath))

    return {result for result in results}



def make_ccd_temperature_timeseries_pickle(sectornum):
    """
    convert the engineering data on MAST to a pickled time-series of on-chip
    CCD temperature.

    outputs:
        '/nfs/phtess1/ar1/TESS/FFI/ENGINEERING/sector0001_ccd_temperature_timeseries.pickle'

        which contains a dictionary, with keys 'S_CAM{:d}_ALCU_sensor_CCD{:d}',
        each of which contains a dictionary of time and temperature.
    """

    engdatadir = '/nfs/phtess1/ar1/TESS/FFI/ENGINEERING/'
    if sectornum==1:
        # identical to the other one on mast. LGB tested this in
        # explore_ccd_temperature_timeseries
        fname = os.path.join(engdatadir,'tess2018323111417_sector01-eng.fits')
    elif sectornum>1:
        fnames = glob(os.path.join(
            engdatadir,'tess2*_sector{:s}-eng.fits'.
            format(str(sectornum).zfill(2))
        ))
        if len(fnames) > 1:
            raise AssertionError('expected one engineering file per sector')
        fname = fnames[0]
    else:
        raise NotImplementedError

    hdulist = fits.open(fname)
    ccds, cams = list(range(1,5)), list(range(1,5))

    # NOTE: unclear if "COOKED" is actually the temperature... plausible
    # but unclear.
    d = {}
    temperature_hdr_names = ['S_CAM{:d}_ALCU_sensor_CCD{:d}'.format(cam, ccd)
                             for cam in cams for ccd in ccds]

    for temperature_hdr_name in temperature_hdr_names:
        this_time = hdulist[temperature_hdr_name].data['TIME']
        this_temperature = hdulist[temperature_hdr_name].data['COOKED']
        d[temperature_hdr_name] = {}
        d[temperature_hdr_name]['time'] = this_time
        d[temperature_hdr_name]['temperature'] = this_temperature

    # note: the read out times on each camera are not identical. therefore we
    # can't make this extremely simple as just temperatures in the 16 readouts,
    # vs time. it needs to be: for each readout, what is the temperature
    # time-series?

    pklpath = os.path.join(
        engdatadir,
        'sector{:s}_ccd_temperature_timeseries.pickle'.
        format(str(sectornum).zfill(4))
    )
    with open(pklpath, 'wb') as f:
        pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
    print('made {}'.format(pklpath))


def read_tess_txt_lightcurve(
    lcfile,
    catalogdf,
    infokeys=['Gaia-ID', 'RA[deg]', 'Dec[deg]', 'Parallax[mas]',
              'PM_RA[mas/yr]', 'PM_Dec[mas/year]', 'Ref_Epoch[yr]',
              'AstExcNoise[mas]', 'AstExcNoiseSig', 'AstPriFlag',
              'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag',
              'radial_velocity', 'phot_variable_flag', 'teff_val', 'a_g_val',
              'e_bp_min_rp_val', 'radius_val', 'lum_val']
    ):
    """
    Read grcollect dumped lightcurve file into a dictionary (which can then be
    put into a FITS table, or some other data structure better than text
    files).

    Args:
        lcfile (str): path to *.grcollectilc, matching format below.

        catalogdf (pd.DataFrame): output from read_object_catalog(catfile)

    catfile

    Format details:

    00 tmid_utc exp mid-time in JD_UTC (from DATE-OBS,DATE-END)
    01 rstfc    Unique frame key ({STID}-{FRAMENUMBER}_{CCDNUM})
    02 starid   GAIA ID of the object
    03 xcc      original X coordinate on CCD on photref frame
    04 ycc      original y coordinate on CCD on photref frame
    05 xic      shifted X coordinate on CCD on subtracted frame
    06 yic      shifted Y coordinate on CCD on subtracted frame
    07 fsv      Measured S value
    08 fdv      Measured D value
    09 fkv      Measured K value
    10 bgv      Background value
    11 bge      Background measurement error

    12 ifl1     Flux in aperture 1 (ADU)
    13 ife1     Flux error in aperture 1 (ADU)
    14 irm1     Instrumental magnitude in aperture 1
    15 ire1     Instrumental magnitude error for aperture 1
    16 irq1     Instrumental magnitude quality flag for aperture 1 (0/G OK, X bad)

    17 ifl2     Flux in aperture 2 (ADU)
    18 ife2     Flux error in aperture 2 (ADU)
    19 irm2     Instrumental magnitude in aperture 2
    20 ire2     Instrumental magnitude error for aperture 2
    21 irq2     Instrumental magnitude quality flag for aperture 2 (0/G OK, X bad)

    22 ifl3     Flux in aperture 3 (ADU)
    23 ife3     Flux error in aperture 3 (ADU)
    24 irm3     Instrumental magnitude in aperture 3
    25 ire3     Instrumental magnitude error for aperture 3
    26 irq3     Instrumental magnitude quality flag for aperture 3 (0/G OK, X bad)
    """

    # read the LC into a numpy recarray
    lccolnames = ['tmid_utc','rstfc','starid',
                  'xcc','ycc','xic','yic','fsv','fdv','fkv',
                  'bgv','bge',
                  'ifl1', 'ife1', 'irm1', 'ire1', 'irq1',
                  'ifl2', 'ife2', 'irm2', 'ire2', 'irq2',
                  'ifl3', 'ife3', 'irm3', 'ire3', 'irq3' ]
    dtypelist = ('f8,U40,U19,'+
                 'f8,f8,f8,f8,f8,f8,f8,'+
                 'f8,f8,'+
                 'f8,f8,f8,f8,U1,'+
                 'f8,f8,f8,f8,U1,'+
                 'f8,f8,f8,f8,U1' )

    recarr = np.genfromtxt(lcfile,
                           usecols=tuple(range(len(lccolnames))),
                           names=lccolnames,
                           dtype=dtypelist)

    # generate the objectid by stripping the filename
    objectid = os.path.basename(lcfile).rstrip(os.path.splitext(lcfile)[-1])

    # look up this object in the catalog
    objind = catalogdf['Gaia-ID'] == int(objectid)

    if objind.size > 0:

        objectinfo = catalogdf[infokeys].ix[objind]
        objectinfo = objectinfo.to_dict('records')[0]

    else:

        objectinfo = {'objectid': objectid}

    # this is the lcdict we need
    lcdict = {'objectid': objectid,
              'objectinfo':objectinfo}
    for col in lccolnames:
        lcdict[col] = recarr[col]

    return lcdict


def read_object_reformed_catalog(reformedcatalogfile, isgaiaid=False):
    """
    you often need an objectid, ra, dec, and/or 2mass magnitudes.  get these
    from the "reformed" catalog with a small sub-read.
    """

    if not isgaiaid:
        columns='id,ra,dec,xi,eta,2massJ,2massK,2massqlt,2massI,2massr,2massi,2massz'
        columns = columns.split(',')

        catarr = np.genfromtxt(reformedcatalogfile,
                               comments='#',
                               usecols=list(range(len(columns))),
                               names=columns,
                               dtype='U15,f8,f8,f8,f8,f8,f8,U3,f8,f8,f8,f8',
                               delimiter=' '
                               )
    else:
        columns='id,ra,dec'
        columns = columns.split(',')

        catarr = np.genfromtxt(reformedcatalogfile,
                               comments='#',
                               usecols=list(range(len(columns))),
                               names=columns,
                               dtype='U19,f8,f8',
                               delimiter=' '
                               )
    return catarr


def read_object_catalog(catalogfile):
    """
    read the Gaia .catalog file into a pandas DataFrame
    """

    df = pd.read_csv(catalogfile, header=0, delim_whitespace=True)

    # columns are read in looking like "teff_val[35]". strip the [integer], and
    # the pre-pended "#" on the first column.
    cols = [re.sub('\[\d+\]', '', c) for c in df.columns]
    cols = [c.lstrip('#') for c in cols]

    df = df.rename(index=str,
                   columns={old:new for old, new in zip(df.columns, cols)})

    # columns are:
    # ['Gaia-ID', 'RA[deg]', 'Dec[deg]', 'RAError[mas]', 'DecError[mas]',
    # 'Parallax[mas]', 'Parallax_error[mas]', 'PM_RA[mas/yr]',
    # 'PM_Dec[mas/year]', 'PMRA_error[mas/yr]', 'PMDec_error[mas/yr]',
    # 'Ref_Epoch[yr]', 'AstExcNoise[mas]', 'AstExcNoiseSig', 'AstPriFlag',
    # 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error',
    # 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs',
    # 'phot_bp_mean_flux', 'phot_bp_mean_flux_error',
    # 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs',
    # 'phot_rp_mean_flux', 'phot_rp_mean_flux_error',
    # 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag',
    # 'phot_bp_rp_excess_factor', 'radial_velocity', 'radial_velocity_error',
    # 'phot_variable_flag', 'teff_val', 'teff_percentile_lower',
    # 'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower',
    # 'a_g_percentile_upper', 'e_bp_min_rp_val',
    # 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper',
    # 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper',
    # 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'xi[deg]',
    # 'eta[deg]']

    return df

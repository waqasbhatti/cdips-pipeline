#!/usr/bin/env python

"""
lcutils.py - Luke Bouma (luke@astro.princeton.edu) - Jan 2019

Tools for working with FITS lightcurves.
"""

import os
import os.path
import sys
import logging
from glob import glob
from datetime import datetime

import numpy as np, pandas as pd, matplotlib.pyplot as plt

from numpy import array as nparr, all as npall
np.seterr(all='ignore')

import multiprocessing as mp

import tessutils as tu, imageutils as iu

from astropy.io import fits

DEBUG = False


def _map_key_to_format(key):
    """
    given key, get FITS format code. Works for arbitrary number of apertures.  see
    http://docs.astropy.org/en/stable/io/fits/usage/table.html for
    details.
    """

    # example keys being mapped for TESS:
    # ['bge', 'bgv', 'fdv', 'fkv', 'fsv',
    #  'ife1', 'ife2', 'ife3', 'ifl1', 'ifl2',
    #  'ifl3', 'ire1', 'ire2', 'ire3', 'irm1', 'irm2', 'irm3', 'irq1', 'irq2',
    #  'irq3', 'tjd', 'rstfc', 'xic', 'yic']

    if key in ['bge','bgv','fdv','fkv','fsv','xic','yic']:
        return 'D'
    elif 'ife' in key: # flux err
        return 'D'
    elif 'ifl' in key: # flux
        return 'D'
    elif 'ire' in key: # raw mag err
        return 'D'
    elif 'irm' in key: # raw mag
        return 'D'
    elif 'irq' in key: # quality flag
        return '1A'
    elif 'tjd' in key: # timestamp
        return 'D'
    elif 'rstfc' in key: # frame id
        return '40A'

def _map_key_to_comment(k):
    kcd = {
        "tjd"   : "Spacecraft JD-2457000 in TDB",
        "rstfc" : "Unique frame key",
        "starid": "GAIA ID of the object",
        "xcc"   : "original X coordinate on CCD on photref frame",
        "ycc"   : "original y coordinate on CCD on photref frame",
        "xic"   : "shifted X coordinate on CCD on subtracted frame",
        "yic"   : "shifted Y coordinate on CCD on subtracted frame",
        "fsv"   : "Measured S value",
        "fdv"   : "Measured D value",
        "fkv"   : "Measured K value",
        "bgv"   : "Background value",
        "bge"   : "Background measurement error",
        "ifl1"  : "Flux in aperture 1 (ADU)",
        "ife1"  : "Flux error in aperture 1 (ADU)",
        "irm1"  : "Instrumental mag in aperture 1",
        "ire1"  : "Instrumental mag error for aperture 1",
        "irq1"  : "Instrumental quality flag ap 1, 0/G OK, X bad",
        "ifl2"  : "Flux in aperture 2 (ADU)",
        "ife2"  : "Flux error in aperture 2 (ADU)",
        "irm2"  : "Instrumental mag in aperture 2",
        "ire2"  : "Instrumental mag error for aperture 2",
        "irq2"  : "Instrumental quality flag ap 2, 0/G OK, X bad",
        "ifl3"  : "Flux in aperture 3 (ADU)",
        "ife3"  : "Flux error in aperture 3 (ADU)",
        "irm3"  : "Instrumental mag in aperture 3",
        "ire3"  : "Instrumental mag error for aperture 3",
        "irq3"  : "Instrumental quality flag ap 3, 0/G OK, X bad",
        "ccdtemp" : "mean CCD temperature S_CAM_ALCU_sensor_CCD",
        "ntemps"  : "number of temperatures avgd to get ccdtemp",
        'dtr_isub': "img subtraction photometry performed",
        'dtr_epd' : "EPD detrending performed",
        'dtr_tfa' : "TFA detrending performed"
    }
    return kcd[k]


def convert_grcollect_to_fits_lc_worker(task):
    """
    Given grcollect lightcurve, make a FITS binary table with the lightcurve
    data. Collect header information from the frames that created the
    lightcurve.
    """

    grclcpath, outpath, catdf, observatory, lcdir, fitsdir = task

    if observatory != 'tess':
        # if not TESS, you'll need to modify the header lists.
        raise NotImplementedError

    lcd = tu.read_tess_txt_lightcurve(grclcpath, catdf)

    times = lcd['tjd']
    earliestframename = lcd['rstfc'][np.argmin(times)]
    earliestframepath = os.path.join(fitsdir, earliestframename+'.fits')

    # read the fits header from the earliest frame in the sequence; inherit
    # various header properties from this frame.
    kwlist=['SIMPLE', 'SIMDATA', 'TELESCOP', 'INSTRUME', 'CAMERA', 'CCD',
            'EXPOSURE', 'TIMEREF', 'TASSIGN', 'TIMESYS', 'BJDREFI', 'BJDREFF',
            'TIMEUNIT', 'TELAPSE', 'LIVETIME', 'TSTART', 'TSTOP', 'DATE-OBS',
            'DATE-END', 'DEADC', 'TIMEPIXR', 'TIERRELA', 'BTC_PIX1',
            'BTC_PIX2', 'BARYCORR', 'INT_TIME', 'READ_TIME', 'FRAMETIM',
            'NUM_FRM', 'TIMEDEL', 'NREADOUT']

    hdict = iu.get_header_keyword_list(earliestframepath, kwlist, ext=0)
    cdict = iu.get_header_comment_list(earliestframepath, kwlist, ext=0)

    # make the primary header
    hdr = fits.Header()
    for k in np.sort(list(hdict.keys())):
        hdr[k.replace('_','')] = ( hdict[k], cdict[k] )

    primary_hdu = fits.PrimaryHDU(header=hdr)

    # make the lightcurve data extension
    datakeys = [k for k in np.sort(list(lcd.keys()))
                if k not in ['starid','xcc','ycc','objectid','objectinfo']
               ]

    outnames = [k.upper() for k in datakeys]

    formats = [_map_key_to_format(k) for k in datakeys]

    datacols = [lcd[k] for k in datakeys]

    fitscollist = [fits.Column(name=n, format=f, array=a) for n,f,a in
                   zip(outnames, formats, datacols)]

    # (TESS-specific): get the CCD temperature timeseries. This is stored
    # in each frame, as identified by the 'rstfc' data key.
    if observatory=='tess':

        fitsimglist = [os.path.join(fitsdir, fc+'.fits')
                       for fc in lcd['rstfc']]

        ccdtemplist, ntempslist = [], []
        for fpath in fitsimglist:
            ccdtemplist.append(iu.get_header_keyword(fpath, 'CCDTEMP'))
            ntempslist.append(iu.get_header_keyword(fpath, 'NTEMPS'))
        ccdtemp = nparr(ccdtemplist)
        ntemps = nparr(ntempslist)

        fitscollist.append(
            fits.Column(name='CCDTEMP', format='D', array=ccdtemp)
        )
        fitscollist.append(
            fits.Column(name='NTEMPS', format='J', array=ntemps)
        )

    hdutimeseries = fits.BinTableHDU.from_columns(fitscollist)

    # update the header comments for the timeseries data. to do this, since the
    # timeseries data is listed in sequential "TTYPE1", "TTYPE2", etc columns
    # (this is forced), need to map those column names to the keys first.
    tfields = hdutimeseries.header['TFIELDS'] # number of table fields
    tshdrkv = {}
    for ind in range(1,tfields+1):
        key = 'TTYPE{}'.format(ind)
        tshdrkv[key] = hdutimeseries.header[key]
    for k,v in tshdrkv.items():
        hdutimeseries.header.comments[k] = _map_key_to_comment(v.lower())

    # ditto, for the primary header data. include xcc and ycc. FITS does not
    # allow nans.
    lcd['objectinfo'] = {
        (k if not pd.isnull(v) else k):
        (v if not pd.isnull(v) else 'NaN')
        for k,v in lcd['objectinfo'].items()
    }

    for k in lcd['objectinfo'].keys():
        primary_hdu.header['HIERARCH '+k] = lcd['objectinfo'][k]

    primary_hdu.header['XCC'] = np.mean(lcd['xcc'])
    primary_hdu.header['YCC'] = np.mean(lcd['ycc'])
    primary_hdu.header['DTR_ISUB'] = True
    primary_hdu.header['DTR_EPD'] = False
    primary_hdu.header['DTR_TFA'] = False
    for k in ['xcc','ycc','DTR_ISUB','DTR_EPD','DTR_TFA']:
        primary_hdu.header.comments[k] = _map_key_to_comment(k.lower())

    hdulist = fits.HDUList([primary_hdu, hdutimeseries])

    hdulist.writeto(outpath, overwrite=True)
    print('wrote {}'.format(outpath))


def parallel_convert_grcollect_to_fits_lc(lcdirectory,
                                          fitsdir,
                                          catfile='/nfs/phtess1/ar1/TESS/FFI/RED_IMGSUB/TUNE/s0002/RED_4-4-1028_ISP/GAIADR2-RA98.4609411225734-DEC-58.5412164069738-SIZE12.catalog',
                                          ilcglob='*.grcollectilc',
                                          nworkers=16,
                                          maxworkertasks=1000,
                                          observatory='tess'):
    """
    Parallelizes convert_grcollect_to_fits_lc_worker.

    the TESS SPOC lc pattern is:
        tessyyyydddhhmmss-ssctr-tid-scid-cr_lc.fits.gz

    we will eventually use:
        projid1234_s000X_gaiaid_llc.fits.gz
    """

    # load the catalog file.
    if observatory=='tess':
        catdf = tu.read_object_catalog(catfile)
    else:
        raise NotImplementedError

    ilclist = glob(os.path.join(lcdirectory, ilcglob))
    outlist = [i.replace(os.path.splitext(ilcglob)[-1], '_llc.fits')
               for i in ilclist]

    print('%sZ: %s files to convert grcollect -> FITS lightcurves' %
          (datetime.utcnow().isoformat(), len(ilclist)))

    path_exists = []
    for outpath in outlist:
        if os.path.exists(outpath):
            path_exists.append(1)
        else:
            path_exists.append(0)

    path_exists = nparr(path_exists)

    if npall(path_exists):
        print(
            'found all {:d} FITS lightcurves, continuing'.
            format(len(path_exists))
        )
        return 0

    else:

        tasks = [(x, y, catdf, observatory, lcdirectory, fitsdir) for x,y in
                 zip(ilclist, outlist)]

        pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)
        results = pool.map(convert_grcollect_to_fits_lc_worker, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()

        return 1


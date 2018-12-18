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

##########################################

def from_CAL_to_fitsh_compatible(fits_list, outdir):
    '''
    This function takes a list of calibrated images from MAST and turns them
    into fitsh-compatible images.

    This means:
        * single extension
        * trimmed to remove virtual columns

    Args:
        fits_list (np.ndarray): list of calibrated full frame images from MAST.
        (Each image is a single CCD).

        outdir (str): directory to write the files to.

    Returns:
        nothing.
    '''

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    path_exists = []
    for ix, fits_name in enumerate(fits_list):
        outname = (
            outdir + os.path.basename(
                fits_name.replace('-s_ffic.fits', '_cal_img.fits'))
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

    for ix, fits_name in enumerate(fits_list):

        outname = (
            outdir + os.path.basename(
                fits_name.replace('-s_ffic.fits', '_cal_img.fits'))
        )

        if os.path.exists(outname):
            print(
                '{:d}/{:d}. Found {:s}, continuing.'.format(
                    ix, len(fits_list), outname
                )
            )
            continue

        data, hdr = iu.read_fits(fits_name, ext=1)

        # FITS header is 1-based counting, but python is 0-based. To convert
        # from FITS->python from slicing you need to subtract 1.
        rows = hdr['SCIROWS']-1 # row start.
        rowe = hdr['SCIROWE']-1+1 # row end (inclusive)
        cols = hdr['SCCSA']-1 # col start
        cole = hdr['SCCED']-1+1 # col end (inclusive)

        trim = data[slice(rows,rowe),slice(cols,cole)]

        assert trim.shape == (2048, 2048)

        fits.writeto(outname, trim, header=hdr)

        print(
            '{:d}/{:d}. Wrote to {:s} to {:s}'.format(
                ix, len(fits_list), fits_name, outname
            )
        )


def from_ete6_to_fitsh_compatible(fits_list, outdir):
    '''
    This function takes a list of ETE6 reduced images and turns them into
    fitsh-compatible images.
    '''

    from_CAL_to_fitsh_compatible(fits_list, outdir)

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




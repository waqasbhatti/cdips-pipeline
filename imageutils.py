#!/usr/bin/env python

'''
imageutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Jan 2013

This contains various utilities for operating on images. This includes
generating stamps for an image, converting an image to JPEGs, and getting
values of certain keywords from FITS headers.

====================
fits-reading functions:
    read_fits
    read_fits_header
    trim_image
    make_superflat
    compressed_fits_ext
    get_header_keyword
    get_header_keyword_list
    get_header_comment_list
    get_data_keyword_list

image-scaling functions:
    zscale_img
    clipped_linscale_img
    logscale_img
    clipped_logscale_img
    extract_img_background

image-sectioning and image-writing functions:
    mplplot_logscale_img_w_colorbar
    mplplot_diffscale_img_w_colorbar
    img_to_stamps
    stamps_background
    stamps_to_jpeg
    fits_to_stamps_jpeg
    fits_to_full_jpeg: make a jpg from a fits image
    frame_radecbox_to_jpeg: cuts out box centered at RA/DEC and width
    fitscoords_to_jpeg
    nparr_to_full_jpeg
    check_frame_warping

movie-making functions:
    make_mp4_from_jpegs
    make_mov_from_jpegs

image-processing diagnostic plots:

====================
'''

import os
import os.path
import sys
import logging
from glob import glob
from datetime import datetime

import numpy as np
np.seterr(all='ignore')

import numpy.ma as npma
import numpy.random as npr

import scipy.misc
import scipy.ndimage
import scipy

from scipy.optimize import leastsq
USE_LEASTSQ = 1

try:
    from scipy.optimize import curve_fit
    USE_LEASTSQ=0
except:
    print('cannot import curve_fit, will use leastsq')
    USE_LEASTSQ=1


import astropy.io.fits as pyfits
from astropy import wcs

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects

from astrobase import imageutils as iu

# get the ImageFont
fontpath = os.path.join(os.path.dirname(__file__), 'DejaVuSans.ttf')

# load the font
if os.path.exists(fontpath):
    fontxxsmall = ImageFont.truetype(fontpath, 4)
    fontxsmall = ImageFont.truetype(fontpath, 8)
    fontsmall = ImageFont.truetype(fontpath, 12)
    fontnormal = ImageFont.truetype(fontpath, 20)
    fontlarge = ImageFont.truetype(fontpath, 28)
else:
    print('could not find bundled '
          'DejaVu Sans font, using ugly defaults...')
    fontsmall = ImageFont.load_default()
    fontnormal = ImageFont.load_default()
    fontlarge = ImageFont.load_default()


# setup a logger
LOGGER = None

def set_logger_parent(parent_name):
    globals()['LOGGER'] = logging.getLogger('%s.imageutils' % parent_name)


## FITS UTILITIES ##


def read_fits(fits_file,ext=0):
    '''
    Shortcut function to get the header and data from a fits file and a given
    extension.

    '''

    hdulist = pyfits.open(fits_file)
    img_header = hdulist[ext].header
    img_data = hdulist[ext].data
    hdulist.close()

    return img_data, img_header


def read_fits_header(fits_file, ext=0):
    '''
    Shortcut function to just read the header of the FITS file and return it.

    '''
    hdulist = pyfits.open(fits_file)
    img_header = hdulist[ext].header
    hdulist.close()

    return img_header



def trim_image(fits_img,
               fits_hdr,
               custombox=None):
    '''
    Returns a trimmed image using the TRIMSEC header of the image header.
    FIXME: check if this does the right thing.

    custombox is a string of the form [Xlo:Xhi,Ylo:Yhi] and will trim the image
    to a custom size.

    '''

    if custombox:

        trimsec = custombox

    else:

        if 'TRIMSEC' in fits_hdr:
            trimsec = fits_hdr['TRIMSEC']
        elif 'DATASEC' in fits_hdr:
            trimsec = fits_hdr['DATASEC']
        else:
            if custombox is None:
                if LOGGER:
                    LOGGER.error('no DATASEC or TRIMSEC in image header')
                else:
                    print('cannot trim image, no DATASEC or '
                          'TRIMSEC in image header')
                return

    if trimsec != '[0:0,0:0]':

        datasec = trimsec.strip('[]').split(',')

        try:
            datasec_y = [int(x) for x in datasec[0].split(':')]
            datasec_x = [int(x) for x in datasec[1].split(':')]

            trimmed_img = fits_img[datasec_x[0]-1:datasec_x[1],
                                   datasec_y[0]-1:datasec_y[1]]
        except ValueError as e:
            if LOGGER:
                LOGGER.error('datasec/trimsec not correctly set in FITS header, '
                             ' not trimming')
            else:
                print('datasec/trimsec not correctly set in FITS header, '
                      ' not trimming')
            trimmed_img = fits_img

    else:
        if LOGGER:
            LOGGER.error('datasec/trimsec not correctly set in FITS header, '
                         ' not trimming')
        else:
            print('datasec/trimsec not correctly set in FITS header, '
                  ' not trimming')
        trimmed_img = fits_img

    return trimmed_img



def make_superflat(image_glob,
                   fits_imagetype_card = 'IMAGETYP',
                   fits_flat_keyword='flat',
                   smoothlevel=11,
                   ext=None,
                   method='mean',
                   saveto=None):
    '''
    This generates a normalized superflat image for a series of flatfield
    images.

    1. finds all flat field images in image_glob
    2. takes their average
    3. normalizes by dividing out the median value (optional)
    4. smooths the flatfield image so that small scale problems still show up
       when this flat field is divided from the object frames (optional)

    '''

    image_flist = sorted(glob(image_glob))

    # go through the images and find all the flats

    flat_imgs = {}
    flat_count = 0

    for fits_image in image_flist:

        compressed_ext = compressed_fits_ext(fits_image)

        if ext is None and compressed_ext:
            img, hdr = read_fits(fits_image,
                                 ext=compressed_ext[0])
        elif (ext is not None):
            img, hdr = read_fits(fits_image,ext=ext)
        else:
            img, hdr = read_fits(fits_image)

        if hdr[fits_imagetype_card] == fits_flat_keyword:
            trimmed_img = trim_image(img, hdr)
            flat_imgs[fits_image] = trimmed_img
            flat_count = flat_count + 1
            print('found flat %s' % fits_image)

    if flat_count > 1:

        all_flats = np.asarray([flat_imgs[k] for k in flat_imgs])
        del flat_imgs

        # get the median/mean of the flats depending on method
        if method == 'mean':
            median_flat = np.mean(all_flats,axis=0)
        elif method == 'median':
            median_flat = np.median(all_flats,axis=0)

        smoothed_flat = scipy.ndimage.median_filter(median_flat,
                                                    size=smoothlevel)

        if saveto:
            pyfits.writeto(saveto,smoothed_flat)
        else:
            return smoothed_flat

    else:

        return None



def compressed_fits_ext(fits_file):
    '''
    Check if a fits file is a compressed FITS file. Return the extension numbers
    of the compressed image as a list if these exist, otherwise, return None.

    '''

    hdulist = pyfits.open(fits_file)

    compressed_img_exts = []

    for i, ext in enumerate(hdulist):
        if isinstance(ext,pyfits.hdu.compressed.CompImageHDU):
            compressed_img_exts.append(i)

    hdulist.close()

    if len(compressed_img_exts) < 1:
        return None
    else:
        return compressed_img_exts


def get_header_keyword(fits_file,
                       keyword,
                       ext=0):

    return iu.get_header_keyword(fits_file, keyword, ext=ext)


def get_data_keyword(fits_file,
                     keyword,
                     ext=1):

    return iu.get_data_keyword(fits_file, keyword, ext=ext)


def get_header_keyword_list(fits_file,
                            keyword_list,
                            ext=0):

    return iu.get_header_keyword_list(fits_file, keyword_list, ext=ext)


def get_header_comment_list(fits_file,
                            keyword_list,
                            ext=0):

    return iu.get_header_comment_list(fits_file, keyword_list, ext=ext)


def get_data_keyword_list(fits_file,
                          keyword_list,
                          ext=1):

    return iu.get_data_keyword_list(fits_file, keyword_list, ext=ext)


## IMAGE SCALING FUNCTIONS ##

def pixel_scale_func(x, m, c):
    return m*x + c


def pixel_scale_func_residual(params, x, y):

    f = pixel_scale_func(x, params[0], params[1])
    return y - f


def zscale_img(img_array,
               cap=255.0,
               fracsamples=0.1):
    '''
    This scales the image pixels in a manner similar to what DS9 does when
    zscale and linear are selected in the scale menu.

    Algorithm found here:

    http://iraf.net/phpBB2/viewtopic.php?t=77998&sid=b5ee7df81074f31fa7086aa1f31a74be

    Quoting the second comment from there:

    - sample the image (1000 points or so depending on size) in a grid covering
      the full frame to get a representative sample of all pixels in the image.

    - sort the sample pixels to get min/max/median values

    - iteratively fit a line to map the sample data to the number of pixels you
      want on output (e.g. 256 for 8-bit data, ximtool uses 200 or so for the
      display). Some of the sample pixels are usually rejected at each iteration
      to get a better fit.

    - from the fitted slope derive the optimal z1/z2 end values for the
      data. When mapping the image, input values outside this range maps to the
      extremes of your output range, everything in between maps linearly. The
      brightness/contrast adjustments are done by changing the offset and slope
      of this linear transformation respectively.  In display servers this
      usually means just rewriting the colormap for the image but it could also
      be used to remap all the image pixels.

      nsamples = fraction of total pixels to use for statistics
      cap = fix the max value to be within range 0-255

    '''

    img_shape = img_array.shape

    total_pixels = img_shape[0]*img_shape[1]
    nsamples = int(np.floor(fracsamples*total_pixels))

    random_index_x = npr.random_integers(0,high=img_shape[1]-1,size=nsamples)
    random_index_y = npr.random_integers(0,high=img_shape[0]-1,size=nsamples)

    # the x values
    img_sample = img_array[random_index_x, random_index_y]

    sample_med, sample_min, sample_max = (np.median(img_sample),
                                          np.nanmin(img_sample),
                                          np.nanmax(img_sample))
    sample_std = np.std(img_sample)

    trimmed_sample_ind = np.where(abs(img_sample - sample_med) < 1.0*sample_std)
    trimmed_sample = img_sample[trimmed_sample_ind]

    trimmed_sample = np.sort(trimmed_sample)

    # the y values: we're mapping our img_sample to a range between 0 and cap
    pixel_scale = np.linspace(0, cap, num=len(trimmed_sample))

    initial_slope = np.median(pixel_scale/trimmed_sample)
    initial_intercept = (np.median(pixel_scale) -
                         initial_slope*np.median(trimmed_sample))

    if USE_LEASTSQ == 1:

        params = leastsq(pixel_scale_func_residual,
                               np.array([initial_slope,
                                         initial_intercept]),
                               args=(trimmed_sample, pixel_scale))

        scale_params = params[0]

    else:
        scale_params, scale_covariance = curve_fit(pixel_scale_func,
                                                   trimmed_sample,
                                                   pixel_scale,
                                                   p0=(initial_slope,
                                                       initial_intercept))

    sample_med, sample_min, sample_max = (np.median(trimmed_sample),
                                          np.nanmin(trimmed_sample),
                                          np.nanmax(trimmed_sample))

    min_scale_param = sample_min*scale_params[0] + scale_params[1]
    max_scale_param = sample_max*scale_params[0] + scale_params[1]

    print(min_scale_param,
          max_scale_param)

    print(np.min(img_array), np.max(img_array))

    clipped_image_array = np.clip(img_array, min_scale_param, max_scale_param)

    return scale_params[0]*clipped_image_array + scale_params[1]


def clipped_linscale_img(img_array,
                         cap=255.0,
                         lomult=2.0,
                         himult=2.0):
    '''
    This clips the image between the values:

    [median(img_array) - lomult*stdev(img_array),
     median(img_array) + himult*stdev(img_array)]

    and returns a linearly scaled image using the cap given.

    '''

    img_med, img_stdev = np.median(img_array), np.std(img_array)
    clipped_linear_img = np.clip(img_array,
                                 img_med-lomult*img_stdev,
                                 img_med+himult*img_stdev)
    return cap*clipped_linear_img/(img_med+himult*img_stdev)


def logscale_img(img_array,
                 cap=255.0,
                 coeff=1000.0):
    '''
    This scales the image according to the relation:

    logscale_img = np.log(coeff*(img/max(img))+1)/np.log(coeff)

    Taken from the DS9 scaling algorithms page at:

    http://hea-www.harvard.edu/RD/ds9/ref/how.html

    According to that page:

    coeff = 1000.0 works well for optical images
    coeff = 100.0 works well for IR images

    '''

    logscaled_img = np.log(coeff*img_array/np.nanmax(img_array)+1)/np.log(coeff)
    return cap*logscaled_img


def clipped_logscale_img(img_array,
                         cap=255.0,
                         lomult=2.0,
                         himult=2.0,
                         loclip=None,
                         hiclip=None,
                         coeff=1000.0):
    '''
    This clips the image between values, and then log-scales it. If lomult and
    himult are passed, the clipping happens between

        [median(img_array) - lomult*stdev(img_array),
         median(img_array) + himult*stdev(img_array)]

    else if loclip and hiclip are passed, the clipping happens between

        [loclip, hiclip].

    The log-scaled image is returned using the cap given.

    logscale_img = np.log(coeff*(img/max(img))+1)/np.log(coeff)
    '''

    img_med, img_stdev = np.median(img_array), np.std(img_array)
    if isinstance(lomult, float) and isinstance(himult, float):
        clipped_linear_img = np.clip(img_array,
                                     img_med-lomult*img_stdev,
                                     img_med+himult*img_stdev)
    elif isinstance(loclip, float) and isinstance(hiclip, float):
        clipped_linear_img = np.clip(img_array, loclip, hiclip)
    else:
        raise AssertionError(
            'expected either (lomult,himult), or (loclip,hiclip) to be passed'
        )

    clipped_linear_img = clipped_linear_img/(img_med+himult*img_stdev)

    # janky
    clipped_linear_img[clipped_linear_img<0] = np.nan

    div = np.nanmax(clipped_linear_img)

    logscaled_img = (
            np.log(coeff*clipped_linear_img/div+1)
            /
            np.log(coeff)
    )

    return cap*logscaled_img


def extract_img_background(img_array,
                           custom_limits=None,
                           median_diffbelow=200.0,
                           image_min=None):
    '''
    This extracts the background of the image array provided:

    - masks the array to only values between the median and the min of flux
    - then returns the median value in 3 x 3 stamps.

    img_array = image to find the background for

    custom_limits = use this to provide custom median and min limits for the
                    background extraction

    median_diffbelow = subtract this value from the median to get the upper
                       bound for background extraction

    image_min = use this value as the lower bound for background extraction

    '''

    if not custom_limits:

        backmax = np.median(img_array)-median_diffbelow
        backmin = image_min if image_min is not None else np.nanmin(img_array)

    else:

        backmin, backmax = custom_limits

    masked = npma.masked_outside(img_array, backmin, backmax)
    backmasked = npma.median(masked)

    return backmasked


## IMAGE SECTION FUNCTIONS ##

def mplplot_logscale_img_w_colorbar(
    img,
    outpath,
    vmin=10, vmax=int(1e3),
    cmap='binary_r',
    titlestr=None):

    if os.path.exists(outpath):
        print('found {}. continue'.format(outpath))
        return 0

    plt.close('all')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, ax = plt.subplots(figsize=(6,4.5))

    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    cset1 = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    #ax.set_xticklabels('')
    #ax.set_yticklabels('')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    #ax.xaxis.set_ticks_position('none')
    #ax.yaxis.set_ticks_position('none')

    cb1 = fig.colorbar(cset1, ax=ax, extend='both')
    #cb2.set_ticks([-1e3,-1e2,-1e1,0,1e1,1e2,1e3])
    #cb2.set_ticklabels(['-$10^3$','-$10^2$','-$10^1$','0',
    #                    '$10^1$','$10^2$','$10^3$'])

    if isinstance(titlestr, str):
        ax.set_title(titlestr, fontsize='x-small')

    fig.tight_layout(pad=0)

    fig.savefig(outpath, bbox_inches='tight', dpi=300)
    print('{}: made {}'.format(datetime.utcnow().isoformat(), outpath))


def mplplot_diffscale_img_w_colorbar(
    img,
    outpath,
    vmin=-1000, vmax=1000,
    cmap='RdBu_r',
    titlestr=None):

    if os.path.exists(outpath):
        print('found {}. continue'.format(outpath))
        return 0

    plt.close('all')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, ax = plt.subplots(figsize=(6,4.5))

    diffnorm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=vmin,
                                 vmax=vmax)

    cset2 = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, norm=diffnorm)

    #ax.set_xticklabels('')
    #ax.set_yticklabels('')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    #ax.xaxis.set_ticks_position('none')
    #ax.yaxis.set_ticks_position('none')

    cb2 = fig.colorbar(cset2, ax=ax, extend='both')

    cb2.set_ticks([-1e3,-1e2,-1e1,0,1e1,1e2,1e3])
    cb2.set_ticklabels(['-$10^3$','-$10^2$','-$10^1$','0',
                        '$10^1$','$10^2$','$10^3$'])

    if isinstance(titlestr, str):
        ax.set_title(titlestr, fontsize='x-small')

    fig.tight_layout(pad=0)

    fig.savefig(outpath, bbox_inches='tight', dpi=300)
    print('{}: made {}'.format(datetime.utcnow().isoformat(), outpath))



def img_to_stamps(img,
                  stampsize=256):
    '''
    Generate stamps for an image of size imgsizex x imgsize y. Stamps not in the
    center of the image will be generated for the edges of the image. This
    generates 3 x 3 stamps for each image.

    top_left_corner = img[:xstampsize,:ystampsize]
    bottom_right_corner = img[-xstampsize:,-ystampsize:]

    top_right_corner = img[imgsizex-xstampsize:,:ystampsize]
    bottom_left_corner = img[:xstampsize,imgsizey-ystampsize:]

    top_center = img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,:ystampsize]
    bottom_center = img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                        imgsizey-ystampsize:]

    center = img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                 imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2]

    right_center = img[imgsizex-xstampsize:,
                       imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2]
    left_center = img[:xstampsize,imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2]


    '''

    imgsizex, imgsizey = img.shape
    xstampsize, ystampsize = stampsize, stampsize

    # get the total number of possible stamps
    n_possible_xstamps = imgsizex/float(xstampsize)
    n_possible_ystamps = imgsizey/float(ystampsize)


    # if we can actually make stamps, then go ahead
    if (n_possible_xstamps >= 3) and (n_possible_ystamps >= 3):

        # FIXME: the coordinate slices should be swapped here, i.e. x,y -> y,x,
        # because the np indexing scheme is y,x instead of x,y
        return {'topleft':img[:xstampsize,:ystampsize],
                'topcenter':img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                                :ystampsize],
                'topright':img[imgsizex-xstampsize:,:ystampsize],
                'midleft':img[:xstampsize,
                               imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2],
                'midcenter':img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                                imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2],
                'midright':img[imgsizex-xstampsize:,
                               imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2],
                'bottomleft':img[:xstampsize,imgsizey-ystampsize:],
                'bottomcenter':img[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                                   imgsizey-ystampsize:],
                'bottomright':img[-xstampsize:,-ystampsize:]}
    else:
        if LOGGER:
            LOGGER.error('stampsize is too large for this image')
        else:
            print('error: stampsize is too large for this image')
        return None


def stamps_background(image_stamps,
                      custom_limits=None,
                      median_diffbelow=200.0,
                      image_min=None):
    '''
    This returns background values for each of the stamps in the image_stamps
    object, using the extract_img_background function above.

    '''

    return dict(
        [(
                key,extract_img_background(
                    image_stamps[key],
                    custom_limits=custom_limits,
                    median_diffbelow=median_diffbelow,
                    image_min=image_min
                    )
                )
         for key in image_stamps]
        )


def stamps_to_jpeg(image_stamps,
                   out_fname,
                   sepwidth=1,
                   scale=False,
                   scale_func=clipped_linscale_img,
                   scale_func_params={'cap':255.0,
                                      'lomult':2,
                                      'himult':2.5}):
    '''
    This turns the stamps returned from the function img_to_stamps above into
    a single 3 x 3 postage stamp image. Uses sepwidth pixels as the separator
    between each row/line of stamps.

    '''

    toprow_xsize, toprow_ysize = image_stamps['topright'].shape
    toprow_separr = np.array([[255.0]*sepwidth]*toprow_ysize)

    # note, these should be topleft, topcenter, topright, but since a[x,y] is
    # actually a[y,x] in np array coordinates, it is backwards.
    # FIXME: fix this by fixing img_to_stamps above

    # get the stamps
    if scale:

        topleft = scale_func(image_stamps['topleft'],
                             **scale_func_params)
        midleft = scale_func(image_stamps['midleft'],
                             **scale_func_params)
        bottomleft = scale_func(image_stamps['bottomleft'],
                             **scale_func_params)

        topcenter = scale_func(image_stamps['topcenter'],
                             **scale_func_params)
        midcenter = scale_func(image_stamps['midcenter'],
                             **scale_func_params)
        bottomcenter = scale_func(image_stamps['bottomcenter'],
                             **scale_func_params)

        topright = scale_func(image_stamps['topright'],
                             **scale_func_params)
        midright = scale_func(image_stamps['midright'],
                             **scale_func_params)
        bottomright = scale_func(image_stamps['bottomright'],
                             **scale_func_params)

    else:

        topleft = image_stamps['topleft']
        midleft = image_stamps['midleft']
        bottomleft = image_stamps['bottomleft']

        topcenter = image_stamps['topcenter']
        midcenter = image_stamps['midcenter']
        bottomcenter = image_stamps['bottomcenter']

        topright = image_stamps['topright']
        midright = image_stamps['midright']
        bottomright = image_stamps['bottomright']


    toprow_stamp = np.hstack((topleft,
                              toprow_separr,
                              midleft,
                              toprow_separr,
                              bottomleft))

    midrow_xsize, midrow_ysize = midright.shape
    midrow_separr = np.array([[255.0]*sepwidth]*midrow_ysize)

    # similarly, these should be midleft, midcenter, midright
    midrow_stamp = np.hstack((topcenter,
                              midrow_separr,
                              midcenter,
                              midrow_separr,
                              bottomcenter))

    bottomrow_xsize, bottomrow_ysize = bottomright.shape
    bottomrow_ysize = bottomright.shape[1]
    bottomrow_separr = np.array([[255.0]*sepwidth]*bottomrow_ysize)

    # similarly, these should be bottomleft, bottomcenter, bottomright
    bottomrow_stamp = np.hstack((topright,
                                 bottomrow_separr,
                                 midright,
                                 bottomrow_separr,
                                 bottomright))

    full_stamp = np.vstack((toprow_stamp,
                            np.array([255.0]*(toprow_xsize*3 + sepwidth*2)),
                            midrow_stamp,
                            np.array([255.0]*(midrow_xsize*3 + sepwidth*2)),
                            bottomrow_stamp))

    scipy.misc.imsave(out_fname,full_stamp)
    return full_stamp


def fits_to_stamps_jpeg(fits_image,
                        out_fname=None,
                        ext=None,
                        stampsize=256,
                        sepwidth=1,
                        scale_func=clipped_linscale_img,
                        scale_func_params={'cap':255.0,
                                           'lomult':2,
                                           'himult':2.5}):
    '''
    This turns a FITS image into a 3 x 3 stamps JPEG.

    '''

    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    trimmed_img = trim_image(img, hdr)
    scaled_img = scale_func(trimmed_img,**scale_func_params)
    stamps = img_to_stamps(scaled_img,stampsize=stampsize)
    if out_fname is None:
        out_fname = fits_image + '.stamp.jpeg'
    stamps_img = stamps_to_jpeg(stamps,out_fname,sepwidth=sepwidth)



def fits_to_full_jpeg(fits_image,
                      out_fname=None,
                      ext=None,
                      resize=False,
                      flip=True,
                      outsizex=800,
                      outsizey=800,
                      annotate=True,
                      fits_jdsrc=None,
                      scale_func=clipped_linscale_img,
                      scale_func_params={'cap':255.0,
                                         'lomult':2,
                                         'himult':2.5},
                      frame_time=None,
                      colorscheme=None):
    '''
    This converts a FITS image to a full frame JPEG.

    kwargs:

        scale_func (function): clipped_linscale_img, clipped_logscale_img

        colorscheme (None or str): name of matplotlib colorscheme to use. e.g.,
        "bwr" looks good for subtracted images.
    '''
    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    #trimmed_img = trim_image(img, hdr)
    trimmed_img = img
    jpegaspect = float(img.shape[1])/float(img.shape[0])
    scaled_img = scale_func(trimmed_img,**scale_func_params)

    if resize:
        resized_img = scipy.misc.imresize(scaled_img,
                                          (int(img.shape[1]/2.2),
                                           int(img.shape[0]/2.2)))
    else:
        resized_img = scaled_img

    if not out_fname:

        out_fname = '%s-%s-%s-%s-proj%s-%s.jpg' % (
            fits_image.rstrip('.fits.fz'),
            hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
            hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
            (hdr['FILTERS'].replace('+','') if
             'FILTERS' in hdr else 'filtunknown'),
            hdr['PROJID'] if 'PROJID' in hdr else 'unknown',
            hdr['OBJECT'] if 'OBJECT' in hdr else 'objectunknown'
            )

    scipy.misc.imsave(out_fname,resized_img)

    # recolor the saved image if told to do so
    if colorscheme:
        cm = mplcm.get_cmap(colorscheme)
        outimg = Image.open(out_fname)
        im = np.array(outimg)
        im = cm(im)
        im = np.uint8( im*255.0 )
        outimg = Image.fromarray(im)
        rgb_outimg = outimg.convert('RGB')
        rgb_outimg.save(out_fname)

    # flip the saved image if told to do so
    if flip:
        outimg = Image.open(out_fname)
        outimg = outimg.transpose(Image.FLIP_TOP_BOTTOM)
        outimg.save(out_fname)

    # annotate the image if told to do so
    if annotate:

        outimg = Image.open(out_fname)
        draw = ImageDraw.Draw(outimg)
        annotation = "%s: %s - %s - %s - PR%s - %s" % (
            os.path.basename(fits_image).rstrip('.fits.fz'),
            hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
            hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
            (hdr['FILTERS'].replace('+','') if
             'FILTERS' in hdr else 'filtunknown'),
            hdr['PROJID'] if 'PROJID' in hdr else 'unknown',
            hdr['OBJECT'] if 'OBJECT' in hdr else 'objectunknown'
        )
        draw.text((10,10),
                  annotation,
                  font=fontnormal,
                  fill=255)

        # now add the time as well

        # if we're supposed to use another file for the JD source, do so
        # this is useful for subtracted images
        if fits_jdsrc is not None and os.path.exists(fits_jdsrc):
            framejd = get_header_keyword(fits_jdsrc, 'JD')
        elif frame_time is not None:
            framejd = frame_time
        else:
            framejd = hdr['JD'] if 'JD' in hdr else None

        if framejd is not None:
            timeannotation = '%.5f' % framejd
            draw.text((10, resized_img.shape[1] - 40),
                      timeannotation,
                      font=fontlarge,
                      fill=255)

        del draw
        outimg.save(out_fname)

    return out_fname


def _given_radecbox_get_xybox(wcsfrom, fits_image, radecbox, radeccenter, hdr,
                              img, do_spoc_trim_shift=False, forcesquare=True,
                              verbose=False):
    """
    Stand-alone helper to frame_radecbox_to_jpeg. requires a FITS image, its
    WCS, and the coords of a box,

        `radeccenter = [ra, dec, boxwidth, boxheight]`.

    This function is relevant for any image box-trimming though.

    args:
        wcsfrom: wcs file
        fits_image: path to fits image
        hdr, img: from reading the fits image

    returns:
        xmin,xmax,ymin,ymax to trim the image.
    """

    try:
        # get the WCS header
        if wcsfrom and os.path.exists(wcsfrom):
            w = wcs.WCS(wcsfrom)
        else:
            w = wcs.WCS(fits_image)
    except:
        print('no WCS found!')
        w = None

    # convert the radecbox into a pixbox
    if w and radecbox and not radeccenter:

        rd = np.array([[radecbox[0], radecbox[2]],
                       [radecbox[1], radecbox[3]]])

        # we use 0 here for the origin because we'll be cutting using np.arrays
        if verbose:
            print('requested coords = %s' % repr(rd))
        pix = w.all_world2pix(rd,0)

    # otherwise, convert the radeccenter into pixcenter
    elif w and radeccenter and not radecbox:

        rd = np.array(
            [
                [radeccenter[0] - (radeccenter[2])/2.0,
                 radeccenter[1] - (radeccenter[3])/2.0],
                [radeccenter[0] + (radeccenter[2])/2.0,
                 radeccenter[1] + (radeccenter[3])/2.0],

            ]
        )

        if verbose:
            print('requested coords = %s' % repr(rd))
        pix = w.all_world2pix(rd,0)

    else:

        if not w:
            print("no suitable WCS found")
        else:
            print("can't specify both radeccenter and "
                  "radecbox at the same time")
        return None

    # do the cutout using a box generated by the radec -> pix bits above
    x1, x2, y1, y2 = pix[0,0], pix[1,0], pix[0,1], pix[1,1]

    if do_spoc_trim_shift:

        x1 -= (hdr['SCCSA']-1)
        y1 -= (hdr['SCIROWS']-1)
        x1 -= 0.5
        y1 -= 0.5

        x2 -= (hdr['SCCSA']-1)
        y2 -= (hdr['SCIROWS']-1)
        x2 -= 0.5
        y2 -= 0.5

    # figure out xmin, xmax, ymin, ymax
    if x1 > x2:
        xmin = x2
        xmax = x1
    else:
        xmin = x1
        xmax = x2

    if y1 > y2:
        ymin = y2
        ymax = y1
    else:
        ymin = y1
        ymax = y2

    # round the pix coords to integers
    xmin, xmax = int(np.round(xmin)), int(np.round(xmax))
    ymin, ymax = int(np.round(ymin)), int(np.round(ymax))

    # make sure we take care of edges
    if xmin < 0:
        xmin = 0
    if xmax >= img.shape[1]:
        xmax = img.shape[1] - 1
    if ymin < 0:
        ymin = 0
    if ymax >= img.shape[0]:
        ymax = img.shape[0] - 1

    if forcesquare:

        ydelta = ymax-ymin
        xdelta = xmax-xmin
        sqdelta = max((xdelta, ydelta))

        ymid = ymin + ydelta/2
        xmid = xmin + xdelta/2

        ymin, ymax = (int(np.round(ymid - sqdelta/2)),
                      int(np.round(ymid + sqdelta/2)))
        xmin, xmax = (int(np.round(xmid - sqdelta/2)),
                      int(np.round(xmid + sqdelta/2)))

    return xmin, xmax, ymin, ymax



def frame_radecbox_to_jpeg(
        fits_image,
        wcsfrom=None,
        radecbox=None,
        radeccenter=None,
        out_fname=None,
        ext=None,
        flip=True,
        annotatejd=True,
        annotate=True,
        jdsrc=None,
        forcesquare=False,
        overplotscalebar=False,
        rescaleimage=False,
        scale_func=clipped_linscale_img,
        scale_func_params={'cap':255.0,
                           'lomult':2,
                           'himult':2.5},
        colorscheme=None,
        verbose=True,
        do_spoc_trim_shift=False):
    '''This cuts out a box centered at RA/DEC and width from the FITS to JPEG.

    wcsfrom indicates that the frame WCS should be taken from the specified file
    (usually a .wcs in our pipeline).

    if radecbox and not radeccenter:
        radecbox = [rmin, rmax, dmin, dmax] of box to cut out of FITS

    elif radeccenter and not radecbox:
        radeccenter = [rcenter, dcenter, rwidth, dwidth]

    else:
        do nothing, since we can't have both at the same time

    Other options:

        forcesquare (bool): forces output image to be square.

        overplotscalebar (bool):

        rescaleimage (bool):
    '''
    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    #trimmed_img = trim_image(img, hdr)
    trimmed_img = img
    jpegaspect = float(img.shape[1])/float(img.shape[0])

    xmin, xmax, ymin, ymax = _given_radecbox_get_xybox(
        wcsfrom, fits_image,
        radecbox, radeccenter,
        hdr, img,
        do_spoc_trim_shift=do_spoc_trim_shift, forcesquare=forcesquare,
        verbose=verbose
    )

    # numpy is y,x so make sure to reverse the order
    trimmed_img = trimmed_img[ymin:ymax, xmin:xmax]

    # do the scaling after the image has been cut so it's right for the objects
    # in the cutout
    scaled_img = scale_func(trimmed_img,**scale_func_params)

    if not out_fname:

        out_fname = '%s-%s-%s-%s-proj%s-%s.jpg' % (
            fits_image.rstrip('.fits.fz'),
            hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
            hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
            (hdr['FILTERS'].replace('+','') if
             'FILTERS' in hdr else 'filtunknown'),
            hdr['PROJID'] if 'PROJID' in hdr else 'unknown',
            hdr['OBJECT'] if 'OBJECT' in hdr else 'objectunknown'
            )

        if radecbox and not radeccenter:
            out_fname = '%s-R%sR%s-D%sD%s.jpg' % (
                out_fname.rstrip('.jpg'),
                radecbox[0], radecbox[1],
                radecbox[2], radecbox[3]
            )

        elif radeccenter and not radecbox:
            out_fname = '%s-RC%sDC%s-RW%sDW%s.jpg' % (
                out_fname.rstrip('.jpg'),
                radeccenter[0], radeccenter[1],
                radeccenter[2], radeccenter[3]
            )

    if flip:
        scaled_img = np.flipud(scaled_img)

    scipy.misc.imsave(out_fname, scaled_img)

    if colorscheme:
        cm = mplcm.get_cmap(colorscheme)
        outimg = Image.open(out_fname)
        im = np.array(outimg)
        im = cm(im)
        im = np.uint8( im*255.0 )
        outimg = Image.fromarray(im)
        rgb_outimg = outimg.convert('RGB')
        rgb_outimg.save(out_fname)

    # annotate the image if told to do so
    if annotatejd and jdsrc and os.path.exists(jdsrc):

        # get the JD header keyword from jdsrc
        framejd = get_header_keyword(jdsrc, 'JD')

        outimg = Image.open(out_fname)
        draw = ImageDraw.Draw(outimg)
        annotation = "JD %.3f" % framejd
        draw.text((4,2),
                  annotation,
                  fill=255,
                  font=fontsmall)

        del draw
        outimg.save(out_fname)

    if annotate:
        outimg = Image.open(out_fname)
        draw = ImageDraw.Draw(outimg)
        if not isinstance(annotate, str):
            annotation = "%s: %s - %s - %s - PR%s - %s" % (
                os.path.basename(fits_image).rstrip('.fits.fz'),
                hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
                hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
                (hdr['FILTERS'].replace('+','') if
                 'FILTERS' in hdr else 'filtunknown'),
                hdr['PROJID'] if 'PROJID' in hdr else 'unknown',
                hdr['OBJECT'] if 'OBJECT' in hdr else 'objectunknown'
            )
        else:
            annotation = annotate
        draw.text((10,10),
                  annotation,
                  font=fontxsmall,
                  fill=255)

        del draw
        outimg.save(out_fname)

    if overplotscalebar:
        outimg = Image.open(out_fname)
        draw = ImageDraw.Draw(outimg)

        linelength = 15 # pixels. for TESS-> ~=5 arcminutes.
        if not forcesquare:
            raise AssertionError
        refpx = int(0.95*np.array(outimg.size)[0])

        x1, x2 = refpx-linelength, refpx
        y1, y2 = refpx, refpx

        draw.line( [(x1,y1),(x2,y2)], fill=255, width=2)
        del draw
        outimg.save(out_fname)

    if rescaleimage:
        outimg = Image.open(out_fname)

        if isinstance(rescaleimage, tuple):
            size = rescaleimage
        else:
            size = (512, 512)

        outimg = outimg.resize(size, resample=Image.BILINEAR)

        outimg.save(out_fname)


    return out_fname

def fitscoords_to_jpeg(fits_image,
                       out_fname=None,
                       ext=None,
                       flip=True,
                       coordbox=None,
                       coordcenter=None,
                       annotatejd=True,
                       jdsrc=None,
                       scale_func=clipped_linscale_img,
                       scale_func_params={'cap':255.0,
                                          'lomult':2,
                                          'himult':2.5}):
    '''
    This converts a FITS image to a full frame JPEG.

    if coordbox and not coordcenter:
        coordbox = [xmin, xmax, ymin, max] of box to cut out of FITS

    elif coordcenter and not coordbox:
        coordcenter = [xcenter, ycenter, xwidth, ywidth]

    else:
        do nothing, since we can't have both at the same time


    '''
    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    trimmed_img = img
    jpegaspect = float(img.shape[1])/float(img.shape[0])
    scaled_img = scale_func(trimmed_img,**scale_func_params)

    if coordbox and not coordcenter:
        # numpy is y,x
        scaled_img = scaled_img[coordbox[2]:coordbox[3],
                                coordbox[0]:coordbox[1]]

    elif coordcenter and not coordbox:
        # numpy is y,x

        x1, x2 = (coordcenter[0] - coordcenter[2]/2.0,
                  coordcenter[0] + coordcenter[2]/2.0)
        y1, y2 = (coordcenter[1] - coordcenter[3]/2.0,
                  coordcenter[1] + coordcenter[3]/2.0)

        # figure out xmin, xmax, ymin, ymax
        if x1 > x2:
            xmin = x2
            xmax = x1
        else:
            xmin = x1
            xmax = x2

        if y1 > y2:
            ymin = y2
            ymax = y1
        else:
            ymin = y1
            ymax = y2

        # round the pix coords to integers
        xmin, xmax = int(np.round(xmin)), int(np.round(xmax))
        ymin, ymax = int(np.round(ymin)), int(np.round(ymax))

        # make sure we take care of edges
        if xmin < 0:
            xmin = 0
        if xmax >= img.shape[1]:
            xmax = img.shape[1] - 1
        if ymin < 0:
            ymin = 0
        if ymax >= img.shape[0]:
            ymax = img.shape[0] - 1

        scaled_img = scaled_img[ymin:ymax, xmin:xmax]

    if not out_fname:

        out_fname = '%s-%s-%s-%s-proj%s-%s.jpg' % (
            fits_image.rstrip('.fits.fz'),
            hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
            hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
            (hdr['FILTERS'].replace('+','') if
             'FILTERS' in hdr else 'filtunknown'),
            hdr['PROJID'] if 'PROJID' in hdr else 'unknown',
            hdr['OBJECT'] if 'OBJECT' in hdr else 'objectunknown'
            )

        if coordbox and not coordcenter:
            out_fname = '%s-X%sX%s-Y%sY%s.jpg' % (
                out_fname.rstrip('.jpg'),
                coordbox[0], coordbox[1],
                coordbox[2], coordbox[3]
            )

        elif coordcenter and not coordbox:
            out_fname = '%s-XC%sYC%s-XW%sYW%s.jpg' % (
                out_fname.rstrip('.jpg'),
                coordcenter[0], coordcenter[1],
                coordcenter[2], coordcenter[3]
            )

    # flip the saved image
    if flip:
        scaled_img = np.flipud(scaled_img)

    scipy.misc.imsave(out_fname, scaled_img)

    # annotate the image if told to do so
    if annotatejd and jdsrc and os.path.exists(jdsrc):

        # get the JD header keyword from jdsrc
        framejd = get_header_keyword(jdsrc, 'JD')

        outimg = Image.open(out_fname)
        draw = ImageDraw.Draw(outimg)
        annotation = "JD %.3f" % framejd
        draw.text((4,2),annotation,fill=255, font=fontsmall)

        del draw
        outimg.save(out_fname)

    return out_fname



def nparr_to_full_jpeg(nparr,
                       out_fname,
                       outsizex=770,
                       outsizey=770,
                       scale=True,
                       scale_func=clipped_linscale_img,
                       scale_func_params={'cap':255.0,
                                          'lomult':2,
                                          'himult':2.5}):
    '''
    This just writes a numpy array to a JPEG.

    '''
    if scale:
        scaled_img = scale_func(nparr,**scale_func_params)
    else:
        scaled_img = nparr

    resized_img = scipy.misc.imresize(scaled_img,
                                      (outsizex,outsizey))
    if out_fname is None:
        out_fname = fits_image + '.jpeg'
    scipy.misc.imsave(out_fname,resized_img)



def check_frame_warping(frame,
                        margins=50,
                        threshold=10.0,
                        showplot=False):
    '''This checks if an image is warped (perhaps by a bad shift/convolution).

    Calculates the median of the rows and columns of the image taking into
    account the margin on either side (as specified by the margins kwarg). Then
    fits a straight line to the trend. If the chi-sq of the fit is above the
    specified threshold, returns False as the image is likely to be
    warped. Otherwise, returns True.

    WARNING: the "threshold" for this to work depends very strongly on the
    image.  Particularly near the galactic plane, you should expect images to
    have some linear trend (i.e. brighter at lower galactic latitude).

    In such cases, a "high threshold" of ~20,000 might be appropriate. In other
    cases further from the galactic plane a "low threshold" of ~15,000 might be
    better. This is obviously heuristic empirical things, and a better method
    should be implemented.

    '''

    hdu = pyfits.open(frame)
    image = hdu[0].data
    hdu.close()

    clippedimage = image[margins:-margins, margins:-margins]
    imagecoordnum = np.arange(len(clippedimage))

    # get the medians in the x and y directions
    medx = np.nanmedian(clippedimage,axis=1)
    medy = np.nanmedian(clippedimage,axis=0)

    # fit a 1-degree polynomial
    xfitcoeffs = np.polyfit(imagecoordnum,medx,1)
    yfitcoeffs = np.polyfit(imagecoordnum,medy,1)

    xfitpoly = np.poly1d(xfitcoeffs)
    yfitpoly = np.poly1d(yfitcoeffs)

    xfit = xfitpoly(imagecoordnum)
    yfit = yfitpoly(imagecoordnum)

    xfit_redchisq = np.sum((medx - xfit)*(medx - xfit))/(len(imagecoordnum) - 2)
    yfit_redchisq = np.sum((medy - yfit)*(medy - yfit))/(len(imagecoordnum) - 2)

    warpinfo = {'medx':medx,
                'medy':medy,
                'xfitpoly':xfitpoly,
                'yfitpoly':yfitpoly,
                'xfit':xfit,
                'yfit':yfit,
                'xfit_redchisq':xfit_redchisq,
                'yfit_redchisq':yfit_redchisq}

    if showplot:

        import matplotlib.pyplot as plt
        from datetime import datetime
        plt.close('all')
        f, ax = plt.subplots()

        ax.plot(imagecoordnum, medx, 'k-')
        ax.plot(imagecoordnum, xfit, 'k-', alpha=0.5)
        ax.plot(imagecoordnum, medy, 'b-')
        ax.plot(imagecoordnum, yfit, 'b-', alpha=0.5)

        ax.set_xlabel('img coord number')
        ax.set_ylabel('x and y medians, and linear fits')

        savename = os.path.basename(frame).strip('.fits')+'_diagnostic.png'
        f.savefig(savename, bbox_inches='tight', dpi=250)

        print('%sZ: wrote diagnostic warp check plot to %s' %
              (datetime.utcnow().isoformat(), savename))


    if (xfit_redchisq > threshold) or (yfit_redchisq > threshold):
        return False, warpinfo
    else:
        return True, warpinfo


def make_mp4_from_jpegs(jpgglob, outmp4path, ffmpegpath='ffmpeg'):
    """
    Make mp4 movie from jpg images. (Codec/preset configured to work with
    ffmpeg v3.4. Fails for v4.X).

    Args:
        jpgglob: e.g.,
        /nfs/phtess1/ar1/TESS/FFI/RED_IMGSUB/FULL/s0001/RED_3-2-1011_ISP/JPEG-SUBTRACTEDCONV-rsub-9ab2774b-tess*cal_img-xtrns.jpg

        outmp4path: e.g.,
        /nfs/phtess1/ar1/TESS/FFI/MOVIES/s0001_full_cam3_ccd2_projid1011_SUBTRACTEDCONV.mp4
    """

    returncode = os.system('which ffmpeg')
    if not returncode==0:
        raise AssertionError(
            '`ffmpeg` must be installed to use make_mp4_from_jpegs')

    # framerate: obvious.
    # libx264: encoding
    # vf: encoding requires even number of pixels. This filter divided original
    # heigh and width by two, rounds up to nearest pixel, multiplies by two,
    # and add white padding pixels.
    FFMPEGCMD = (
        "{ffmpegpath} -framerate 24 "
        "-pattern_type "
        "glob -i '{jpgglob}' "
        "-c:v libx264 "
        "-preset fast "
        "-vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white\" "
        "{outmp4path}"
    )

    cmdtorun = FFMPEGCMD.format(jpgglob=jpgglob,
                                outmp4path=outmp4path,
                                ffmpegpath=ffmpegpath)

    returncode = os.system(cmdtorun)

    if returncode == 0:
        print('%sZ: made movie %s' %
              (datetime.utcnow().isoformat(), outmp4path))
        return 0
    else:
        print('ERR! %sZ: failed to make movie %s' %
              (datetime.utcnow().isoformat(), outmp4path))
        print('ERR! command was %s' % cmdtorun)
        return 256


def make_mov_from_jpegs(jpgglob, outmovpath, ffmpegpath='ffmpeg'):
    """
    Similar to above, but makes .mov (a format that is compatible with e.g.,
    keynote) (Codec/preset configured to work with ffmpeg v3.4. Fails for
    v4.X).


    Args:
        jpgglob: e.g.,
        /nfs/phtess1/ar1/TESS/FFI/RED_IMGSUB/FULL/s0001/RED_3-2-1011_ISP/JPEG-SUBTRACTEDCONV-rsub-9ab2774b-tess*cal_img-xtrns.jpg

        outmovpath: e.g.,
        /nfs/phtess1/ar1/TESS/FFI/MOVIES/s0001_full_cam3_ccd2_projid1011_SUBTRACTEDCONV.mov
    """

    returncode = os.system('which ffmpeg')
    if not returncode==0:
        raise AssertionError(
            '`ffmpeg` must be installed to use make_mov_from_jpegs')

    FFMPEGCMD = (
        "{ffmpegpath} -framerate 24 "
        "-pattern_type "
        "glob -i '{jpgglob}' "
        "-c:v libx264 "
        "-pix_fmt yuv420p "
        "-preset fast "
        "{outmovpath}"
    )

    cmdtorun = FFMPEGCMD.format(
        ffmpegpath=ffmpegpath,
        jpgglob=jpgglob,
        outmovpath=outmovpath
    )

    returncode = os.system(cmdtorun)

    if returncode == 0:
        print('%sZ: made movie %s' %
              (datetime.utcnow().isoformat(), outmovpath))
        return 0
    else:
        print('ERR! %sZ: failed to make movie %s' %
              (datetime.utcnow().isoformat(), outmovpath))
        print('ERR! command was %s' % cmdtorun)
        return 256


def plot_stages_of_img_proc_sector_cam_ccd(
    sector=6, cam=1, ccd=2, projid=1501, overwrite=0, outdir=None,
    slicebounds=[slice(300,812), slice(300,812)]
    ):

    if not isinstance(outdir, str):
        raise NotImplementedError
    if not os.path.exists(outdir):
        raise ValueError('did not find {}'.format(outdir))

    datadir = (
        '/nfs/phtess2/ar0/TESS/FFI/RED/sector-{}/cam{}_ccd{}/'.
        format(sector, cam, ccd)
    )
    diffdir = (
        '/nfs/phtess2/ar0/TESS/FFI/RED_IMGSUB/FULL/s{}/RED_{}-{}-{}_ISP'.
        format(str(sector).zfill(4), cam, ccd, projid)
    )

    bkgdfiles = glob(os.path.join(
        datadir,
        'tess20*-s{}-{}-{}-*_cal_img_bkgd.fits'.
        format(str(sector).zfill(4), cam, ccd))
    )
    calfiles = glob(os.path.join(
        datadir,
        'tess20*-s{}-{}-{}-*_cal_img.fits'.
        format(str(sector).zfill(4), cam, ccd))
    )
    difffiles = glob(os.path.join(
        diffdir,
        'rsub-*-*-s{}-{}-{}-*_cal_img_bkgdsub-xtrns.fits'.
        format(str(sector).zfill(4), cam, ccd))
    )

    for b,c,d in zip(bkgdfiles, calfiles, difffiles):

        outpath = os.path.join(
            outdir,
            os.path.basename(c).replace('_cal_img.fits',
                                        '_img_proc_stages.png')
        )

        if os.path.exists(outpath) and not overwrite:
            print('found {} and not overwrite; continue'.format(outpath))
            continue
        else:
            plot_stages_of_img_proc(b, c, d, outpath, slicebounds=slicebounds)


def plot_stages_of_img_proc(
    bkgdfile, calfile, difffile, outpath,
    slicebounds=[slice(300,812), slice(300,812)]
    ):

    vmin, vmax = 10, int(1e3)

    bkgd_img, _ = read_fits(bkgdfile)
    cal_img, _ = read_fits(calfile)
    diff_img, _ = read_fits(difffile)

    plt.close('all')
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    fig, axs = plt.subplots(ncols=2, nrows=3)

    # top right: log of calibrated image
    lognorm = colors.LogNorm(vmin=vmin, vmax=vmax)

    cset1 = axs[0,1].imshow(cal_img, cmap='binary_r', vmin=vmin, vmax=vmax,
                            norm=lognorm)
    #txt = axs[0,1].text(0.02, 0.96, 'image', ha='left', va='top',
    #                    fontsize='small', transform=axs[0,1].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])


    diff_vmin, diff_vmax = -1000, 1000

    diffnorm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=diff_vmin,
                                 vmax=diff_vmax)

    # top left: background map
    axs[0,0].imshow(bkgd_img - np.median(cal_img), cmap='RdBu_r',
                    vmin=diff_vmin, vmax=diff_vmax, norm=diffnorm)
    #txt = axs[0,0].text(0.02, 0.96, 'background', ha='left', va='top',
    #                    fontsize='small', transform=axs[0,0].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])

    # middle left: calibrated - background
    cset2 = axs[1,0].imshow(cal_img - bkgd_img, cmap='RdBu_r', vmin=diff_vmin,
                            vmax=diff_vmax, norm=diffnorm)
    #txt = axs[1,0].text(0.02, 0.96, 'image - background', ha='left', va='top',
    #                    fontsize='small', transform=axs[1,0].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])

    # middle right: calibrated - median
    axs[1,1].imshow(cal_img - np.median(cal_img), cmap='RdBu_r',
                    vmin=diff_vmin, vmax=diff_vmax, norm=diffnorm)
    #txt = axs[1,1].text(0.02, 0.96, 'image - median(image)', ha='left', va='top',
    #                    fontsize='small', transform=axs[1,1].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])


    # lower left:  difference image (full)
    toplen = 57
    top = mplcm.get_cmap('Oranges_r', toplen)
    bottom = mplcm.get_cmap('Blues', toplen)
    newcolors = np.vstack((top(np.linspace(0, 1, toplen)),
                           np.zeros(((256-2*toplen),4)),
                           bottom(np.linspace(0, 1, toplen))))
    newcmp = ListedColormap(newcolors, name='lgb_cmap')

    cset3 = axs[2,0].imshow(diff_img, cmap='RdBu_r', vmin=diff_vmin,
                            vmax=diff_vmax, norm=diffnorm)
    #txt = axs[2,0].text(0.02, 0.96, 'difference', ha='left', va='top',
    #                    fontsize='small', transform=axs[2,0].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])

    xmin, xmax, ymin, ymax = (slicebounds[0].start, slicebounds[0].stop,
                              slicebounds[1].start, slicebounds[1].stop)
    width = ymax - ymin
    height = xmax - xmin

    rect = patches.Rectangle((ymin, xmin), width, height, linewidth=0.6,
                             edgecolor='black', facecolor='none',
                             linestyle='--')
    axs[2,0].add_patch(rect)



    # lower right: difference image (zoom)
    sel = slicebounds
    axs[2,1].imshow(diff_img[sel], cmap='RdBu_r',
                    vmin=diff_vmin, vmax=diff_vmax, norm=diffnorm)
    #txt = axs[2,1].text(0.02, 0.96, 'difference (zoom)', ha='left', va='top',
    #                    fontsize='small', transform=axs[2,1].transAxes,
    #                    color='black')
    #txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
    #                      path_effects.Normal()])

    for ax in axs.flatten():
        ax.set_xticklabels('')
        ax.set_yticklabels('')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')

    divider0 = make_axes_locatable(axs[0,1])
    divider1 = make_axes_locatable(axs[1,1])
    divider2 = make_axes_locatable(axs[2,1])

    cax0 = divider0.append_axes('right', size='5%', pad=0.05)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)

    cb1 = fig.colorbar(cset1, ax=axs[0,1], cax=cax0, extend='both')
    cb2 = fig.colorbar(cset2, ax=axs[1,1], cax=cax1, extend='both')
    cb3 = fig.colorbar(cset3, ax=axs[2,1], cax=cax2, extend='both')

    cb2.set_ticks([-1e3,-1e2,-1e1,0,1e1,1e2,1e3])
    cb2.set_ticklabels(['-$10^3$','-$10^2$','-$10^1$','0',
                        '$10^1$','$10^2$','$10^3$'])
    cb3.set_ticks([-1e3,-1e2,-1e1,0,1e1,1e2,1e3])
    cb3.set_ticklabels(['-$10^3$','-$10^2$','-$10^1$','0',
                        '$10^1$','$10^2$','$10^3$'])

    fig.tight_layout(h_pad=0, w_pad=-14, pad=0)

    fig.savefig(outpath, bbox_inches='tight', dpi=400)
    print('{}: made {}'.format(datetime.utcnow().isoformat(), outpath))

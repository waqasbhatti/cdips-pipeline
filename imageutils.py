#!/usr/bin/env python

'''
imageutils.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Jan 2013

This contains various utilities suitable for quick-processing HAT
images. Intended for use with the real-time image diagnostics functionality of
webcontrol interface.

Contains functions for:

1. reading quickastrom files, and getting out the tracking correction, number of
   sources, etc.
2. reading hatphot output and generating FWHM statistics for a grid placed on
   the image
3. generating stamps for an image, converting an image to JPEGs
4. getting the value of a certain keyword from the FITS header for a series of
   FITS files

Important FITS header keywords:

FOCUS (steps)
BJD
MIDTIME (HH:MM:SS.SSS - middle of exposure)
MIDDATE (YYYY-MM-DD - middle of exposure)
TIMESYS (should be UTC)
OBJECT (field name)
JD (JD of middle exposure)
HA (hour angle)
Z (zenith distance)
ABORTED (0 = exp not aborted, 1 = aborted exp)
IMAGETYP

'''

import os
import os.path
import sys
import logging
import bz2
import base64
import glob

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

import pyfits

from PIL import Image
from PIL import ImageDraw

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

    image_flist = sorted(glob.glob(image_glob))

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
    '''
    Get the value of a header keyword in a fits file optionally using an
    extension.

    '''
    hdulist = pyfits.open(fits_file)

    if keyword in hdulist[ext].header:
        val = hdulist[ext].header[keyword]
    else:
        val = None

    hdulist.close()
    return val


def get_header_keyword_list(fits_file,
                            keyword_list,
                            ext=0):

    hdulist = pyfits.open(fits_file)

    out_dict = {}

    for keyword in keyword_list:

        if keyword in hdulist[ext].header:
            out_dict[keyword] = hdulist[ext].header[keyword]
        else:
            out_dict[keyword] = None

    hdulist.close()
    return out_dict


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
        backmin = image_min if image_min is not None else np.min(img_array)

    else:

        backmin, backmax = custom_limits

    masked = npma.masked_outside(img_array, backmin, backmax)
    backmasked = npma.median(masked)

    return backmasked


## IMAGE SECTION FUNCTIONS ##

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
                      outsizex=770,
                      outsizey=770,
                      scale_func=clipped_linscale_img,
                      scale_func_params={'cap':255.0,
                                         'lomult':2,
                                         'himult':2.5}):
    '''
    This converts a FITS image to a full frame JPEG.

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
    
    resized_img = scipy.misc.imresize(scaled_img,
                                      (int(img.shape[1]/2.2),int(img.shape[0]/2.2)))
    out_fname = '%s-%s-%s-%s-%s.jpg' % (fits_image.rstrip('.fits.fz'),
                                        hdr['IMAGETYP'].lower() if 'IMAGETYP' in hdr else 'typeunknown',
                                        hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown',
                                        hdr['FILTERS'].replace('+','') if 'FILTERS' in hdr else 'filtunknown',
                                        hdr['PROJID'] if 'PROJID' in hdr else 'projunknown')
    
    print('fits: %s, xsize = %s, ysize = %s, x/y = %s, out = %s' % (fits_image, 
                                                                    img.shape[1], 
                                                                    img.shape[0], 
                                                                    jpegaspect,
                                                                    out_fname))
    scipy.misc.imsave(out_fname,resized_img)
    # flip the saved image
    outimg = Image.open(out_fname)
    outimg = outimg.transpose(Image.FLIP_TOP_BOTTOM)
    outimg.save(out_fname)


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


def stamp_flatten(object_img,
                  flat_img):
    '''
    This flattens the image stamp-by-stamp.

    '''

    # stamp both images
    stamped_object = img_to_stamps(object_img)
    stamped_flat = img_to_stamps(flat_img)

    stamped_output = {}

    for k in stamped_object:

        flattened_stamp = stamped_object[k]/np.median(stamped_flat[k])
        stamped_output[k] = flattened_stamp

    return stamped_output


def image_flatten(object_img,
                  flat_img,
                  ext=None,
                  normflat=False):

    compressed_ext = compressed_fits_ext(object_img)

    if ext is None and compressed_ext:
        objimg, objhdr = read_fits(object_img,
                                   ext=compressed_ext[0])
    elif (ext is not None):
        objimg, objhdr = read_fits(object_img,ext=ext)
    else:
        objimg, objhdr = read_fits(object_img)

    trimmed_object_img = trim_image(objimg, objhdr)

    compressed_ext = compressed_fits_ext(flat_img)

    if ext is None and compressed_ext:
        flatimg, flathdr = read_fits(flat_img,
                                     ext=compressed_ext[0])
    elif (ext is not None):
        flatimg, flathdr = read_fits(flat_img,ext=ext)
    else:
        flatimg, flathdr = read_fits(flat_img)

    if normflat:
        norm_superflat = flatimg/np.median(flatimg)
        flattened_image = trimmed_object_img/norm_superflat
    else:
        flattened_image = trimmed_object_img/flatimg


    return flattened_image



def flattened_stamps_to_jpeg(flattened_stamp,
                             out_fname,
                             altstamping=False,
                             scale=True,
                             scale_func=clipped_linscale_img,
                             scale_func_params={'cap':255.0,
                                                'lomult':2,
                                                'himult':2.5}):
    '''
    This uses PIL to pull stamped images together into a single image.

    Assumes number of stamps is 3 x 3 and size is 256 x 256.

    '''

    # this pulls together stamps into one 768 x 768 image, scales it, stamps it
    # again, puts in the 3 x 3 grid lines, and finally writes the file out to
    # JPEG
    if altstamping:

        full_frame = np.zeros((768,768))
        xstampsize, ystampsize = 256, 256

        # make up the full_frame to be scaled from the stamp regions
        full_frame[:xstampsize,:ystampsize] = flattened_stamp['topleft']
        full_frame[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                   :ystampsize] = flattened_stamp['topcenter']
        full_frame[imgsizex-xstampsize:,:ystampsize] = (
            flattened_stamps['topright']
            )

        full_frame[:xstampsize,
                    imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2] = flattened_stamps['midleft']
        full_frame[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                   imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2] = flattened_stamps['midcenter']
        full_frame[imgsizex-xstampsize:,
                   imgsizey/2-ystampsize/2:imgsizey/2+ystampsize/2] = flattened_stamp['midright']

        full_frame[:xstampsize,imgsizey-ystampsize:] = flattened_stamp['bottomleft']
        full_frame[imgsizex/2-xstampsize/2:imgsizex/2+xstampsize/2,
                   imgsizey-ystampsize:] = flattened_stamp['bottomcenter']
        full_frame[-xstampsize:,-ystampsize:] = flattened_stamp['bottomright']

        # now scale this full image
        scaled_frame = scale_func(full_frame,**scale_func_params)

        # restamp the scaled_frame
        restamped_frame = img_to_stamps(scaled_frame)

        # write the stamped image out to JPEG
        stamps_to_jpeg(restamped_frame,
                       out_fname,
                       scale=False)

    else:

        full_image = Image.new('L',(770,770))

        subframe_coords = {'topleft':(0,0),
                           'topcenter':(0,256),
                           'topright':(0,512),
                           'bottomleft':(512,0,),
                           'bottomcenter':(512,256),
                           'bottomright':(512,512),
                           'midleft':(256,0),
                           'midcenter':(256,256),
                           'midright':(256,512),
                           }


        # paste in the subframes
        for subframe in flattened_stamp:

            stamp = flattened_stamp[subframe]
            if scale:
                scaled_stamp = scale_func(stamp,**scale_func_params)
            else:
                scaled_stamp = stamp

            stamp_img = scipy.misc.toimage(scaled_stamp)
            full_image.paste(stamp_img,subframe_coords[subframe])

        # put the lines to distinguish between stamps
        full_image_draw = ImageDraw.Draw(full_image)
        full_image_draw.line([(255,0),(255,770)],fill=255)
        full_image_draw.line([(511,0),(511,770)],fill=255)
        full_image_draw.line([(0,255),(770,255)],fill=255)
        full_image_draw.line([(0,511),(770,511)],fill=255)

        full_image.save(out_fname)


def flattened_fits_to_stamps_jpeg(fits_image,
                                  smoothed_flat_image,
                                  out_fname=None,
                                  ext=None,
                                  stampsize=256,
                                  sepwidth=1,
                                  scale_func=clipped_linscale_img,
                                  scale_func_params={'cap':255.0,
                                                     'lomult':2,
                                                     'himult':2.5}):
    '''
    This applies an existing smoothed_flat_image to the fits_image as a flat
    field, then returns the image as stamps.

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

    # read and apply the flat image
    flat_img, flat_hdr = read_fits(smoothed_flat_image)

    stamped_flattened_img = stamp_flatten(trimmed_img,
                                          flat_img)

    # generate the stamp JPEG
    if out_fname is None:
        out_fname = fits_image + '.stamp.jpeg'
    flattened_stamps_to_jpeg(stamped_flattened_img,
                             out_fname,
                             scale=True,
                             scale_func=scale_func,
                             scale_func_params=scale_func_params)


def flattened_fits_to_jpeg(fits_image,
                           flat_image,
                           out_fname=None,
                           ext=None,
                           outsizex=770,
                           outsizey=770,
                           stampsize=256,
                           sepwidth=1,
                           scale_func=clipped_linscale_img,
                           scale_func_params={'cap':255.0,
                                              'lomult':2,
                                              'himult':2.5}):
    '''
    This is a version of the function below that uses flattened images.

    '''

    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    try:
        image_type = hdr['IMAGETYP']
    except KeyError:
        LOGGER.warning('no imagetype in FITS header, not making stamps')
        image_type = 'unknown'

    # get the naxis of the frame to see if this is a multi-dim FITS
    frame_naxis = hdr['NAXIS'] if 'NAXIS' in hdr else 1

    # if there's more than one image plane, then collapse all of them into
    # one
    if frame_naxis > 2:
        img = img[0]

    trimmed_img = trim_image(img, hdr)

    if 'object' in image_type or 'focus' in image_type:

        # read and apply the flat image
        flat_img, flat_hdr = read_fits(flat_image)
        stamped_flattened_img = stamp_flatten(trimmed_img,
                                              flat_img)

        # generate the stamp JPEG
        if out_fname is None:
            out_fname = fits_image + '.stamp.jpeg'

        flattened_stamps_to_jpeg(stamped_flattened_img,
                                 out_fname,
                                 scale=True,
                                 scale_func=scale_func,
                                 scale_func_params=scale_func_params)
    else:

        scaled_img = scale_func(trimmed_img,**scale_func_params)
        resized_img = scipy.misc.imresize(scaled_img,
                                          (outsizex,outsizey))
        if out_fname is None:
            out_fname = fits_image + '.jpeg'
        scipy.misc.imsave(out_fname,resized_img)

    return out_fname, trimmed_img.shape



def alt_flattened_fits_to_jpeg(fits_image,
                               flat_image,
                               out_fname=None,
                               ext=None,
                               outsizex=770,
                               outsizey=770,
                               stampsize=256,
                               sepwidth=1,
                               normflat=True,
                               scale_func=clipped_linscale_img,
                               scale_func_params={'cap':255.0,
                                                  'lomult':2,
                                                  'himult':2.5}):
    '''
    This is a version of the function below that uses flattened images and the
    alternative flatting method.

    '''

    flattened_image = image_flatten(fits_image,
                                    flat_image,
                                    normflat=normflat)


    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    # get the naxis of the frame to see if this is a multi-dim FITS
    frame_naxis = hdr['NAXIS'] if 'NAXIS' in hdr else 1

    # if there's more than one image plane, then collapse all of them into
    # one
    if frame_naxis > 2:
        img = img[0]

    try:
        image_type = hdr['IMAGETYP']
    except KeyError:
        LOGGER.warning('no imagetype in FITS header, not making stamps')
        image_type = 'unknown'

    scaled_flattened_image = scale_func(flattened_image,
                                        **scale_func_params)

    if 'object' in image_type or 'focus' in image_type:

        # generate the stamp JPEG
        if out_fname is None:
            out_fname = fits_image + '.stamp.jpeg'


        flattened_stamps = img_to_stamps(scaled_flattened_image,
                                         stampsize=stampsize)
        stamps_to_jpeg(flattened_stamps, out_fname)
        return out_fname, scaled_flattened_image.shape

    else:

        resized_img = scipy.misc.imresize(scaled_flattened_image,
                                          (outsizex,outsizey))
        if out_fname is None:
            out_fname = fits_image + '.jpeg'
        scipy.misc.imsave(out_fname,resized_img)

    return out_fname, trimmed_img.shape



def alt2_flattened_fits_to_jpeg(fits_image,
                                flat_image,
                                out_fname=None,
                                ext=None,
                                outsizex=770,
                                outsizey=770,
                                stampsize=256,
                                sepwidth=1,
                                normflat=True,
                                scale_func=clipped_linscale_img,
                                scale_func_params={'cap':255.0,
                                                   'lomult':5.0,
                                                   'himult':2.0}):
    '''
    This is another version of the function above that uses flattened images and
    the alternative flattening method. This one does flattening per stamp, but
    then puts all stamps into one image, scales that, restamps it and then
    writes to JPEG (yes, it's convoluted).

    '''

    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    try:
        image_type = hdr['IMAGETYP']
    except KeyError:
        LOGGER.warning('no imagetype in FITS header, not making stamps')
        image_type = 'unknown'

    # get the naxis of the frame to see if this is a multi-dim FITS
    frame_naxis = hdr['NAXIS'] if 'NAXIS' in hdr else 1

    # if there's more than one image plane, then collapse all of them into
    # one
    if frame_naxis > 2:
        img = img[0]

    trimmed_img = trim_image(img, hdr)

    if 'object' in image_type or 'focus' in image_type:

        # read and apply the flat image
        flat_img, flat_hdr = read_fits(flat_image)
        stamped_flattened_img = stamp_flatten(trimmed_img,
                                              flat_img)
        # generate the stamp JPEG
        if out_fname is None:
            out_fname = fits_image + '.stamp.jpeg'

        flattened_stamp_to_jpeg(stamped_flatttened_img,
                                out_fname,
                                altstamping=True,
                                scale=False)

        return out_fname, flattened_image.shape

    else:

        scaled_flattened_image = scale_func(trimmed_img,
                                            **scale_func_params)
        resized_img = scipy.misc.imresize(scaled_flattened_image,
                                          (outsizex,outsizey))
        if out_fname is None:
            out_fname = fits_image + '.jpeg'
        scipy.misc.imsave(out_fname,resized_img)

    return out_fname, trimmed_img.shape



def fits_to_jpeg(fits_image,
                 out_fname=None,
                 ext=None,
                 outsizex=770,
                 outsizey=770,
                 stampsize=256,
                 sepwidth=1,
                 scale_func=clipped_linscale_img,
                 scale_func_params={'cap':255.0,
                                    'lomult':2,
                                    'himult':2.5}):
    '''
    This figures out what kind of FITS image it is reading, and turns it into a
    stamp or full-frame JPEG based on the IMAGETYP FITS header keyword.

    '''

    compressed_ext = compressed_fits_ext(fits_image)

    if ext is None and compressed_ext:
        img, hdr = read_fits(fits_image,
                             ext=compressed_ext[0])
    elif (ext is not None):
        img, hdr = read_fits(fits_image,ext=ext)
    else:
        img, hdr = read_fits(fits_image)

    # get the naxis of the frame to see if this is a multi-dim FITS
    frame_naxis = hdr['NAXIS'] if 'NAXIS' in hdr else 1

    # if there's more than one image plane, then collapse all of them into
    # one
    if frame_naxis > 2:
        img = img[0]

    trimmed_img = trim_image(img, hdr)
    scaled_img = scale_func(trimmed_img,**scale_func_params)

    try:
        image_type = hdr['IMAGETYP']
    except KeyError:
        LOGGER.warning('no imagetype in FITS header, not making stamps')
        image_type = 'unknown'


    if 'object' in image_type or 'focus' in image_type:
        stamps = img_to_stamps(scaled_img,stampsize=stampsize)
        if out_fname is None:
            out_fname = fits_image + '.stamp.jpeg'
        out_img = stamps_to_jpeg(stamps,out_fname,sepwidth=sepwidth)
    else:
        resized_img = scipy.misc.imresize(scaled_img,
                                          (outsizex,outsizey))
        if out_fname is None:
            out_fname = fits_image + '.jpeg'
        scipy.misc.imsave(out_fname,resized_img)

    return out_fname, trimmed_img.shape


def jpeg_to_zb64(image_fname):
    '''
    This is used to convert a JPEG snapshot image to a compressed base-64 text
    representation suitable for transfer via JSON over a ZMQ wire and for
    storage in an SQite DB.

    '''

    img = open(image_fname,'rb')
    jpeg_str = img.read()
    img.close()

    b64_jpeg = base64.b64encode(jpeg_str)

    if LOGGER:
        LOGGER.info('uncompressed base64 image size: %s' % len(b64_jpeg))
    else:
        print('imagedb.jpeg_to_zb64: uncompressed base64 image size: %s' %
              len(b64_jpeg))

    bz2_string = bz2.compress(jpeg_str)

    zipped_img_str = base64.b64encode(bz2_string)

    if LOGGER:
        LOGGER.info('compressed base64 image size: %s' % len(zipped_img_str))
    else:
        print('imagedb.jpeg_to_zb64: compressed base64 image size: %s' %
              len(zipped_img_str))

    return zipped_img_str


def zb64_to_jpeg(z_string, out_fname):
    '''
    This converts a compressed base-64 text representation of a JPEG snapshot
    back into a JPEG. Writes out to the file with path out_fname.

    '''

    b64_str = bz2.decompress(base64.b64decode(z_string))
    outf = open(out_fname,'wb')
    outf.write(b64_str)
    outf.close()


def jpeg_to_b64(image_fname):
    '''
    This is used to convert a JPEG snapshot image to a base-64 text
    representation suitable for transfer via JSON over a ZMQ wire and for
    storage in an SQLite DB.

    '''

    img = open(image_fname,'rb')
    b64_str = base64.b64encode(img.read())
    img.close()
    return b64_str


def b64_to_jpeg(b64_str, out_fname):
    '''
    This converts a base-64 text representation of a JPEG snapshot back into a
    JPEG. Writes out to the file with path out_fname.

    '''

    jpeg_str = base64.b64decode(b64_str)
    outf = open(out_fname,'wb')
    outf.write(jpeg_str)
    outf.close()


## DIAGNOSTIC UTILITY FUNCTIONS ##

def read_hatlog(hatlog_path):
    '''
    This reads a HAT.LOG file, and splits it by line. Use it to read the HAT.LOG
    file for a specific night's image directory.

    '''

    hatlog_path = os.path.abspath(hatlog_path)
    hatlog = open(hatlog_path,'rb')
    hatlog_lines = hatlog.readlines()
    hatlog.close()
    hatlog_lines = [line.strip() for line in hatlog_lines]
    return hatlog_lines[1:]


def analyze_qastrom(hatlog_lines, framenumber=None, singleccd=False):
    '''
    This gets the qastrom corrections from the HAT.LOG.

    '''

    qastrom_runs = [x.split() for x in hatlog_lines if 'qastrom::run' in x]

    # handle the usual case of guiding on 4-CCD overlap
    if singleccd is False:

        # if a framenumber is provided, use it to cut down on the lines we have
        # to process
        if framenumber:
            # get correction values only if corrections were successful
            qastrom_corrections = [x.split() for x in hatlog_lines
                                   if (('Corr:' in x) and
                                       ('wrn' not in x) and
                                       ('err' not in x) and
                                       (framenumber in x))]

        else:
            # get correction values only if corrections were successful
            qastrom_corrections = [x.split() for x in hatlog_lines
                                   if (('Corr:' in x) and
                                       ('wrn' not in x) and
                                       ('err' not in x))]

    # handle the special case of single CCD guiding
    elif singleccd is True:

        # if a framenumber is provided, use it to cut down on the lines we have
        # to process
        if framenumber:
            # get correction values only if corrections were successful
            qastrom_corrections = [x.split() for x in hatlog_lines
                                   if (('Corr:' in x) and ('guiding' in x) and
                                       ('wrn' not in x) and
                                       ('err' not in x) and
                                       (framenumber in x))]

        else:
            # get correction values only if corrections were successful
            qastrom_corrections = [x.split() for x in hatlog_lines
                                   if (('Corr:' in x) and ('guiding' in x) and
                                       ('wrn' not in x) and
                                       ('err' not in x))]

    n_runs = len(qastrom_runs)
    n_corrs = len(qastrom_corrections)

    # split and strip the corrections
    if n_corrs > 0:

        qastrom_time = [' '.join([x[0],x[1]]) for x in qastrom_corrections]
        qastrom_corr = [x[-1] for x in qastrom_corrections]
        qastrom_corr = [x.strip('()') for x in qastrom_corr]
        qastrom_corr = [x.split(',') for x in qastrom_corr]

        # transpose and get out the x and y direction corrections
        qastrom_corr = zip(*qastrom_corr)
        try:

            qastrom_xcorr = [float(x) for x in qastrom_corr[0]]
            qastrom_ycorr = [float(x) for x in qastrom_corr[1]]

        except ValueError:

            if LOGGER:
                LOGGER.error('qastrom correction parsing failed')
                LOGGER.error('offending qastrom_corr was = %s' % repr(qastrom_corr))
            else:
                print('qastrom correction parsing failed')
                print('offending qastrom_corr was = %s' % repr(qastrom_corr))

            qastrom_xcorr = [0.0]
            qastrom_ycorr = [0.0]

    else:

        qastrom_time = []
        qastrom_xcorr = []
        qastrom_ycorr = []

    return {'n_run':n_runs,
            'n_corr':n_corrs,
            'qastrom_time':qastrom_time,
            'qastrom_xcorr':qastrom_xcorr,
            'qastrom_ycorr':qastrom_ycorr}


def analyze_focus(hatlog_lines):
    '''
    This gets the focus corrections from the HAT.LOG.

    '''

    focus_corrections = [x.split()[-1] for x in hatlog_lines if
                         'Focus correction for ' in x]

    if len(focus_corrections) > 0:
        focus_corrections = [float(x) for x in focus_corrections]
        n_focus_corr = len(focus_corrections)

    else:
        focus_corrections = []
        n_focus_corr = 0

    return {'n_corr':n_focus_corr,
            'focus_corr':focus_corrections}


def read_hatphot(hatphot_file):
    '''
    Reads a hatphot file, returns columns as arrays.

    '''

    hatphot = open(hatphot_file,'r')
    hatphot_lines = hatphot.readlines()
    hatphot.close()

    hatphot_lines = [x.strip('\n') for x in hatphot_lines if x[0] != '#']
    hatphot_lines = [x.split() for x in hatphot_lines]
    hatphot_cols = zip(*hatphot_lines)

    if len(hatphot_cols) >= 5:

        return {'framenum':os.path.basename(hatphot_file).strip('.hatphot'),
                'x':hatphot_cols[1],
                'y':hatphot_cols[2],
                'fwhm':hatphot_cols[5],
                'ellip':hatphot_cols[7],
                'ndet':len(hatphot_lines)}
    else:

        return None


def process_hatphot(hatphot_coldict,
                    imgsizex=4096,
                    imgsizey=4096,
                    stampsize=256):
    '''
    This takes the result of read_hatphot above and figures out the FWHM and
    ellipticity values in a 3x3 stamp of the image.

    Returned FWHM and ellipticity stamp arrays are in order:

    [topleft, topcenter, topright,
     midleft, midcenter, midright,
     botleft, botcenter, botright]

    '''

    xstampsize, ystampsize = stampsize, stampsize

    x, y = hatphot_coldict['x'], hatphot_coldict['y']
    fwhm, ellip = hatphot_coldict['fwhm'], hatphot_coldict['ellip']

    x = np.array([float(a) for a in x])
    y = np.array([float(a) for a in y])
    fwhm = np.array([float(a) for a in fwhm])
    ellip = np.array([float(a) for a in ellip])

    topleft_i = np.where( (x < xstampsize) &
                          (y < ystampsize) )
    topcenter_i = np.where( (x > (imgsizex/2-xstampsize/2)) &
                            (x < (imgsizex/2+xstampsize/2)) &
                            (y < ystampsize) )
    topright_i = np.where( (x > (imgsizex-xstampsize)) &
                           (y < ystampsize) )
    midleft_i = np.where( (x < xstampsize) &
                          (y > (imgsizey/2-ystampsize/2)) &
                          (y < (imgsizey/2+ystampsize/2)) )
    midcenter_i = np.where( (x > (imgsizex/2-xstampsize/2)) &
                            (x < (imgsizex/2+xstampsize/2)) &
                            (y > (imgsizey/2-ystampsize/2)) &
                            (y < (imgsizey/2+ystampsize/2)) )
    midright_i = np.where( (x > (imgsizex-xstampsize)) &
                           (y > (imgsizey/2-ystampsize/2)) &
                           (y < (imgsizey/2+ystampsize/2)) )
    bottomleft_i = np.where( (x < xstampsize) &
                             (y > (imgsizey - ystampsize)) )
    bottomcenter_i = np.where( (x > (imgsizex/2-xstampsize/2)) &
                               (x < (imgsizex/2+xstampsize/2)) &
                               (y > (imgsizey-ystampsize)) )
    bottomright_i = np.where( (x > (imgsizex-xstampsize)) &
                              (y > (imgsizey-ystampsize)) )

    fwhm_median = [
        np.median(x) for x in [fwhm[topleft_i],
                               fwhm[topcenter_i],
                               fwhm[topright_i],
                               fwhm[midleft_i],
                               fwhm[midcenter_i],
                               fwhm[midright_i],
                               fwhm[bottomleft_i],
                               fwhm[bottomcenter_i],
                               fwhm[bottomright_i]]
        ]

    ellip_median = [
        np.median(x) for x in [ellip[topleft_i],
                               ellip[topcenter_i],
                               ellip[topright_i],
                               ellip[midleft_i],
                               ellip[midcenter_i],
                               ellip[midright_i],
                               ellip[bottomleft_i],
                               ellip[bottomcenter_i],
                               ellip[bottomright_i]]
        ]



    return {
        'framenum':hatphot_coldict['framenum'],
        'fwhm_stamp':fwhm_median,
        'ellip_stamp':ellip_median,
        'ndet':hatphot_coldict['ndet'],
        }


def analyze_hatphot(hatphot_file,
                    imgsizex=4096,
                    imgsizey=4096,
                    stampsize=256):
    '''
    This puts together the two functions above.

    Returned FWHM and ellipticity stamp arrays are in order:

    [topleft, topcenter, topright,
     midleft, midcenter, midright,
     botleft, botcenter, botright]

    '''

    coldict = read_hatphot(hatphot_file)
    if coldict:
        hatphot_stamp = process_hatphot(coldict,
                                        imgsizex=imgsizex,
                                        imgsizey=imgsizey,
                                        stampsize=stampsize)

        return hatphot_stamp
    else:
        return None


def process_strans(strans_file,arcsec=True):
    '''
    This reads the strans file and returns the values of the following fields:

    Pointing offset relative to header (dRA, dDEC)
    Field center offset relative to header (dRA, dDEC)

    The first of these is taken as the nominal quickastrom generated correction
    to the pointing of the telescope.

    '''

    strans = open(strans_file,'r')
    strans_contents = strans.readlines()
    strans.close()

    p_offset = [x.strip('#\n').split(':')[-1].strip() for x in strans_contents if
                'Pointing offset' in x][0]
    fc_offset = [x.strip('#\n').split(':')[-1].strip() for x in strans_contents if
                 'Field center offset' in x][0]

    p_offset = p_offset.split(',')
    fc_offset = fc_offset.split(',')


    if arcsec:
        # get the offsets in arcseconds
        p_offset = [float(x)*3600.0 for x in p_offset]
        fc_offset = [float(x)*3600.0 for x in fc_offset]
    else:
        # get the offsets (they're already in arcsec)
        p_offset = [float(x) for x in p_offset]
        fc_offset = [float(x) for x in fc_offset]

    return {'p_offset':p_offset,
            'fc_offset':fc_offset}

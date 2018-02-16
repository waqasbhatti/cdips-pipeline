#!/usr/bin/env python

'''
fits2jpeg.py - Waqas Bhatti (waqas@astro.princeton.edu) - Jan 2013

This converts a FITS file or a list of FITS to JPEGs. Be default, tries to be
smart about the FITS image type. For FITS with imagetypes of 'object' or
'focus', it makes 3 x 3 stamp JPEGs of the image. For FITS of any other
imagetype, it makes full frame JPEGS.

'''

import optparse
import os
import os.path
import sys

import imageutils

if __name__ == '__main__':

    # parse the command line
    usage = 'Usage: %prog [options] [FITS file to convert] [output filename]'
    options = optparse.OptionParser(usage=usage)
    options.add_option('-t','--type',
                       dest='out_type',
                       action='store',
                       type='string',
                       default=None,
                       help=('type of output to create: is one of [stamp, full]. '
                             'by default, fits2jpeg figures this out automatically '
                             'based on the IMAGETYP keyword in the FITS header.'))
    opts, args = options.parse_args()

    if len(args) < 1:
        options.error('need a FITS file to work on')
    elif len(args) < 2:
        options.error('need an output filename')

    in_file, out_file = args

    if opts.out_type is None:
        imageutils.fits_to_full_jpeg(in_file,out_fname=out_file)
    elif opts.out_type == 'stamp':
        imageutils.fits_to_stamps_jpeg(in_file,out_fname=out_file)
    elif opts.out_type == 'full':
        imageutils.fits_to_full_jpeg(in_file,out_fname=out_file)
    else:
        print('unknown output type %s' % opts.out_type)
        sys.exit(1)

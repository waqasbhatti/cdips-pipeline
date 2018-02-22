#!/usr/bin/env python

'''
fitshdr.py - Waqas Bhatti (wbhatti@astro.princeton.edu) - Jan 2013

Simple utility to read and display the image header for a FITS file. Gets rid of
blank lines in the header by default, set the --full command line option to show
these as well.

'''

import optparse
import os
import os.path
import sys

import astropy.io.fits as pyfits


if __name__ == '__main__':

    # parse the command line
    usage = 'Usage: %prog [options] [FITS file to read]'
    options = optparse.OptionParser(usage=usage)
    options.add_option('-f','--full',
                       dest='full_header',
                       action='store_true',
                       default=False,
                       help='show the full header, including all blank lines')
    options.add_option('-x','--ext',
                       dest='read_ext',
                       action='store',
                       type='int',
                       default=None,
                       metavar='EXT',
                       help='show header for extension EXT')
    options.add_option('-k','--keyword',
                       dest='keyword',
                       action='store',
                       type='string',
                       default=None,
                       metavar='KEYWORD',
                       help='show value of FITS keyword KEYWORD')
    opts, args = options.parse_args()
    if len(args) < 1:
        options.error('need at least one FITS file to work on')

    # go through each file in the argument list
    for fits_file in args:

        hdu = pyfits.open(fits_file)
        n_ext = len(hdu)

        if opts.read_ext is not None:

            ext = hdu[opts.read_ext]

            if opts.keyword is not None:

                try:
                    keyword_val = ext.header[opts.keyword]
                    print('%s: %s'% (fits_file, keyword_val))

                except KeyError:
                    print('keyword %s not found' % opts.keyword)

            else:

                header = ext.header.tostring(padding=False,
                                             sep='\n',
                                             endcard=False)

                header = header.split('\n')
                header = [x.strip() for x in header]

                if not opts.full_header:
                    header = [x for x in header if len(x) > 0]

                print('\nfile: %s, header for extension %i:' % (fits_file,
                                                                opts.read_ext))
                for card in header:
                    print(card)

        else:

            for i,ext in enumerate(hdu):

                data_length = ext.fileinfo()['datSpan']

                if data_length > 0:

                    if opts.keyword is not None:

                        try:
                            keyword_val = ext.header[opts.keyword]
                            print('%s: %s'% (fits_file, keyword_val))

                        except KeyError:
                            print('keyword %s not found' % opts.keyword)


                    else:

                        header = ext.header.tostring(padding=False,
                                                     sep='\n',
                                                     endcard=False)

                        header = header.split('\n')
                        header = [x.strip() for x in header]

                        if not opts.full_header:
                            header = [x for x in header if len(x) > 0]

                        print('\nfile: %s, header for extension %i:' %
                              (fits_file,i))
                        for card in header:
                            print(card)

        hdu.close()


#!/usr/bin/env python

"""
epdlightcurves.py

Run EPD on lightcurves.

Uses a flexible model to allow use of additional external parameters instead
of hard-baking the exact form of the EPD function.
"""

import os
import sys
from glob import glob
from datetime import datetime
from collections import OrderedDict

import numpy as np
import scipy
from scipy.stats import pearsonr
from scipy.signal import medfilt
from scipy.stats import sigmaclip
from scipy.linalg import lstsq

import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import psycopg2 as pg

import multiprocessing as mp
from functools import partial

# DB Login Info
PGPASSFILE = os.path.expanduser('~/.pgpass')
PGUSER = 'hpx2'
PGDATABASE = 'hpx2'
PGHOST = 'xphtess1'

# Fields to get from DB
DEFAULT_DBPARAMS = ['ha', 'z', 'moondist', 'moonelev', 'moonph', 'wind',
                    'humidity', 'skytdiff', 'ambtemp']

# EPD parameters. Should be dict-like object with keys = column name,
# values = None or a function to generate that column.
DEFAULT_EPDPARAMS = OrderedDict([
    ('fsv', None),
    ('fdv', None),
    ('fkv', None),
    ('fsv2', lambda x: x.fsv**2),
    ('fdv2', lambda x: x.fdv**2),
    ('fkv2', lambda x: x.fkv**2),
    ('fsd', lambda x: x.fsv*x.fdv),
    ('fsk', lambda x: x.fsv*x.fkv),
    ('fdk', lambda x: x.fdv*x.fkv),
    ('xis2', lambda x: np.sin(2*np.pi*x.xic)),
    ('xic2', lambda x: np.cos(2*np.pi*x.xic)),
    ('yis2', lambda x: np.sin(2*np.pi*x.yic)),
    ('yic2', lambda x: np.cos(2*np.pi*x.yic)),
    ('moonelev', None),
    ('cosmoonelev', lambda x: np.cos(x.moonelev * np.pi/180)),
    ('sinmoonelev', lambda x: np.sin(x.moonelev * np.pi/180)),
    ('cos4moondist', lambda x: np.cos(4*x.moondist * np.pi/180)),
    ('sin4moondist', lambda x: np.sin(4*x.moondist * np.pi/180)),
    ('wind', None),
    ('wind2', lambda x: x.wind**2),
    ('humidity', None),
    ('hum2', lambda x: x.humidity**2),
    ('skytdiff', None),
    ('cosha', lambda x: np.cos(x.ha * np.pi/12)),
    ('sinha', lambda x: np.sin(x.ha * np.pi/12)),
    ('cosz', lambda x: np.cos(x.z * np.pi/180)),
    ('sinz', lambda x: np.sin(x.z * np.pi/180)),
])


#######################
#  Utility Functions  #
#######################
def read_lc(lcfile, cols=None, colnames=None,
            epmags=False, tfmags=False):
    '''
    Reads the requested columns in the lightcurve file.

    Args:
        cols: List of columns to read, or None (reads all columns)
        colnames: List of column names
    '''
    # Default columns in grcollectilc file
    LCCOLS = ['rjd', 'rstfc', 'id',
              'xcc', 'ycc', 'xic', 'yic',
              'fsv', 'fdv', 'fkv', 'bgv', 'bev',
              'fl1', 'fe1', 'rm1', 're1', 'rq1',
              'fl2', 'fe2', 'rm2', 're2', 'rq2',
              'fl3', 'fe3', 'rm3', 're3', 'rq3']
    if epmags:
        LCCOLS += ['ep1', 'ep2', 'ep3']
    if tfmags:
        LCCOLS += ['tf1', 'tf2', 'tf3']
    LCCOLNUMS = [0, 1, 5, 6, 7, 8, 9,
                 14, 15, 19, 20, 24, 25]

    if cols is None:
        cols = LCCOLNUMS
        colnames = np.asarray(LCCOLS)[LCCOLNUMS]
    elif cols == 'all':
        cols = list(range(len(LCCOLS)))
        colnames = LCCOLS

    lc = pd.read_csv(lcfile, delim_whitespace=True, comment='#',
                     header=None, usecols=cols, names=colnames)
    lc = lc.reset_index(drop=True)

    return lc


def connect_db(pgpassfile=PGPASSFILE, pguser=PGUSER, pgdatabase=PGDATABASE,
               pghost=PGHOST):
    '''
    Connect to the database.
    '''
    if os.path.exists(PGPASSFILE):
        with open(PGPASSFILE) as infd:
            pgpass_contents = infd.readlines()
            pgpass_contents = [x.split(':') for x in pgpass_contents if not
                               x.startswith('#')]
            PGPASSWORD = [x[-1] for x in pgpass_contents
                          if (x[0] == PGHOST and
                              x[2] == PGDATABASE and
                              x[3] == PGUSER)]
            PGPASSWORD = PGPASSWORD[0].strip('\n')
    else:
        print('ERR! %sZ: could not find postgres database credientials at %s' %
              (datetime.utcnow().isoformat(), PGPASSFILE))

    return pg.connect(user=PGUSER, password=PGPASSWORD,
                      database=PGDATABASE, host=PGHOST)


def get_dbcols(paramlist, fieldinfo, inheader=True):
    '''
    Get additional external parameters from database.

    Args:
        paramlist: List of parameters to get.
        inheader: If true, gets params from the fitsheader json object.
    '''
    if inheader:
        paramstr = ', '.join(["fitsheader->'%s'" % p.upper() for p in paramlist])
    else:
        paramstr = ', '.join(paramlist)

    # Connect to database
    db = connect_db()
    cursor = db.cursor()
    query = "select fits, " + paramstr + \
            (" from calibratedframes where (fitsheader->'PROJID' = '%s')"
             % fieldinfo['projectid'])
    cursor.execute(query)
    rows = cursor.fetchall()
    db.close()

    dbcols = pd.DataFrame(rows, columns=['fits'] + paramlist)
    dbcols.loc[:, 'rstfc'] = dbcols.fits.map(
        lambda f: os.path.splitext(os.path.basename(f))[0])

    return dbcols

def generate_newcols(lc, funcs):
    '''
    Add new columns to the dataframe
    '''
    newseries = []
    for n, f in funcs.items():
        if f is None:
            continue
        if n in lc:
            lc.loc[:, n] = f(lc)
        else:
            newseries.append(pd.Series(f(lc), name=n, index=lc.index))
    if len(newseries) > 0:
        lc = pd.concat([lc] + newseries, axis=1)
    return lc

# Note: this should probably be in a different file
def get_new_filelist(infilelist, outext=None, outfilefunc=None,
                     outfileglob=None, outdir=None):
    '''
    Utility function that takes an input list of files, finds the expected
    output files, determines which have already been generated,
    and returns those input files for which the task needs to be run on.

    Args:
        infilelist: Input file list
        outext: Extension of output file, or
        outfilefunc: Callable that generates output filename from each
            input file. Takes precedence over outext.
            outext and outfilefunc cannot both be None.
        outfileglob: glob for outfile. If None, assumes it is '*.outext'.
        outdir: Output file directory. If None, assumes files are in same
            directory as infilelist.

    Returns:
        newfilelist: Files that task should be run on.
    '''
    if outdir is None:
        outdir = os.path.dirname(infilelist[0])

    infiles = list(map(os.path.basename, infilelist))

    # Get list of requested output files
    if outfilefunc is not None:
        outfilelist = list(map(outfilefunc, infiles))
    elif outext is not None:
        outfilelist = list(map(lambda f: os.path.splitext(f)[0] + outext,
                               infiles))
    else:
        raise ValueError

    # Get list of already existing files
    if outfileglob is None:
        if outext is not None:
            outfileglob = '*' + outext
        else:
            outfileglob = '*' + os.path.splitext(outfilelist[0])[-1]
    alreadyexists = glob(os.path.join(outdir, outfileglob))
    alreadyexists = list(map(os.path.basename, alreadyexists))

    # Find differences between the two lists
    ae_mask = np.isin(outfilelist, alreadyexists)
    filelist = [f[0] for f in zip(infilelist, ae_mask) if not f[1]]

    return filelist


#############
#  Run EPD  #
#############
def run_epd(lcfile, epdparams=None, lccols=None, lccolnames=None, addcols=None,
            framecol='rstfc', sigclip=5.0, smooth=21, minndet=200,
            num_aps=3, outdir=None, outext='.epdlc'):
    '''
    Runs EPD on the given lcfile

    Args:
        lcfile (str): Path to lcfile.
        epdparams (dict or OrderedDict): keys are columns,
            values are either None or a callable that operates on a lightcurve
        lccols (list): Columns to read from the lc
        lccolnames (list): Names of columns read from the lc
        addcols (DataFrame): Additional external parameters
        framecol (str): Column to merge lc and addcols on

        sigclip (float): Sigma-clipping to be done on lc before fitting
        smooth (int or False): Smoothing parameter to apply to lc
        minndet (int): Minimum number of observations in lc
        num_aps (int): Number of apertures
        outext (str): Extension of output file.

    Returns:
        outfile if successful, None otherwise.
    '''
    # Defaults
    if epdparams is None:
        epdparams = DEFAULT_EPDPARAMS

    # Read LC
    lc = read_lc(lcfile, lccols, lccolnames)
    # Join with additional columns
    if addcols is not None:
        lc = pd.merge(lc, addcols, how='left', on=framecol)
    # Generate additional columns required for EPD
    lc = generate_newcols(lc, epdparams)
    epdcols = list(epdparams.keys())

    # Only use rows where all params are finite
    finite_params = np.ones(len(lc), dtype=bool)
    for col in epdcols:
        finite_params &= np.isfinite(lc[col])

    failedall = True
    # Run for each aperture
    for ap in range(1, num_aps+1):
        magcol = 'rm%d' % ap
        finiteind = np.isfinite(lc[magcol]) & finite_params
        if np.sum(np.ones_like(finiteind)[finiteind]) < minndet:
            failedall &= True
            continue
        else:
            failedall = False

        mag_median = np.nanmedian(lc[magcol])
        mag_stdev = np.nanstd(lc[magcol])

        # Sigma clip
        if sigclip:
            excludeind = abs(lc[magcol] - mag_median) < (sigclip * mag_stdev)
            finalind = finiteind & excludeind
        else:
            finalind = finiteind

        # Smoothing
        if smooth:
            mags = medfilt(lc[magcol][finalind], smooth)
        else:
            mags = lc[magcol][finalind]

        ### Perform fitting
        # Construct the matrix from the EPD columns (and a constant term)
        epdmatrix = np.c_[lc[epdcols], np.ones(len(lc))]
        # Use least squares fitting to only chosen rows
        coeffs, residuals, rank, singulars = lstsq(epdmatrix[finalind],
                                                   mags)
        # Now compute full EPD mags
        lc.loc[:,'ep%d' % ap] = (lc[magcol] - np.dot(coeffs, epdmatrix.T)
                                 + mag_median)

    # Save lcfile if at least one epdcolumn was generated
    lcid = os.path.splitext(os.path.basename(lcfile))[0]
    if failedall:
        print('%sZ: EPD failed for %s.' % (datetime.utcnow().isoformat(),
                                           lcid + '.grcollectilc'))
        return None
    else:
        if outdir is None:
            outdir = os.path.dirname(lcfile)
        outfile = os.path.join(outdir, lcid + outext)

        # Proper file output by re-reading lines from input file
        inf = open(lcfile, 'rb')
        inflines = inf.readlines()
        inf.close()

        outf = open(outfile, 'wb')
        epcols = ['ep%d' % (i+1) for i in range(num_aps)]
        for line, row in zip(inflines, lc.iterrows()):
            epmags = [('%.6f' % row[1][col]) for col in epcols]
            inline = line.decode('utf-8').rstrip('\n')
            outline = ' '.join([inline] + epmags) + '\n'
            outf.write(outline.encode('utf-8'))
        outf.close()

        print('%sZ: EPD OK for %s.' % (datetime.utcnow().isoformat(),
                                       lcid + '.grcollectilc'))
        return outfile


#################
#  Parallelize  #
#################
def parallel_run_epd(lcdir, fieldinfo, lcglob='*.grcollectilc', outdir=None,
                     epdparams=None, dbparams='default', framecol='rstfc',
                     lccols=None, lccolnames=None,
                     sigclip=5.0, smooth=21, minndet=200, num_aps=3,
                     overwrite=False, nworkers=2, maxworkertasks=1000):
    '''
    Run EPD in parallel. Gets columns from DB only once.
    '''
    lcfiles  = glob(os.path.join(lcdir, lcglob))
    if overwrite:
        lcfiles = get_new_filelist(lcfiles, outext=outext)

    if outdir is not None and not os.path.isdir(outdir):
        os.mkdir(outdir)

    # For testing
    lcfiles = lcfiles[:5]

    # Get additional parameter columns from database
    if dbparams == 'default':
        dbparams = DEFAULT_DBPARAMS
    if dbparams is not None:
        addcols = get_dbcols(dbparams, fieldinfo)
    else:
        addcols = None

    run_epd_worker = partial(run_epd, epdparams=epdparams, lccols=lccols,
                             lccolnames=lccolnames, addcols=addcols,
                             framecol=framecol, sigclip=sigclip, smooth=smooth,
                             minndet=minndet, num_aps=num_aps,
                             outdir=outdir, outext='.epdlc')

    pool = mp.Pool(nworkers, maxtasksperchild=maxworkertasks)
    results = pool.map(run_epd_worker, lcfiles)
    pool.close()
    pool.join()

    print('%sZ: EPD done. %s lightcurves processed, %s successful.' %
          (datetime.utcnow().isoformat(), len(lcfiles),
           sum(l is not None for l in results)))

    return results




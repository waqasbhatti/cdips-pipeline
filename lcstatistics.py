# -*- coding: utf-8 -*-
from __future__ import division, print_function
"""
functions for assessing statistics of lightcurves you've made.

Reading functions
    read_tfa_lc:
        read TFA lightcurve

    read_acf_stat_files:
        read stack of csv files with autocorrelation functions evaluated at
        selected time lags.

    read_stats_file:
        read compiled epd or tfa stat files (originally from aperturephot)

Stat calculation functions
    compute_acf_statistics_worker:
        worker to compute autocorrelation function statistics

    parallel_compute_acf_statistics:
        compute autocorrelation function stats for many lightcurves

    compute_dacf:
        compute the discrete autocorrelation function for a lightcurve

    compute_correction_coeffs:
        fit catalog mags to fluxes to obtain a corrected catalog mag

    compute_lc_statistics:
        worker to compute median, MAD, mean, stdev for each lightcurve

    parallel_compute_lc_statistics:
        compute above statistics for a directory of lightcurves and compile

    compute_lc_statistics_fits:
        alternative to compute_lc_statistics, but only for FITS LCs. (has extra
        flexibility in column choice)

Plotting functions
    percentiles_RMSorMAD_stats_and_plots:
        make csv files and (optionally) percentiles plots of MAD vs magnitude.

    acf_percentiles_stats_and_plots:
        make csv files & plots summarizing ACF statistics for many stars

    plot_raw_epd_tfa:
        Plot a 2 (or 3) row, 1 column plot with rows of:
            * raw mags vs time
            (* EPD mags vs time)
            * TFA mags vs time.

    plot_lightcurve_and_ACF:
        Plot a 3 row, 2 column plot with rows of:
            * raw mags vs time (and ACF)
            * EPD mags vs time (and ACF)
            * TFA mags vs time. (and ACF)

    plot_raw_tfa_bkgd:
        Plot a 3 row, 1 column plot with rows of:
            * raw mags vs time
            * TFA mags vs time.
            * background vs time

    plot_raw_tfa_bkgd_fits:
        Wrapper to above, given a *llc.fits lightcurve.

TODO:
move items from aperturephot.py here, and update all respective calls.

"""

from astrobase.periodbase import macf
from astrobase.varbase.autocorr import autocorr_magseries

import os, pickle, itertools
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import aperturephot as ap
from glob import glob
from astropy import units as units, constants as constants
from datetime import datetime

from numpy import array as nparr
from scipy.stats import sigmaclip as stats_sigmaclip

from datetime import datetime
import multiprocessing as mp
from functools import partial

from astropy.io import fits

#####################
# READING FUNCTIONS #
#####################

def read_tfa_lc(tfafile,
                jdcol=0,
                timename='btjd',
                rmcols=[14,19,24],
                epcols=[27,28,29],
                tfcols=[30,31,32],
                recols=[15,20,25]):
    """
    Get the times and mags from a .tfalc file. For speed, don't read any of the
    rest. The full contents of a .tfalc file are written below.

    If error columns are None, they are not read.

    Args:
        tfafile (str): path to TFA lightcurve file

        jdcol (int): integer index of time column

        timename (str): name to assign time column

        rmcols (list of ints): integer indices of RAW magnitude columns.
        epcols (list of ints): integer indices of EPD magnitude columns.
        tfcols (list of ints): integer indices of TFA magnitude columns.

        recols (list of ints, or None): indices of RAW magnitude error columns.
        (EPD and TFA magnitudes currently do not have errors; they are
        presumably inherited from RAW magnitudes).

    Returns:
        np.ndarray of times and magnitudes.

    ------------------------------------------
    00 rjd    Barycentric TESS Julian Date (BTJD := BJD - 2457000) if TESS.
              Reduced Julian Date (RJD = JD - 2400000.0) if HATPI.
    01 rstfc  Unique frame key ({STID}-{FRAMENUMBER}_{CCDNUM})
    02 hat    HAT ID of the object
    03 xcc    Original X coordinate on CCD on photref frame
    04 ycc    Original y coordinate on CCD on photref frame
    05 xic    Shifted X coordinate on CCD on subtracted frame
    06 yic    Shifted Y coordinate on CCD on subtracted frame
    07 fsv    Measured S value
    08 fdv    Measured D value
    09 fkv    Measured K value
    10 bgv    Background value
    11 bge    Background measurement error

    12 ifl1   Flux in aperture 1 (ADU)
    13 ife1   Flux error in aperture 1 (ADU)
    14 irm1   Instrumental magnitude in aperture 1
    15 ire1   Instrumental magnitude error for aperture 1
    16 irq1   Instrumental magnitude quality flag for aperture 1 (0/G OK, X bad)

    17 ifl2   Flux in aperture 2 (ADU)
    18 ife2   Flux error in aperture 2 (ADU)
    19 irm2   Instrumental magnitude in aperture 2
    20 ire2   Instrumental magnitude error for aperture 2
    21 irq2   Instrumental magnitude quality flag for aperture 2 (0/G OK, X bad)

    22 ifl3   Flux in aperture 3 (ADU)
    23 ife3   Flux error in aperture 3 (ADU)
    24 irm3   Instrumental magnitude in aperture 3
    25 ire3   Instrumental magnitude error for aperture 3
    26 irq3   Instrumental magnitude quality flag for aperture 3 (0/G OK, X bad)

    27 ep1    EPD magnitude for aperture 1
    28 ep2    EPD magnitude for aperture 2
    29 ep3    EPD magnitude for aperture 3

    30 tf1    TFA magnitude for aperture 1
    31 tf2    TFA magnitude for aperture 2
    32 tf3    TFA magnitude for aperture 3
    ------------------------------------------
    """

    lcmagcols = [pre+str(ix)
                 for pre in ['RM','EP','TF'] for ix in range(1,4)]
    lcerrcols = [pre+'ERR'+str(ix)
                 for pre in ['RM'] for ix in range(1,4)]

    if (
        not isinstance(recols,list)
    ):
        lcdata = np.genfromtxt(tfafile,
                               usecols=tuple(
                                   [jdcol] + rmcols + epcols + tfcols),
                               names=[timename] + lcmagcols)
    elif (
        isinstance(recols,list)
    ):
        lcdata = np.genfromtxt(tfafile,
                               usecols=tuple(
                                   [jdcol]
                                   + rmcols + epcols + tfcols
                                   + recols),
                               names=[timename] + lcmagcols + lcerrcols)

    else:
        raise ValueError(
            'Make consistent choice of including or not including errors')

    return lcdata


def read_acf_stat_files(acfstatfiles, N_acf_types=9):
    """
    return dataframe of length (N_stars*N_lag_times per star) with evaluated
    autocorrelation functions for each aperure and detrending step.
    """

    if not acfstatfiles:
        raise AssertionError('cannot assess run if there are no ACF statfiles')

    iscsvs = nparr([f.endswith('.csv') for f in acfstatfiles])
    if not np.all(iscsvs):
        raise ValueError('expected only csv files')

    # read acf files
    df = pd.concat([pd.read_csv(f) for f in acfstatfiles],
                   ignore_index = True)

    return df

def read_stats_file(statsfile, fovcathasgaiaids=False):
    """
    Reads the stats file into a numpy recarray.
    """

    if fovcathasgaiaids:
        idstrlength = 19
    else:
        idstrlength = 17

    cols = []
    fmts = ''
    with open(statsfile) as f:
        # Find start of comments
        line = f.readline()
        while 'column' not in line:
            line = f.readline()
        line = f.readline()

        # Get column names
        while line.startswith('#'):
            colnames = line.split(':')[1].strip()
            if colnames.startswith('sigma'):
                sigclip = True
                colnames = colnames[14:]
            else:
                sigclip = False
            colnames = [c.strip() for c in colnames.split(',')]

            for c in colnames:
                if 'object' in c:
                    cols.append('lcobj')
                    fmts += 'U{:d},'.format(idstrlength)
                elif 'catalog' in c:
                    cols.append('catmag')
                    fmts += 'f8,'
                elif 'corrected' in c:
                    cols.append('corrmag_' + c[-3:].lower())
                    fmts += 'f8,'
                else:
                    if sigclip:
                        cols.append(c.lower().replace(' ', '_sigclip_'))
                    else:
                        cols.append(c.lower().replace(' ', '_'))
                    if 'ndet' in c:
                        fmts += 'i8,'
                    else:
                        fmts += 'f8,'
            line = f.readline()

    stats = np.genfromtxt(statsfile, dtype=fmts, names=cols)

    #  stats = np.genfromtxt(
    #      statsfile,
    #      dtype=(
    #          'U{:d},f8,'
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM1
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM2
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM3
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP1
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP2
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP3
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF1
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF2
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF3
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF1
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF2
    #          'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF3
    #          'f8,f8,f8'.format(idstrlength)    # corrmags
    #      ),
    #      names=[
    #          'lcobj','catmag',
    #          'median_rm1','mad_rm1','mean_rm1','stdev_rm1','ndet_rm1',
    #          'median_sigclip_rm1','mad_sigclip_rm1','mean_sigclip_rm1',
    #          'stdev_sigclip_rm1','ndet_sigclip_rm1',
    #          'median_rm2','mad_rm2','mean_rm2','stdev_rm2','ndet_rm2',
    #          'median_sigclip_rm2','mad_sigclip_rm2','mean_sigclip_rm2',
    #          'stdev_sigclip_rm2','ndet_sigclip_rm2',
    #          'median_rm3','mad_rm3','mean_rm3','stdev_rm3','ndet_rm3',
    #          'median_sigclip_rm3','mad_sigclip_rm3','mean_sigclip_rm3',
    #          'stdev_sigclip_rm3','ndet_sigclip_rm3',
    #          'median_ep1','mad_ep1','mean_ep1','stdev_ep1','ndet_ep1',
    #          'median_sigclip_ep1','mad_sigclip_ep1','mean_sigclip_ep1',
    #          'stdev_sigclip_ep1','ndet_sigclip_ep1',
    #          'median_ep2','mad_ep2','mean_ep2','stdev_ep2','ndet_ep2',
    #          'median_sigclip_ep2','mad_sigclip_ep2','mean_sigclip_ep2',
    #          'stdev_sigclip_ep2','ndet_sigclip_ep2',
    #          'median_ep3','mad_ep3','mean_ep3','stdev_ep3','ndet_ep3',
    #          'median_sigclip_ep3','mad_sigclip_ep3','mean_sigclip_ep3',
    #          'stdev_sigclip_ep3','ndet_sigclip_ep3',
    #          'median_tf1','mad_tf1','mean_tf1','stdev_tf1','ndet_tf1',
    #          'median_sigclip_tf1','mad_sigclip_tf1','mean_sigclip_tf1',
    #          'stdev_sigclip_tf1','ndet_sigclip_tf1',
    #          'median_tf2','mad_tf2','mean_tf2','stdev_tf2','ndet_tf2',
    #          'median_sigclip_tf2','mad_sigclip_tf2','mean_sigclip_tf2',
    #          'stdev_sigclip_tf2','ndet_sigclip_tf2',
    #          'median_tf3','mad_tf3','mean_tf3','stdev_tf3','ndet_tf3',
    #          'median_sigclip_tf3','mad_sigclip_tf3','mean_sigclip_tf3',
    #          'stdev_sigclip_tf3','ndet_sigclip_tf3',
    #          'median_rf1','mad_rf1','mean_rf1','stdev_rf1','ndet_rf1',
    #          'median_sigclip_rf1','mad_sigclip_rf1','mean_sigclip_rf1',
    #          'stdev_sigclip_rf1','ndet_sigclip_rf1',
    #          'median_rf2','mad_rf2','mean_rf2','stdev_rf2','ndet_rf2',
    #          'median_sigclip_rf2','mad_sigclip_rf2','mean_sigclip_rf2',
    #          'stdev_sigclip_rf2','ndet_sigclip_rf2',
    #          'median_rf3','mad_rf3','mean_rf3','stdev_rf3','ndet_rf3',
    #          'median_sigclip_rf3','mad_sigclip_rf3','mean_sigclip_rf3',
    #          'stdev_sigclip_rf3','ndet_sigclip_rf3',
    #          'corrmag_ap1','corrmag_ap2','corrmag_ap3',
    #      ]
    #  )

    return stats

#####################################
# FUNCTIONS TO CALCULATE STATISTICS #
#####################################

def compute_acf_statistics_worker(task, n_apertures=3, timename='TMID_BJD',
                                  filterwindow=7, istessffi=True,
                                  isfitslc=True):
    # dtrtypes is e.g., ['IRM','PCA','TFA']

    #NOTE : might want to check a couple smoothing values ("filterwindows")...

    try:
        if not istessffi:
            raise NotImplementedError(
            'this function assumes an ffi cadence of 30 minutes to evaluate ACFs. '
            'to generalize this you just need to interpolate, but I am lazy.'
            )

        tfafile, outdir, eval_times_hr, dtrtypes = task

        if not isfitslc:
            lcdata = read_tfa_lc(tfafile)
        else:
            hdulist = fits.open(tfafile)
            lcdata = hdulist[1].data

        outpickle = os.path.join(
            outdir,
            os.path.basename(tfafile).replace('.fits','_acf_stats.pickle')
        )
        outcsv = os.path.join(
            outdir,
            os.path.basename(tfafile).replace('.fits','_acf_stats.csv')
        )
        if os.path.exists(outpickle) and os.path.exists(outcsv):
            return 1

        d_pkl, outdf = {}, pd.DataFrame({})
        flux_d = {}
        acfs = {}
        for ap in range(1,n_apertures+1):

            for dtrtype in dtrtypes:

                apname = dtrtype+'{}'.format(ap)
                errname = 'IRE{}'.format(ap)

                time = lcdata[timename]

                flux_d[dtrtype] = lcdata[apname]

                err = lcdata[errname]

                acfs[dtrtype] = autocorr_magseries(time, flux_d[dtrtype], err,
                                                   maxlags=None,
                                                   fillgaps='noiselevel',
                                                   sigclip=5,
                                                   magsarefluxes=True,
                                                   filterwindow=filterwindow)

            if not np.isclose(acfs[dtrtype]['cadence']*24*60, 30, atol=1e-3):
                raise NotImplementedError(
                'this function assumes an ffi cadence of 30 minutes to evaluate ACFs. '
                'to generalize this you just need to interpolate, but I am lazy.'
                )

            apstr = 'AP{}'.format(ap)
            d_pkl[apstr] = {}

            for k, v in acfs.items():

                dtr = k.lower()
                acf = v

                d_pkl[apstr]['acf_{}'.format(dtr)] = acf['acf']
                d_pkl[apstr]['lag_time_{}'.format(dtr)] = acf['lags']*acf['cadence']
                d_pkl[apstr]['itimes_{}'.format(dtr)] = acf['itimes']
                d_pkl[apstr]['ifluxs_{}'.format(dtr)] = acf['imags']
                d_pkl[apstr]['ierrs_{}'.format(dtr)] = acf['ierrs']
                d_pkl[apstr]['cadence'] = acf['cadence']

            # assuming FFI cadence of 30 minutes, evalute ACF at desired lags.
            n_acf_vals = len(acfs[k]['acf'])

            # in "TUNE" mode, might want short-timescale ACFs (otherwise,
            # fails)
            if max(2*eval_times_hr)>n_acf_vals:
                eval_times_hr = np.array([1,2,6,12])
            else:
                pass

            eval_times_hr = nparr(eval_times_hr)
            outdf['LAG_TIME_HR'] = eval_times_hr

            colstrs = ['{}{}'.format(dtr, ap) for dtr in dtrtypes]

            for colstr, dtr in zip(colstrs, dtrtypes):

                # wonky indexing scheme. 30 minute cadence -> e.g., 1hr lag is
                # index number 2.
                these_acf_vals = acfs[dtr]['acf'][2*eval_times_hr]

                outdf[colstr+"_ACF"] = these_acf_vals

        with open(outpickle, 'wb') as f:
            pickle.dump(d_pkl, f, pickle.HIGHEST_PROTOCOL)
        print('wrote {}'.format(outpickle))

        outdf.to_csv(outcsv, index=False)
        print('wrote {}'.format(outcsv))

        return 1

    except Exception as e:
        print('{} failed, error was {}'.format(repr(task), repr(e)))


def parallel_compute_acf_statistics(
    tfafiles, outdir, eval_times_hr=[1,2,6,12,24,48,60,96,120,144,192],
    nworkers=16, maxworkertasks=1000,
    dtrtypes=['IRM','PCA','TFA']
):
    """
    Given list of TFA lightcurves, calculate autocorrelation functions and
    evaluate them at specific times, e.g., 1hr, 2hr, 6hr, 24hr.

    Gives two outputs:
        1) a pickle file with the autocorrelation function results.
        2) a csv file with the summarized ACF values at the desired lag times.

    Args:
        tfafiles (list of strs): list of paths to TFA files

        outdir (str): directory where output csv files with ACF statistics
        are written.

        eval_times_hr (list of ints): times at which to evaluate the
        autocorrelation functions, in hours.

        dtrtypes (list of strs): stages of detrending for which to evaluate the
        ACFs.
    """

    print('%sZ: %s files to compute ACF statistics for' %
          (datetime.utcnow().isoformat(), len(tfafiles)))

    tasks = [(x, outdir, eval_times_hr, dtrtypes) for x in tfafiles]

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    # fire up the pool of workers
    results = pool.map(compute_acf_statistics_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {result for result in results}


def compute_dacf(times, mags, errs, lagstep=0.1,
                 lagmin=0.0, lagmax=10.0, laglist=None):
    """
    Compute the discrete autocorrelation function for a lightcurve.

    Follows the definition of DACF in Edelson & Krolik (1988).
    This works better for ground-based (HATPI) data with large data gaps,
    avoiding the need to fill gaps with white noise as in astrobase.
    Equivalent to the VARTOOLS implementation, but can run on multiple
    LC columns simultaneously.

    Args:
        times (array):  Array of times
        mags (array):   Magnitude column(s)
        errs (array):   Error column(s). There should either by one errcol
                        for each magcol, or mags should be a list of magcols
                        where each sublist corresponds to one errcol.
        lagmin (float): Minimum lag to evaluate
        lagmax (float): Maximum lag
        lagstep (float):Size of lag bins
        laglist (list): Only evaluate the lag at these points. Takes priority
                        over lagmin, lagmax.
    """
    times = np.asarray(times)
    mags = np.asarray(mags)
    if type(errs) != str:
        errs = np.asarray(errs)

    # Check dimensions
    if mags.ndim == 1:
        mags = np.asarray([mags])
    # List of magcols provided
    if mags.ndim == 2:
        # Use stdev as errs
        if type(errs) == str and errs == 'correcterrs':
            errs = [np.nanstd(magcol) for magcol in mags]
            errs = np.asarray([np.full(mags.shape[1], err) for err in errs])
            err_inds = [i for i in range(mags.shape[0])]
        # One err col
        elif errs.ndim == 1:
            errs = np.asarray([errs])
            err_inds = [0 for i in range(mags.shape[0])]
        # List of errcols
        elif errs.ndim == 2:
            assert errs.shape[0] == mags.shape[0]
            err_inds = [i for i in range(mags.shape[0])]
    # List of list of magcols
    if mags.ndim == 3:
        if type(errs) == str and errs == 'correcterrs':
            errs = [np.nanstd(magcol) for maglvl in mags for magcol in maglvl]
            errs = np.asarray([np.full(mags.shape[2], err) for err in errs])
            err_inds = list(range(mags.shape[0] * mags.shape[1]))
        # Must provide list of errcols
        else:
            assert errs.ndim == 2 and errs.shape[0] == mags.shape[0]
            err_inds = []
            for i in range(errs.shape[0]):
                err_inds += [i for j in range(len(mags[i]))]
        # Flatten magcols
        mags = mags.reshape(-1, mags.shape[2])

    # Note: will ignore complete rows that have NaNs.
    finite_inds = np.full(len(times), True)
    for i in range(mags.shape[0]):
        finite_inds &= np.isfinite(mags[i])
    for i in range(errs.shape[0]):
        finite_inds &= np.isfinite(errs[i])
    times = times[finite_inds]
    mags = mags[:, finite_inds]
    errs = errs[:, finite_inds]

    if laglist is None:
        if lagmin is None:
            lagmin = 0.0
        if lagmax is None:
            lagmax = np.max(times) - np.min(times)
        nbins = int(np.ceil((lagmax - lagmin)/lagstep) + 1)
        lagarr = lagmin + np.asarray(range(nbins))*lagstep
    else:
        lagmin = 0.0
        nbins = len(laglist)
        lagarr = np.asarray(laglist)

    # Create matrix of bins
    bins = np.empty((len(times), len(times)))
    for i in range(len(times)):
        bins[i, :] = times - lagmin - times[i] + lagstep/2.
    bins = np.floor(bins / lagstep)
    ones = np.ones_like(bins, dtype=int)

    # Pre compute values
    # Error-weighted average
    mags_avg = np.average(mags, weights=1/errs[err_inds]**2, axis=1)
    # Error-weighted mean-centered mag
    mags_ce = (mags - mags_avg[:, None]) / errs[err_inds]
    mags_cesq = mags_ce**2
    # Compute all point-wise correlations and squared errors
    pwc = [mags_ce[i] * mags_ce[i,:, None] for i in range(mags.shape[0])]
    sqe = [mags_cesq[i] + mags_cesq[i,:,None] for i in range(mags.shape[0])]
    # Other empty matrices
    udcf = [np.empty(nbins) for i in range(mags.shape[0])]
    eudcf = [np.empty(nbins) for i in range(mags.shape[0])]
    Nudcf = np.zeros(nbins, dtype=int)

    # Actually compute the cross-correlation
    for (l, lag) in enumerate(lagarr):
        # Find indices which fall into particular bins
        b = np.rint(lag / lagstep)
        flag = (bins == b)
        Nudcf[l] = np.sum(ones[flag])
        for i in range(mags.shape[0]):
            udcf[i][l] = np.sum(pwc[i][flag])
            eudcf[i][l] = np.sum(sqe[i][flag])

    # Normalize by number of points in each bin
    udcf = np.asarray(udcf) / Nudcf
    eudcf = np.sqrt(np.asarray(eudcf)) / Nudcf

    # Only keep lags with non-zero points
    if laglist is None:
        nonzero = Nudcf > 0
        udcf = udcf[:,nonzero]
        eudcf = eudcf[:,nonzero]
        Nudcf = Nudcf[nonzero]
        lagarr = lagarr[nonzero]
        nbins = len(lagarr)

    if udcf.shape[0] == 1:
        udcf = udcf.reshape(udcf.shape[1])
        eudcf = eudcf.reshape(eudcf.shape[1])

    return {
        'udcf': udcf,
        'eudcf': eudcf,
        'Nudcf': Nudcf,
        'nbins': nbins,
        'timestep': lagstep,
        'lags': lagarr
    }


def write_dacf(dacf, outfile, header=None, magtypes=['rm'], num_aps=3):
    """
    Writes a DACF to file.
    """
    outf = open(outfile, 'w')

    # Write column headers
    if header is None:
        header = '#Lag Npairs '
        cols = ['AC_{0:s}{1:d} AC_ERR_{0:s}{1:d}'.format(mag.upper(), apnum)
                for apnum in range(1, num_aps+1) for mag in magtypes]
        header += ' '.join(cols)
    outf.write(header.strip() + '\n')

    for r in range(dacf['nbins']):
        row = '{:.3f} '.format(dacf['lags'][r])
        row += '{:d} '.format(dacf['Nudcf'][r])
        for col in range(dacf['udcf'].shape[0]):
            row += '{:.6f} {:.6f} '.format(
                dacf['udcf'][col, r], dacf['eudcf'][col, r])
        outf.write(row.strip() + '\n')

    outf.close()

    return outfile


def read_dacf(dacffile):
    """
    Reads a dacf object from file.
    """
    dacfarr = np.genfromtxt(dacffile, comments='#', names=True)
    lagarr = dacfarr['Lag']
    nbins = len(lagarr)
    timestep = np.round(np.min(lagarr[1:] - lagarr[:-1]), 3)
    Nudcf = dacfarr['Npairs']
    accols = [c for c in dacfarr.dtype.names if (c.startswith('AC') and
                                                  'ERR' not in c)]
    acerrcols = [c for c in dacfarr.dtype.names if c.startswith('AC_ERR')]
    udcf = []
    eudcf = []
    for ac, err in zip(accols, acerrcols):
        udcf.append(dacfarr[ac])
        eudcf.append(dacfarr[err])
    return {
        'udcf': np.asarray(udcf),
        'eudcf': np.asarray(eudcf),
        'Nudcf': Nudcf,
        'nbins': nbins,
        'timestep': timestep,
        'lags': lagarr
    }


def resample_dacf(dacf, lagstep=0.1, lagmin=0.0, lagmax=10.0):
    """
    Resamples a dacf object
    """
    newlags = np.arange(lagmin, lagmax+1e-5, lagstep)
    idxs = np.searchsorted(newlags, dacf['lags'])
    if dacf['udcf'].ndim == 2:
        newudcf = np.full((dacf['udcf'].shape[0], len(newlags)), np.nan)
        neweudcf = np.full((dacf['eudcf'].shape[0], len(newlags)), np.nan)
        newudcf[:,idxs] = dacf['udcf']
        neweudcf[:,idxs] = dacf['eudcf']
    else:
        newudcf = np.full(len(newlags), np.nan)
        neweudcf = np.full(len(newlags), np.nan)
        newudcf[idxs] = dacf['udcf']
        neweudcf[idxs] = dacf['eudcf']
    newnudcf = np.zeros(len(newlags), dtype=int)
    newnudcf[idxs] = dacf['Nudcf']

    return {
        'udcf': newudcf,
        'eudcf': neweudcf,
        'Nudcf': newnudcf,
        'lags': newlags,
        'timestep': lagstep,
        'nbins': len(newlags)
    }


def compute_dacf_worker(tfafile, outdir='./', lagstep=0.1,
                        lagmin=0.0, lagmax=10.0, laglist=None,
                        n_apertures=3, skipepd=False,
                        writeoutfile=True):
    """
    Wrapper around compute_dacf.
    """
    lcext = os.path.splitext(tfafile)[1]
    if lcext == '.fits':
        raise NotImplementedError
    else:
        lcdata = read_tfa_lc(tfafile)
        if lcdata.size == 0:
            print('%sZ: WRN! %s is empty.' %
                  (datetime.utcnow().isoformat(), tfafile))
            return None

    outfile = os.path.join(outdir,
                           os.path.basename(tfafile).replace(lcext, '.autocorr'))
    header = 'Lag Npairs '

    times = lcdata['btjd']
    mags = []
    errs = []
    for ap in range(1, n_apertures+1):
        rawap = 'RM{:d}'.format(ap)
        epdap = 'EP{:d}'.format(ap)
        tfaap = 'TF{:d}'.format(ap)
        errap = 'RMERR{:d}'.format(ap)

        if not skipepd:
            mags.append([lcdata[rawap], lcdata[epdap], lcdata[tfaap]])
            header += 'AC_RM{0:d} AC_ERR_RM{0:d} AC_EP{0:d} AC_ERR_EP{0:d} '
            header += 'AC_TF{0:d} AC_ERR_TF{0:d} '
        else:
            mags.append(list(lcdata[rawap], lcdata[tfaap]))
            header += 'AC_RM{0:d} AC_ERR_RM{0:d} AC_TF{0:d} AC_ERR_TF{0:d} '
        header = header.format(ap)
        errs.append(lcdata[errap])

    results = compute_dacf(times, mags, errs, lagstep=lagstep,
                           lagmin=lagmin, lagmax=lagmax, laglist=laglist)

    # Write results to output
    if writeoutfile:
        write_dacf(results, outfile, header=header)
    else:
        results['lcobj'] = os.path.splitext(os.path.basename(tfafile))[0]
    return results


def parallel_compute_dacf(lcfiles, outdir=None, lagmin=0.0, lagmax=10.0,
                          lagstep=0.1, n_apertures=3, skipepd=False,
                          compilestats=True, overwrite=True,
                          nworkers=20, workerntasks=1000):
    """
    Computes the ACF for a list of lightcurves.

    Args:
        lcfiles: List of lcfiles.
    """
    if outdir is None:
        outdir = os.path.dirname(lcfiles[0])
    elif not os.path.isdir(outdir):
        os.mkdir(outdir)

    if not overwrite:
        lcfiles = [fn for fn in lcfiles \
                   if os.path.exists(os.path.splitext(fn)[0] + '.autocorr')]
        if len(lcfiles) == 0:
            print('%sZ: autocorrelation already computed for all lcfiles...' %
                  datetime.utcnow().isoformat())
            return None

    print('%sZ: %s lightcurves to calculate complete ACF.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    dacf_worker = partial(compute_dacf_worker, outdir=outdir,
                          lagstep=lagstep, lagmin=lagmin, lagmax=lagmax,
                          n_apertures=n_apertures, skipepd=skipepd)
    results = pool.map(dacf_worker, lcfiles)
    pool.close()
    pool.join()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    return results


def parallel_compute_dacf_statistics(lcdir='./', lcglob='*.tfalc',
                                     lcfiles = None,
                                     outfile=None, lagstep=0.1,
                                     laglist=[0.0, 0.1, 0.3, 1.0, 5.0, 10.0],
                                     n_apertures=3, skipepd=False,
                                     nworkers=20, workerntasks=1000,
                                     outheader=None):
    """
    Evaluates the ACF at fixed positions for lightcurves and compiles results.

    Args:
        lcfiles: List of lcfiles.
    """
    if lcfiles is None:
        lcfiles = glob(os.path.join(lcdir, lcglob))
    print('%sZ: %s lightcurves to process.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    if outfile is None:
        outfile = './acf-statistics.txt'
    outf = open(outfile, 'wb')

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    dacf_stat_worker = partial(compute_dacf_worker, lagstep=lagstep,
                               laglist=laglist, n_apertures=n_apertures,
                               skipepd=skipepd, writeoutfile=False)
    results = pool.map(dacf_stat_worker, lcfiles)
    pool.close()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    # Write header for stats.
    if outheader is not None:
        outf.write(outheader.encode('utf-8'))
    outcolkey = '# object '

    if not skipepd:
        apcols = ['rm', 'ep', 'tf']
    else:
        apcols = ['rm', 'tf']
    num_cols = n_apertures * len(apcols)
    for (lag, apnum, ap) in itertools.product(laglist, range(1, n_apertures+1),
                                              apcols):
        outcolkey += '%s%d_lag_%.2f_day ' % (ap, apnum, lag)
    outcolkey = outcolkey[:-1] + '\n'
    outf.write(outcolkey.encode('utf-8'))

    for stat in results:
        if stat is None:
            continue
        outstr = stat['lcobj'] + ' '
        for i in range(len(laglist)):
            outstr += ('%f ' * num_cols)
            outstr %= tuple(stat['udcf'][:, i])
        outstr = outstr[:-1] + '\n'
        outf.write(outstr.encode('utf-8'))

    outf.close()
    print('%sZ: wrote statistics to file %s' %
          (datetime.utcnow().isoformat(), outfile))

    return results


def compute_correction_coeffs(catmags, fluxes):
    """
    Fit raw fluxes to the catalog mags to correct for too-faint catalog mags,
    according to:
        catmag = -2.5 * log10(flux/c1) + c2

    This function sets c1 = 1 and fits only for c2, essentially the new zero-point.
    Fit is only performed in the bright limit, 8.0 < catmag < 12.0.

    Args:
        catmags (np.array): Catalog mags
        fluxes (1d or 2d np.array): (list of) raw fluxes for each aperture

    Returns:
        correctioncoeffs (2d np.array)
    """
    bright_limit = np.logical_and(catmags < 12.0, catmags > 8.0)
    finite_ind = np.isfinite(catmags) & bright_limit
    fluxes = pd.DataFrame(fluxes)[finite_ind]
    catmags = catmags[finite_ind]

    correctioncoeffs = []
    for ap_col in fluxes:
        try:
            idxs = fluxes[ap_col] > 0.0
            catmags_cut = catmags[idxs]
            c2 = np.polyfit(np.zeros_like(catmags_cut),
                            catmags_cut +
                            2.5*np.log10(fluxes[ap_col][idxs]), 0)[0]
            correctioncoeffs.append([1.0, c2])
        except Exception as e:
            print('ERR! %sZ: could not compute correction coefficients' %
                  datetime.utcnow().isoformat())
            raise e
            return None

    return correctioncoeffs


def compute_lc_statistics(lcfile,
                          rmcols=[14,19,24],
                          epcols=[27,28,29],
                          tfcols=[30,31,32],
                          rfcols=[12,17,22],
                          sigclip=4.0,
                          num_aps=3,
                          epdlcrequired=True,
                          tfalcrequired=False,
                          fitslcnottxt=None,
                          verbose=True):
    """
    Compute mean, median, MAD, stdev magnitudes for a given lightcurve.

    Args:
        lcfile: Either a text .epdlc, .tfalc file or a fits LC.

        rmcols, epcols, tfcols, rfcols: list of magnitude columns for text LC files.

        num_aps (int): Number of apertures. Required for fits LC files.

        epdlcrequired (bool): Use epdlc columns.

        tfalcrequired (bool): Search for tfalc file (if only epdlc provided).

        fitslcnottxt: Does nothing, for maintaining backward compatibility.

    Returns:
        result (dict): Dictionary containing the following fields
            'lcfile': full path to lightcurve
            'lcobj': name of object
            Stat cols:
            Raw magnitudes:
                'median_rm[1-3]', 'mad_rm[1-3]', 'mean_rm[1-3]', 'stdev_rm[1-3]',
                'ndet_rm[1-3]',
                'median_sigclip_rm[1-3]', 'mad_sigclip_rm[1-3]', ...
            And similarly for EPD mags (_ep), TFA mags (_tf), raw fluxes (_rf)
    """
    lcext = os.path.splitext(lcfile)[1]
    # Use dictionary to collect stats for extensibility
    mags = {}

    #################
    # read lc files #
    #################
    # Using fits lc
    if lcext == '.fits':
        hdulist = fits.open(lcfile)

        if hdulist[0].header['DTR_EPD']:
            mags['rm'] = []
            for i in range(num_aps):
                mags['rm'].append(hdulist[1].data['IRM%d' % (i+1)])
            if epdlcrequired:
                mags['ep'] = []
                for i in range(num_aps):
                    mags['ep'].append(hdulist[1].data['EP%d' % (i+1)])

        elif not hdulist[0].header['DTR_EPD'] and not epdlcrequired:
            # a hack. some code has been written to rely on "EPD"
            # statistics. however, if you're skipping EPD, you want to
            # instead rely on IRM statistics. populating the EPD statistics
            mags['rm'] = []
            mags['ep'] = []
            for i in range(num_aps):
                mags['rm'].append(hdulist[1].data['IRM%d' % (i+1)])
                mags['ep'].append(hdulist[1].data['IRM%d' % (i+1)])
        else:
            print('expected DTR_EPD to be true in compute_lc_statistics')
            raise AssertionError

        if hdulist[0].header['DTR_TFA']:
            mags['tf'] = []
            for i in range(num_aps):
                mags['tf'].append(hdulist[1].data['TFA%d' % (i+1)])
        elif not hdulist[0].header['DTR_TFA'] and tfalcrequired:
            print(
                '{:s}Z: no TFA for {:s} and TFA is required, skipping...'.
                format(datetime.utcnow().isoformat(), lcfile))
            return None
    # Using text lightcurve files
    else:
        # Check for tfalc
        if tfalcrequired and lcext != '.tfalc':
            tfalcfile = lcfile.replace(lcext, '.tfalc')
            if not os.path.exists(tfalcfile):
                print('%sZ: no TFA mags available for %s and '
                      'TFALC is required, skipping...' %
                      (datetime.utcnow().isoformat(), lcfile))
                raise FileNotFoundError
            lcfile = tfalcfile

        if epdlcrequired and (lcext != '.epdlc' or '.tfalc'):
            epdlcfile = lcfile.replace(lcext, '.epdlc')
            if not os.path.exists(epdlcfile):
                print('%sZ: no EPD mags available for %s and '
                      'EPDLC is required, skipping ...' %
                      (datetime.utcnow().isoformat(), lcfile))
                raise FileNotFoundError
            lcfile = epdlcfile

        # Check that file is not empty
        if os.stat(lcfile).st_size == 0:
            print('%sZ: no mags available for %s, skipping...' %
                  (datetime.utcnow().isoformat(), lcfile))
            return None

        # RM
        mags['rm'] = np.genfromtxt(lcfile, usecols=rmcols, unpack=True)
        keys = ['rm']
        # EP
        if epdlcrequired:
            try:
                mags['ep'] = np.genfromtxt(lcfile, usecols=epcols, unpack=True)
                keys.append('ep')
            except:
                print('%sZ: no EPD mags available for %s!' %
                      (datetime.utcnow().isoformat(), lcfile))
        # TF
        if tfalcrequired:
            try:
                mags['tf'] = np.genfromtxt(lcfile, usecols=tfcols, unpack=True)
                keys.append('tf')
            except:
                print('%sZ: no TFA mags available for %s!' %
                      (datetime.utcnow().isoformat(), lcfile))
        # RF
        if rfcols:
            mags['rf'] = np.genfromtxt(lcfile, usecols=rfcols, unpack=True)
            keys.append('rf')


    ##################################
    # get statistics for each column #
    ##################################
    # Create results dict
    results = {'lcfile': lcfile,
               'lcobj': os.path.splitext(os.path.basename(lcfile))[0]}

    # Now, get statistics for each column
    statkeys = ['median', 'mad', 'mean', 'stdev', 'ndet']
    cols = [(x, ap) for x in keys for ap in range(num_aps)]

    for (k, ap) in cols:
        suf = '_%s%d' % (k, ap+1)
        if k not in mags:
            for s in statkeys:
                results[s+suf] = np.nan
                results[s+'_sigclip'+suf] = np.nan
            continue

        mag = mags[k][ap]

        if len(mag) > 4:
            finiteind = np.isfinite(mag)
            mag = mag[finiteind]
            med = np.median(mag)
            results['median'+suf] = med
            results['mad'+suf] = np.median(np.fabs(mag - med))
            results['mean'+suf] = np.mean(mag)
            results['stdev'+suf] = np.std(mag)
            results['ndet'+suf] = len(mag)

            if sigclip:
                sigclip_mag, lo, hi = stats_sigmaclip(mag,
                                                      low=sigclip,
                                                      high=sigclip)
                med_sigclip = np.median(sigclip_mag)
                results['median'+'_sigclip'+suf] = med_sigclip
                results['mad'+'_sigclip'+suf] = np.median(np.fabs(sigclip_mag -
                                                                  med_sigclip))
                results['mean'+'_sigclip'+suf] = np.mean(sigclip_mag)
                results['stdev'+'_sigclip'+suf] = np.std(sigclip_mag)
                results['ndet'+'_sigclip'+suf] = len(sigclip_mag)
            else:
                for s in ['median', 'mad', 'mean', 'stdev', 'ndet']:
                    results[s+'_sigclip'+suf] = np.nan
        else:
            for s in ['median', 'mad', 'mean', 'stdev', 'ndet']:
                results[s+suf] = np.nan
                results[s+'_sigclip'+suf] = np.nan
    if verbose:
        print('%sZ: done with statistics for %s' %
              (datetime.utcnow().isoformat(), lcfile))

    return results


def compute_lc_statistics_fits(lcfile,
                               sigclip=4.0,
                               num_aps=3,
                               desired_stages=['IRM','PCA','TFA'],
                               verbose=True,
                               istessandmaskedges=False):
    """
    Compute mean, median, MAD, stdev magnitudes for a given FITS lightcurve.

    Args:
        lcfile: A fits LC.

        num_aps (int): Number of apertures. Required for fits LC files.

        desired_stages (list): A list of detrending stages for which you want
        to assess statistics, e.g., ['IRM','PCA','TFA'], or ['EPD'], or
        whatever, depending on whatever stages of detrending you have done, and
        the keys you have used for them in the FITS data extension.

        istessandmaskedges (bool): If you're working with TESS data, dropping
        the edges near orbit start/end is standard.

    Returns:
        result (dict): Dictionary containing the following fields
            'lcfile': full path to lightcurve
            'lcobj': name of object
            Stat cols:
            Raw magnitudes:
                'median_rm[1-3]', 'mad_rm[1-3]', 'mean_rm[1-3]', 'stdev_rm[1-3]',
                'ndet_rm[1-3]',
                'median_sigclip_rm[1-3]', 'mad_sigclip_rm[1-3]', ...

            And similarly for EPD mags (_ep), TFA mags (_tf), raw fluxes (_rf),
            or PCA mags (_pm).
    """

    if istessandmaskedges:
        from tessutils import mask_orbit_start_and_end

    lcext = os.path.splitext(lcfile)[1]
    # Use dictionary to collect stats for extensibility
    mags = {}

    stage_to_key_dict = {
        'IRM':'rm',
        'TFA':'tf',
        'IFL':'rf',
        'PCA':'pm',
        'EPD':'ep'
    }

    ######################################
    # read lc files into mags dictionary #
    ######################################
    keys = []
    if lcext == '.fits':
        hdulist = fits.open(lcfile)

        for stage in desired_stages:

            if stage not in stage_to_key_dict.keys():
                raise ValueError('got detrending stage that is not recognized')

            mags[ stage_to_key_dict[stage] ] = []

            for i in range(num_aps):
                mags[ stage_to_key_dict[stage] ].append(
                    hdulist[1].data['{}{}'.format(stage, i+1)]
                )

            keys.append(stage_to_key_dict[stage])

    else:
        raise ValueError('expected a fits light curve')

    ##################################
    # get statistics for each column #
    ##################################
    # Create results dict
    results = {'lcfile': lcfile,
               'lcobj': os.path.splitext(os.path.basename(lcfile))[0]}

    # Now, get statistics for each column
    statkeys = ['median', 'mad', 'mean', 'stdev', 'ndet']
    cols = [(x, ap) for x in keys for ap in range(num_aps)]

    for (k, ap) in cols:
        suf = '_%s%d' % (k, ap+1)
        if k not in mags:
            for s in statkeys:
                results[s+suf] = np.nan
                results[s+'_sigclip'+suf] = np.nan
            continue

        mag = mags[k][ap]

        if istessandmaskedges:

            time = hdulist[1].data['TMID_BJD']

            orbitgap = 1
            expected_norbits = 2
            orbitpadding = 6/24
            raise_error = False

            _, mag = mask_orbit_start_and_end(time, mag, orbitgap=orbitgap,
                                              expected_norbits=expected_norbits,
                                              orbitpadding=orbitpadding,
                                              raise_error=raise_error)

        finiteind = np.isfinite(mag)
        mag = mag[finiteind]

        if len(mag) > 4:

            med = np.median(mag)
            results['median'+suf] = med
            results['mad'+suf] = np.median(np.fabs(mag - med))
            results['mean'+suf] = np.mean(mag)
            results['stdev'+suf] = np.std(mag)
            results['ndet'+suf] = len(mag)

            if sigclip:
                sigclip_mag, lo, hi = stats_sigmaclip(mag,
                                                      low=sigclip,
                                                      high=sigclip)
                med_sigclip = np.median(sigclip_mag)
                results['median'+'_sigclip'+suf] = med_sigclip
                results['mad'+'_sigclip'+suf] = np.median(np.fabs(sigclip_mag -
                                                                  med_sigclip))
                results['mean'+'_sigclip'+suf] = np.mean(sigclip_mag)
                results['stdev'+'_sigclip'+suf] = np.std(sigclip_mag)
                results['ndet'+'_sigclip'+suf] = len(sigclip_mag)
            else:
                for s in ['median', 'mad', 'mean', 'stdev', 'ndet']:
                    results[s+'_sigclip'+suf] = np.nan

        else:
            for s in ['median', 'mad', 'mean', 'stdev', 'ndet']:
                results[s+suf] = np.nan
                results[s+'_sigclip'+suf] = np.nan
    if verbose:
        print('%sZ: done with statistics for %s' %
              (datetime.utcnow().isoformat(), lcfile))

    return results


def parallel_compute_lc_statistics(lcdir, lcglob,
                                   fovcatalog, fovcathasgaiaids=False,
                                   tfalcrequired=False, num_aps=3,
                                   fovcatcols=(0,5), fovcatmaglabel='G',
                                   verbose=True, outfile=None,
                                   nworkers=16, workerntasks=500,
                                   rmcols=[14,19,24], epcols=[27,28,29],
                                   tfcols=[30,31,32], rfcols=[12,17,22],
                                   correctioncoeffs=None,
                                   sigclip=4.0, epdlcrequired=True,
                                   outheader=None):
    """
    This calculates statistics on all lc files in lcdir.

    Args:
        lcdir (str): directory containing lightcurves

        lcglob (str): glob to epd lcs, inside lcdir. E.g., '*.epdlc'. These
        contain the rlc, and are used to derive the filenames of the tfalcs.

        fovcatalog (str): path to the REFORMED fov catalog, which gets the
        catalog magnitude corresponding to canonical magnitude for any star.

        fovcathasgaiaids (bool): if the reformed FOV catalog has Gaia ids, set
        this to be true. The default is to assume HAT-IDs, which have different
        string lengths & and are read differently.

        num_aps (int): Number of apertures (default=3)

        epdlcrequired (bool): a variety of aperturephot.py tools assume that if
        you are creating statistics, you have run EPD. This isn't necessarily
        true (you may wish to get statistics on the instrumental raw
        magnitudes). If you set this to False, the statistics file will,
        hackily, populate the "EPD statistics" with IRM values.

    Output:

        Puts the results in text file outfile.
        outfile contains the following columns:

            object, catmag,
            median RM[1-3], MAD RM[1-3], mean RM[1-3], stdev RM[1-3], ndet RM[1-3],
            median EP[1-3], MAD EP[1-3], mean EP[1-3], stdev EP[1-3], ndet EP[1-3],
            median TF[1-3], MAD TF[1-3], mean TF[1-3], stdev TF[1-3], ndet TF[1-3],
            median RF[1-3], MAD RF[1-3], mean RF[1-3], stdev RF[1-3], ndet RF[1-3],
            corr_catmag AP[1-3]

        if a value is missing, it will be np.nan.

    Notes:
        For ISM, consider using correctioncoeffs as well. These are c1, c2
        resulting from a fit to the catalogmag-flux relation using the
        expression:

        catrmag = -2.5 * log10(flux/c1) + c2

        where the fit is done in the bright limit (8.0 < r < 12.0). this
        corrects for too-faint catalog mags because of crowding and blending.

        correctioncoeffs is like:
            [[ap1_c1,ap1_c2],[ap2_c1,ap2_c2],[ap3_c1,ap3_c2]]

    """
    lcfiles = glob(os.path.join(lcdir, lcglob))
    print('%sZ: %s %s lightcurves to process.' %
          (datetime.utcnow().isoformat(), len(lcfiles), lcglob))

    lc_statistics_worker = partial(compute_lc_statistics, rmcols=rmcols,
                                   epcols=epcols, tfcols=tfcols, rfcols=rfcols,
                                   sigclip=sigclip, num_aps=num_aps,
                                   epdlcrequired=epdlcrequired,
                                   tfalcrequired=tfalcrequired,
                                   verbose=verbose)
    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    results = pool.map(lc_statistics_worker, lcfiles)
    pool.close()
    pool.join()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    if not outfile:
        outfile = os.path.join(lcdir, 'lightcurve-statistics.txt')

    outf = open(outfile, 'wb')

    # Write header for stats.
    if outheader is not None:
        outheader += '\n'
        outf.write(outheader.encode('utf-8'))
    default_header = '# total objects: %s, sigmaclip used: %s\n' % (
        len(lcfiles), sigclip)
    outf.write(default_header.encode('utf-8'))

    apkeys=  ['rm']
    if epdlcrequired:
        apkeys.append('ep')
    if tfalcrequired:
        apkeys.append('tf')
    if rfcols:
        apkeys.append('rf')
    # tuples of column names and format string
    statkeys = [('median', '%.6f'), ('mad', '%.6f'), ('mean', '%.6f'),
                ('stdev', '%.6f'), ('ndet', '%s')]
    stat_cols = [{'name':'object', 'format':'%s', 'key':'lcobj'},
                 {'name':'catalog mag %s' % fovcatmaglabel, 'format':'%.3f',
                  'key':'catmag'}]
    outcolumnkey = ('# columns are:\n'
                    '# 0,1: object, catalog mag %s\n') % fovcatmaglabel
    # Generate list of columns
    for (ap, apnum) in itertools.product(apkeys, range(1, num_aps+1)):
        outcolkey = ('# ' + ('%d,'*(len(statkeys)-1)) + '%d: ')
        outcolkey %= tuple(range(len(stat_cols),
                                 len(stat_cols) + len(statkeys)))

        for stat in statkeys:
            stat_cols.append({'name':stat[0], 'format': stat[1],
                              'key':'%s_%s%s' % (stat[0], ap, apnum)})
            outcolkey = outcolkey + ('%s %s%d, ' % (stat[0], ap.upper(), apnum))
        outcolumnkey = outcolumnkey + outcolkey[:-2] + '\n'

        outcolkey = ('# ' + ('%d,'*(len(statkeys)-1)) + '%d: sigma-clipped ')
        outcolkey %= tuple(range(len(stat_cols),
                                 len(stat_cols) + len(statkeys)))
        for stat in statkeys:
            stat_cols.append({'name':stat[0], 'format': stat[1],
                              'key':'%s_sigclip_%s%s' % (stat[0], ap, apnum)})
            outcolkey = outcolkey + ('%s %s%d, ' % (stat[0], ap.upper(), apnum))
        outcolumnkey = outcolumnkey + outcolkey[:-2] + '\n'

    # Columns for corrected cat mag
    outcolkey = ('# ' + ('%d,'*(num_aps-1)) + '%d: ')
    outcolkey %= tuple(range(len(stat_cols),
                             len(stat_cols) + num_aps))
    for apnum in range(1,num_aps+1):
        stat_cols.append({'name':'corrmag', 'format':'%.3f',
                          'key':'corrmag_ap%s' % apnum})
        outcolkey += ('corrected cat mag AP%d, ' % apnum)
    outcolumnkey = outcolumnkey + outcolkey[:-2] + '\n'

    # Write header column to file
    outf.write(outcolumnkey.encode('utf-8'))

    # Generate format string
    formatstr = ''
    for col in stat_cols:
        formatstr = formatstr + col['format'] + ' '
    formatstr = formatstr[:-1] + '\n'

    # open the fovcatalog and read in the column magnitudes and hatids
    if not fovcathasgaiaids:
        # assume HAT-IDs, HAT-123-4567890, 17 character strings
        fovcat = np.genfromtxt(fovcatalog,
                               usecols=fovcatcols,
                               dtype='U17,f8',
                               names=['objid','mag'])
    else:
        # assume GAIA-IDs. From gaia2read, with "GAIA" id option, this is just
        # 19 character integers.
        fovcat = np.genfromtxt(fovcatalog,
                               usecols=fovcatcols,
                               dtype='U19,f8',
                               names=['objid','mag'])
    # Using a dictionary leads to ~ 300x speedup
    fovdict = dict(fovcat)

    # Cross-match objects with catalog to get catmags
    for stat in results:
        if stat is None:
            continue
        # find catalog mag for this object
        if stat['lcobj'] in fovdict:
            stat['catmag'] = fovdict[stat['lcobj']]
        elif 'median_tf3' in stat and not np.isnan(stat['median_tf3']):
            print('no catalog mag for %s, using median TF3 mag' % stat['lcobj'])
            stat['catmag'] = stat['median_tf3']
        elif 'median_ep3' in stat and not np.isnan(stat['median_ep3']):
            print('no catalog or median_tf3 mag for %s, using median EP3 mag' %
                  stat['lcobj'])
            stat['catmag'] = stat['median_ep3']
        else:
            print('WRN! no catalog, TF3 or EP3 mag for {:s}. using nan'.
                  format(stat['lcobj']))
            stat['catmag'] = np.nan

    # Compute correction coefficients
    if correctioncoeffs is True:
        stat_df = pd.DataFrame(results)
        fluxcols = ''.join('median_rf%d ' % (ap+1) for ap in range(num_aps))
        correctioncoeffs = compute_correction_coeffs(stat_df['catmag'],
                                                     stat_df[fluxcols.split()])

    for stat in results:
        if stat is None:
            continue
        # calculate corrected mags if present
        if (correctioncoeffs and len(correctioncoeffs) == 3 and
            rfcols and len(rfcols) == 3):

            ap1_c1, ap1_c2 = correctioncoeffs[0]
            ap2_c1, ap2_c2 = correctioncoeffs[1]
            ap3_c1, ap3_c2 = correctioncoeffs[2]

            stat['corrmag_ap1'] =-2.5*np.log10(stat['median_rf1']/ap1_c1)+ap1_c2
            stat['corrmag_ap2'] =-2.5*np.log10(stat['median_rf2']/ap2_c1)+ap2_c2
            stat['corrmag_ap3'] =-2.5*np.log10(stat['median_rf3']/ap3_c1)+ap3_c2

        else:
            stat['corrmag_ap1'] = stat['catmag']
            stat['corrmag_ap2'] = stat['catmag']
            stat['corrmag_ap3'] = stat['catmag']

        statvals = []
        for col in stat_cols:
            statvals.append(stat[col['key']])
        outline = formatstr % tuple(statvals)

        outf.write(outline.encode('utf-8'))

    outf.close()

    print('%sZ: wrote statistics to file %s' %
          (datetime.utcnow().isoformat(), outfile))

    return results


######################
# PLOTTING FUNCTIONS #
######################
def acf_percentiles_stats_and_plots(statdir, outprefix, make_plot=True,
                                    percentiles=[2,25,50,75,98], skipepd=False):
    """
    make csv files and (optionally) plots of ACF values at various time lags,
    evaluated at percentiles.
    """

    # get and read ACF files
    acfdir = os.path.join(statdir,'acf_stats')
    acfstatfiles = glob(os.path.join(acfdir,'*_acf_stats.csv'))

    df = read_acf_stat_files(acfstatfiles)

    # ACF vs time-lag statistics, summarized as percentiles at 2%, 25%, 50%,
    # 75%, and 98% percentiles.
    apstrlist = (['TFA1','TFA2','TFA3','RAW1','RAW2','RAW3']
                 if skipepd else
                 ['TFA1','TFA2','TFA3', 'EPD1','EPD2','EPD3',
                  'RAW1','RAW2','RAW3'])

    for apstr in apstrlist:
        try:
            csvname = ( os.path.join(
                statdir,'acf_percentiles_stats_{:s}.csv'.format(apstr.upper())
            ))
            if os.path.exists(csvname):
                continue

            timelags = np.sort(np.unique(df['LAG_TIME_HR']))

            percentile_dict = {}
            for timelag in timelags:

                percentile_dict[timelag] = {}

                sel = df['LAG_TIME_HR']==timelag

                for percentile in percentiles:
                    val = np.nanpercentile(df[sel][apstr+'_ACF'], percentile)
                    percentile_dict[timelag][percentile] = np.round(val,7)

            pctile_df = pd.DataFrame(percentile_dict)

            if make_plot:

                plt.close('all')
                fig, ax = plt.subplots(figsize=(4,3))

                markers = itertools.cycle(('o', 'v', '>', 'D', 's', 'P'))

                for ix, row in pctile_df.iterrows():
                    pctile = row.name
                    label = '{}%'.format(str(pctile))

                    timelags = nparr(row.index)
                    vals = nparr(row)

                    ax.plot(timelags, vals, label=label, marker=next(markers))

                ax.legend(loc='best', fontsize='xx-small')

                ax.set_yscale('linear')
                ax.set_xscale('log')
                ax.set_xlabel('ACF time lag [hr]')
                ax.set_ylabel('{:s} ACF value'.format(apstr.upper()))

                titlestr = '{:s} - {:d} ACFs - {:s}'.format(
                    outprefix,
                    len(acfstatfiles),
                    '{:s} percentiles'.format(repr(percentiles))
                )
                ax.set_title(titlestr, fontsize='small')

                plt.gca().grid(color='#a9a9a9',
                               alpha=0.9,
                               zorder=0,
                               linewidth=1.0,
                               linestyle=':')

                ax.set_ylim((-1,1))

                savname = ( os.path.join(
                    statdir,'acf_percentiles_stats_{:s}.png'.format(apstr.upper())
                ))
                fig.tight_layout()
                fig.savefig(savname, dpi=250)
                print('%sZ: made %s plot: %s' %
                      (datetime.utcnow().isoformat(), titlestr, savname))

            outdf = pctile_df.T
            outdf.index.name='lag_time_hr'
            outdf.to_csv(csvname)
            print('%sZ: wrote %s' %
                  (datetime.utcnow().isoformat(), csvname))

        except Exception as e:
            print('%sZ: failed to make percentiles for %s, err was %s' %
                  (datetime.utcnow().isoformat(), apstr, e))


def percentiles_RMSorMAD_stats_and_plots(statdir, outprefix, binned=False,
                                         make_percentiles_plot=True,
                                         percentiles_xlim=[4,17],
                                         percentiles_ylim=[1e-5,1e-1],
                                         percentiles=[2,25,50,75,98],
                                         yaxisval='RMS'):
    """
    make csv files and (optionally) plots of RMS (or MAD) vs magnitude
    statistics, evaluated at percentiles.
    """

    # get and read the most complete stats file (the TFA one!)
    tfastatfile = glob(os.path.join(statdir,'*.tfastats'))

    if not tfastatfile:
        raise AssertionError('cannot assess run if there is no tfa statfile')
    if len(tfastatfile)==1:
        tfastatfile = tfastatfile[0]
    else:
        raise AssertionError('something wrong with tfa statfile')

    if binned:
        stats = ap.read_binnedlc_stats_file(tfastatfile)
    else:
        stats = ap.read_stats_file(tfastatfile)

    if yaxisval=='RMS':
        yaxisstr='stdev_'
    elif yaxisval=='MAD':
        yaxisstr='mad_'
    else:
        raise ValueError('yaxisval must be RMS or MAD')

    # RMS vs magnitude statistics, at 1-mag differences, summarized as percentiles
    # at 2%, 25%, 50%, 75%, and 98% percentiles.
    for apstr in ['tf1','tf2','tf3',
                  'rm1','rm2','rm3',
                  'ep1','ep2','ep3']:

        try:
            medstr = 'med_'+apstr
            yvalstr = yaxisstr+apstr

            minmag = np.floor(np.nanmin(stats[medstr])).astype(int)
            maxmag = np.ceil(np.nanmax(stats[medstr])).astype(int)
            magdiff = 0.5
            mag_bins = [
                (me, me+magdiff) for me in np.arange(minmag, maxmag, magdiff)
            ]

            percentile_dict = {}
            for mag_bin in mag_bins:

                thismagmean = np.round(np.mean(mag_bin),2)
                percentile_dict[thismagmean] = {}

                thismin, thismax = min(mag_bin), max(mag_bin)
                sel = (stats[medstr] > thismin) & (stats[medstr] <= thismax)

                for percentile in percentiles:
                    val = np.nanpercentile(stats[sel][yvalstr], percentile)
                    percentile_dict[thismagmean][percentile] = np.round(val,7)

            pctile_df = pd.DataFrame(percentile_dict)

            if make_percentiles_plot:

                plt.close('all')
                fig, ax = plt.subplots(figsize=(4,3))

                # this plot is lines for each percentile in [2,25,50,75,98],
                # binned for every magnitude interval.
                markers = itertools.cycle(('o', 'v', '>', 'D', 's', 'P'))

                for ix, row in pctile_df.iterrows():
                    pctile = row.name
                    label = '{}%'.format(str(pctile))

                    midbins = nparr(row.index)
                    vals = nparr(row)

                    ax.plot(midbins, vals, label=label, marker=next(markers))

                # show the sullivan+2015 interpolated model
                if yaxisval=='RMS':
                    # overplot toy model I interpolated from that paper
                    Tmag = np.linspace(6, 16, num=200)
                    lnA = 3.29685004771
                    B = 0.8500214657
                    C = -0.2850416324
                    D = 0.039590832137
                    E = -0.00223080159
                    F = 4.73508403525e-5
                    ln_sigma_1hr = (
                        lnA + B*Tmag + C*Tmag**2 + D*Tmag**3 +
                        E*Tmag**4 + F*Tmag**5
                    )
                    sigma_1hr = np.exp(ln_sigma_1hr)
                    sigma_30min = sigma_1hr * np.sqrt(2)

                    ax.plot(Tmag, sigma_30min/1e6, 'k-', zorder=3, lw=2,
                            label='S+15 $\sigma_{\mathrm{30\,min}}$ (interp)')

                ax.legend(loc='best', fontsize='xx-small')

                ax.set_yscale('log')
                ax.set_xlabel('{:s} median instrument magnitude'.
                              format(apstr.upper()))
                ax.set_ylabel('{:s} {:s}'.
                              format(apstr.upper(), yaxisval))
                if percentiles_xlim:
                    ax.set_xlim(percentiles_xlim)
                if percentiles_ylim:
                    ax.set_ylim(percentiles_ylim)

                titlestr = '{:s} - {:d} LCs - {:s}'.format(
                    outprefix,
                    len(stats),
                    '{:s} percentiles'.format(repr(percentiles))
                )
                ax.set_title(titlestr, fontsize='small')

                plt.gca().grid(color='#a9a9a9',
                               alpha=0.9,
                               zorder=0,
                               linewidth=1.0,
                               linestyle=':')

                savname = ( os.path.join(
                    statdir,'percentiles_{:s}_vs_med_mag_{:s}.png'.
                    format(yaxisval, apstr.upper())
                ))
                fig.tight_layout()
                fig.savefig(savname, dpi=250)
                print('%sZ: made %s plot: %s' %
                      (datetime.utcnow().isoformat(), titlestr, savname))

            csvname = ( os.path.join(
                statdir,'percentiles_{:s}_vs_med_mag_{:s}.csv'.
                format(yaxisval, apstr.upper())
            ))
            pctile_df.to_csv(csvname, index=False)
            print('%sZ: wrote %s' %
                  (datetime.utcnow().isoformat(), csvname))
        except Exception as e:
            print('%sZ: failed to make percentiles for %s, err was %s' %
                  (datetime.utcnow().isoformat(), apstr, e))


def plot_raw_epd_tfa(time, rawmag, epdmag, tfamag, ap_index, savpath=None,
                     xlabel='BTJD = BJD - 2457000', skipepd=False):
    """
    Plot a 2 (or 3) row, 1 column plot with rows of:
        * raw mags vs time
        * EPD mags vs time (optional)
        * TFA mags vs time.

    args:
        time, rawmag, epdmag, tfamag (np.ndarray)

        ap_index (int): integer, e.g., "2" for aperture #2.

        skipepd (bool): if true, skips the EPD row.
    """

    plt.close('all')
    nrows = 2 if skipepd else 3
    fig, axs = plt.subplots(nrows=nrows, ncols=1, sharex=True, figsize=(6,4))

    axs = axs.flatten()

    apstr = 'AP{:d}'.format(ap_index)
    stagestrs = (['RM{:d}'.format(ap_index), 'TF{:d}'.format(ap_index)]
                 if skipepd else
                 ['RM{:d}'.format(ap_index), 'EP{:d}'.format(ap_index),
                  'TF{:d}'.format(ap_index)]
                )
    mags = [rawmag,tfamag] if skipepd else [rawmag,epdmag,tfamag]

    for ax, mag, txt in zip(axs, mags, stagestrs):

        ax.scatter(time, mag, c='black', alpha=0.9, zorder=2, s=3,
                   rasterized=True, linewidths=0)

        txt_x, txt_y = 0.99, 0.02
        t = ax.text(txt_x, txt_y, txt, horizontalalignment='right',
                verticalalignment='bottom', fontsize='small', zorder=3,
                transform=ax.transAxes)

        txt_x, txt_y = 0.99, 0.98
        stdmmag = np.nanstd(mag)*1e3
        if stdmmag > 0.1:
            stattxt = '$\sigma$ = {:.1f} mmag'.format(stdmmag)
            ndigits = 2
        elif stdmmag > 0.01:
            stattxt = '$\sigma$ = {:.2f} mmag'.format(stdmmag)
            ndigits = 3
        else:
            stattxt = '$\sigma$ = {:.3f} mmag'.format(stdmmag)
            ndigits = 4
        _ = ax.text(txt_x, txt_y, stattxt, horizontalalignment='right',
                verticalalignment='top', fontsize='small', zorder=3,
                transform=ax.transAxes)
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    for ax in axs:
        ylim = ax.get_ylim()
        ax.set_ylim((max(ylim), min(ylim)))

    axs[-1].set_xlabel(xlabel, fontsize='small')

    # make the y label
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('Magnitude', fontsize='small', labelpad=5)

    if not savpath:
        savpath = 'temp_{:s}.png'.format(apstr)

    fig.tight_layout(h_pad=-0.5)
    fig.savefig(savpath, dpi=250, bbox_inches='tight')
    print('%sZ: made plot: %s' % (datetime.utcnow().isoformat(), savpath))


def plot_tfa_templates(tfa_templates_path, savdir):

    df = pd.read_csv(tfa_templates_path, sep=' ', header=None,
                     names=['fname','x','y'])

    for fname in df['fname']:
        plot_raw_tfa_bkgd_fits(fname, savdir)


def plot_raw_tfa_bkgd_fits(fitspath, savdir, maskbadtimes=True, ap_index=2,
                           customstr=''):

    # define savpath
    savpath = os.path.join(
        savdir,
        os.path.basename(fitspath).replace('.fits','_raw_tfa_bkgd.png')
    )
    if os.path.exists(savpath):
        print('found {}, return'.format(savpath))
        return

    # load lc
    hdulist = fits.open(fitspath)
    hdr = hdulist[0].header
    time_bjd = hdulist[1].data['TMID_BJD']
    barycorr = hdulist[1].data['BJDCORR']
    tfa_mag = hdulist[1].data['TFA{}'.format(ap_index)]
    raw_mag = hdulist[1].data['IRM{}'.format(ap_index)]
    bkgd = hdulist[1].data['BGV']
    hdulist.close()

    if maskbadtimes:
        from tessutils import badtimewindows
        # BJD = TJD + LTT_corr
        tjd = time_bjd - barycorr
        toffset = 2457000
        tjd -= toffset

        good_inds = np.ones_like(tjd).astype(bool)
        for window in badtimewindows:
            bad_window_inds = (tjd > min(window)) & (tjd < max(window))
            good_windows_inds = ~bad_window_inds
            good_inds &= good_windows_inds

        time = time_bjd[good_inds]
        tfa_mag = tfa_mag[good_inds]
        raw_mag = raw_mag[good_inds]
        bkgd = bkgd[good_inds]

    else:
        time = time_bjd

    plot_raw_tfa_bkgd(time, raw_mag, tfa_mag, bkgd, ap_index, savpath,
                      customstr=customstr)



def plot_raw_tfa_bkgd(time, rawmag, tfamag, bkgdval, ap_index, savpath=None,
                      xlabel='BJDTDB', customstr='', tfatime=None):
    """
    Plot a 2 (or 3) row, 1 column plot with rows of:
        * raw mags vs time
        * EPD mags vs time (optional)
        * TFA mags vs time.

    args:
        time, rawmag, epdmag, tfamag (np.ndarray)

        ap_index (int): integer, e.g., "2" for aperture #2.

        skipepd (bool): if true, skips the EPD row.

    kwargs:
        customstr: string that goes on top right of plot

        tfatime: if passed, "time" is used to plot rawmag and bkgdval,
        "tfatime" is used to plot tfamag. Otherwise "time" is used for all of
        them.
    """

    plt.close('all')
    nrows = 3
    fig, axs = plt.subplots(nrows=nrows, ncols=1, sharex=True, figsize=(6,4))

    axs = axs.flatten()

    apstr = 'AP{:d}'.format(ap_index)
    stagestrs = ( ['RM{:d}'.format(ap_index),
                   'TF{:d}'.format(ap_index),
                   'BKGD{:d}'.format(ap_index)] )

    yvals = [rawmag,tfamag,bkgdval]
    nums = list(range(len(yvals)))

    for ax, yval, txt, num in zip(axs, yvals, stagestrs, nums):

        if isinstance(tfatime, np.ndarray) and 'TF' in txt:
            ax.scatter(tfatime, yval, c='black', alpha=0.9, zorder=2, s=3,
                       rasterized=True, linewidths=0)
        else:
            ax.scatter(time, yval, c='black', alpha=0.9, zorder=2, s=3,
                       rasterized=True, linewidths=0)

        #txt_x, txt_y = 0.99, 0.02
        #t = ax.text(txt_x, txt_y, txt, horizontalalignment='right',
        #            verticalalignment='bottom', fontsize='small', zorder=3,
        #            transform=ax.transAxes)

        if num in [0]:
            txt_x, txt_y = 0.99, 0.98
            mag = yval
            stdmmag = np.nanstd(mag)*1e3
            if stdmmag > 0.1:
                stattxt = '$\sigma$ = {:.1f} mmag{}'.format(stdmmag, customstr)
                ndigits = 2
            elif stdmmag > 0.01:
                stattxt = '$\sigma$ = {:.2f} mmag{}'.format(stdmmag, customstr)
                ndigits = 3
            else:
                stattxt = '$\sigma$ = {:.3f} mmag{}'.format(stdmmag, customstr)
                ndigits = 4
            _ = ax.text(txt_x, txt_y, stattxt, horizontalalignment='right',
                    verticalalignment='top', fontsize='small', zorder=3,
                    transform=ax.transAxes)
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    for ax in axs:
        ylim = ax.get_ylim()
        ax.set_ylim((max(ylim), min(ylim)))

    axs[-1].set_xlabel(xlabel, fontsize='small')

    # make the y label
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('IRM / TFA / bkgdval', fontsize='small', labelpad=5)

    if not savpath:
        savpath = 'temp_{:s}.png'.format(apstr)

    fig.tight_layout(h_pad=-0.3)
    fig.savefig(savpath, dpi=250, bbox_inches='tight')
    print('%sZ: made plot: %s' % (datetime.utcnow().isoformat(), savpath))





def plot_lightcurve_and_ACF(
    rawlags, rawacf, epdlags, epdacf, tfalags, tfaacf,
    rawtime, rawflux, epdtime, epdflux, tfatime, tfaflux,
    ap, savpath=None,
    xlabeltime='BJDTDB',xlabelacf='ACF lag [days]',
    skipepd=False):
    """
    Plot a 3 row, 2 column plot with rows of:
        * raw mags vs time (and ACF)
        * EPD mags vs time (and ACF)
        * TFA mags vs time. (and ACF)

    Args:
        lags and ACF: computed with e.g., `compute_acf_statistics_worker`, or
        `autocorr_magseries`.

        ap (int): integer identifying aperture size
    """

    from matplotlib.ticker import FormatStrFormatter

    plt.close('all')
    nrows = 2 if skipepd else 3
    fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(7,4),
                            gridspec_kw= {'width_ratios':[4,1],
                                          'height_ratios':[1]*nrows})
    if skipepd:
        a0,a1,a2,a3 = axs.flatten()
    else:
        a0,a1,a2,a3,a4,a5 = axs.flatten()

    apstr = 'AP{:d}'.format(ap)
    stagestrs = (['RM{:d}'.format(ap),
                  'TF{:d}'.format(ap)]
                 if skipepd else
                 ['RM{:d}'.format(ap),
                  'EP{:d}'.format(ap),
                  'TF{:d}'.format(ap)]
                )

    axlist = (a0,a2) if skipepd else (a0,a2,a4)
    fluxlist = [rawflux,tfaflux] if skipepd else [rawflux,epdflux,tfaflux]
    timelist = [rawtime,tfatime] if skipepd else [rawtime,epdtime,tfatime]

    for ax, mag, time, txt in zip(axlist, fluxlist, timelist, stagestrs):

        ax.scatter(time, mag, c='black', alpha=0.9, zorder=2, s=3,
                   rasterized=True, linewidths=0)

        ax.set_ylim(
            (np.mean(mag)  - 4*np.std(mag),
             np.mean(mag)  + 4*np.std(mag))
        )

        txt_x, txt_y = 0.99, 0.98
        t = ax.text(txt_x, txt_y, txt, horizontalalignment='right',
                verticalalignment='top', fontsize='small', zorder=3,
                transform=ax.transAxes)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    axlist = (a1,a3) if skipepd else (a1,a3,a5)
    acflist = [rawflux,tfaflux] if skipepd else [rawacf, epdacf, tfaacf]
    laglist = [rawlags,tfalags] if skipepd else [rawlags, epdlags, tfalags]

    for ax, acf, lag, txt in zip(axlist, acflist, laglist, stagestrs):

        ax.plot(lag, acf, c='k', linestyle='-', marker='o',
                markerfacecolor='k', markeredgecolor='k', ms=1, lw=0.5,
                zorder=1, rasterized=True)

        ax.set_xscale('log')
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    for ax in (a0,a1):
        ax.set_xticklabels([])

    fig.text(0.45,0.05, xlabeltime, ha='center')
    fig.text(0.88,0.05, xlabelacf, ha='center')

    # make the y label
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('instrument mag (left), ACF (right)',
                         fontsize='small', labelpad=5)

    if not savpath:
        savpath = 'temp_{:s}.png'.format(apstr)

    fig.tight_layout(h_pad=0., w_pad=0)
    fig.savefig(savpath, dpi=250, bbox_inches='tight')
    print('%sZ: made plot: %s' % (datetime.utcnow().isoformat(), savpath))

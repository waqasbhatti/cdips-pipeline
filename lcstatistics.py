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

    compute_lc_statistics:
        worker to compute median, MAD, mean, stdev for each lightcurve

    parallel_compute_lc_statistics:
        compute above statistics for a directory of lightcurves and compile

Plotting functions
    percentiles_RMSorMAD_stats_and_plots
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

TODO:
move items from aperturephot.py here, and update all respective calls.

"""

import os, pickle, itertools
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import aperturephot as ap
from glob import glob
from astropy import units as units, constants as constants
from datetime import datetime

from numpy import array as nparr
from scipy.stats import sigmaclip as stats_sigmaclip

from astrobase.periodbase import macf
from astrobase.varbase.autocorr import autocorr_magseries

from datetime import datetime
import multiprocessing as mp

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

    # open the statfile and read all the columns
    stats = np.genfromtxt(
        statsfile,
        dtype=(
            'U{:d},f8,'
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM1
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM2
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RM3
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP1
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP2
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # EP3
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF1
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF2
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # TF3
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF1
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF2
            'f8,f8,f8,f8,i8,f8,f8,f8,f8,i8,'  # RF3
            'f8,f8,f8'.format(idstrlength)    # corrmags
        ),
        names=[
            'lcobj','catmag',
            'median_rm1','mad_rm1','mean_rm1','stdev_rm1','ndet_rm1',
            'median_sigclip_rm1','mad_sigclip_rm1','mean_sigclip_rm1',
            'stdev_sigclip_rm1','ndet_sigclip_rm1',
            'median_rm2','mad_rm2','mean_rm2','stdev_rm2','ndet_rm2',
            'median_sigclip_rm2','mad_sigclip_rm2','mean_sigclip_rm2',
            'stdev_sigclip_rm2','ndet_sigclip_rm2',
            'median_rm3','mad_rm3','mean_rm3','stdev_rm3','ndet_rm3',
            'median_sigclip_rm3','mad_sigclip_rm3','mean_sigclip_rm3',
            'stdev_sigclip_rm3','ndet_sigclip_rm3',
            'median_ep1','mad_ep1','mean_ep1','stdev_ep1','ndet_ep1',
            'median_sigclip_ep1','mad_sigclip_ep1','mean_sigclip_ep1',
            'stdev_sigclip_ep1','ndet_sigclip_ep1',
            'median_ep2','mad_ep2','mean_ep2','stdev_ep2','ndet_ep2',
            'median_sigclip_ep2','mad_sigclip_ep2','mean_sigclip_ep2',
            'stdev_sigclip_ep2','ndet_sigclip_ep2',
            'median_ep3','mad_ep3','mean_ep3','stdev_ep3','ndet_ep3',
            'median_sigclip_ep3','mad_sigclip_ep3','mean_sigclip_ep3',
            'stdev_sigclip_ep3','ndet_sigclip_ep3',
            'median_tf1','mad_tf1','mean_tf1','stdev_tf1','ndet_tf1',
            'median_sigclip_tf1','mad_sigclip_tf1','mean_sigclip_tf1',
            'stdev_sigclip_tf1','ndet_sigclip_tf1',
            'median_tf2','mad_tf2','mean_tf2','stdev_tf2','ndet_tf2',
            'median_sigclip_tf2','mad_sigclip_tf2','mean_sigclip_tf2',
            'stdev_sigclip_tf2','ndet_sigclip_tf2',
            'median_tf3','mad_tf3','mean_tf3','stdev_tf3','ndet_tf3',
            'median_sigclip_tf3','mad_sigclip_tf3','mean_sigclip_tf3',
            'stdev_sigclip_tf3','ndet_sigclip_tf3',
            'median_rf1','mad_rf1','mean_rf1','stdev_rf1','ndet_rf1',
            'median_sigclip_rf1','mad_sigclip_rf1','mean_sigclip_rf1',
            'stdev_sigclip_rf1','ndet_sigclip_rf1',
            'median_rf2','mad_rf2','mean_rf2','stdev_rf2','ndet_rf2',
            'median_sigclip_rf2','mad_sigclip_rf2','mean_sigclip_rf2',
            'stdev_sigclip_rf2','ndet_sigclip_rf2',
            'median_rf3','mad_rf3','mean_rf3','stdev_rf3','ndet_rf3',
            'median_sigclip_rf3','mad_sigclip_rf3','mean_sigclip_rf3',
            'stdev_sigclip_rf3','ndet_sigclip_rf3',
            'corrmag_ap1','corrmag_ap2','corrmag_ap3',
        ]
    )

    return stats

#####################################
# FUNCTIONS TO CALCULATE STATISTICS #
#####################################

def compute_acf_statistics_worker(task, n_apertures=3, timename='TMID_BJD',
                                  filterwindow=7, istessffi=True,
                                  isfitslc=True):

    #NOTE : might want to check a couple smoothing values ("filterwindows")...

    try:
        if not istessffi:
            raise NotImplementedError(
            'this function assumes an ffi cadence of 30 minutes to evaluate ACFs. '
            'to generalize this you just need to interpolate, but I am lazy.'
            )

        tfafile, outdir, eval_times_hr, skipepd = task

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
        for ap in range(1,n_apertures+1):

            rawap = 'IRM{:d}'.format(ap)
            epdap = 'EP{:d}'.format(ap)
            tfaap = 'TFA{:d}'.format(ap)
            errap = 'IRE{:d}'.format(ap)

            time = lcdata[timename]
            flux_raw = lcdata[rawap]
            if not skipepd:
                flux_epd = lcdata[epdap]
            flux_tfa = lcdata[tfaap]
            err = lcdata[errap]

            acf_raw = autocorr_magseries(time, flux_raw, err, maxlags=None,
                                         fillgaps='noiselevel', sigclip=5,
                                         magsarefluxes=True,
                                         filterwindow=filterwindow)
            if not skipepd:
                acf_epd = autocorr_magseries(time, flux_epd, err, maxlags=None,
                                             fillgaps='noiselevel', sigclip=5,
                                             magsarefluxes=True,
                                             filterwindow=filterwindow)
            acf_tfa = autocorr_magseries(time, flux_tfa, err, maxlags=None,
                                         fillgaps='noiselevel', sigclip=5,
                                         magsarefluxes=True,
                                         filterwindow=filterwindow)

            if not np.isclose(acf_raw['cadence']*24*60, 30, atol=1e-3):
                raise NotImplementedError(
                'this function assumes an ffi cadence of 30 minutes to evaluate ACFs. '
                'to generalize this you just need to interpolate, but I am lazy.'
                )

            apstr = 'AP{}'.format(ap)
            d_pkl[apstr] = {}
            acfs = [acf_raw, acf_epd, acf_tfa] if not skipepd else [acf_raw, acf_tfa]
            dtrs = ['raw','epd','tfa'] if not skipepd else ['raw','tfa']
            for acf, dtr in zip(acfs, dtrs):

                d_pkl[apstr]['acf_{}'.format(dtr)] = acf['acf']
                d_pkl[apstr]['lag_time_{}'.format(dtr)] = acf['lags']*acf['cadence']
                d_pkl[apstr]['itimes_{}'.format(dtr)] = acf['itimes']
                d_pkl[apstr]['ifluxs_{}'.format(dtr)] = acf['imags']
                d_pkl[apstr]['ierrs_{}'.format(dtr)] = acf['ierrs']
                d_pkl[apstr]['cadence'] = acf['cadence']

            # assuming FFI cadence of 30 minutes, evalute ACF at desired lags.
            n_acf_vals = len(acf_tfa['acf'])

            # in "TUNE" mode, might want short-timescale ACFs (otherwise,
            # fails)
            if max(2*eval_times_hr)>n_acf_vals:
                eval_times_hr = np.array([1,2,6,12])
            else:
                pass

            eval_times_hr = nparr(eval_times_hr)
            outdf['LAG_TIME_HR'] = eval_times_hr

            colstrs = (
                ['RAW{:d}'.format(ap),'EPD{:d}'.format(ap),'TFA{:d}'.format(ap)]
                if not skipepd else
                ['RAW{:d}'.format(ap),'TFA{:d}'.format(ap)]
            )
            for acf, colstr in zip(acfs, colstrs):

                # wonky indexing scheme. 30 minute cadence -> e.g., 1hr lag is
                # index number 2.
                these_acf_vals = acf['acf'][2*eval_times_hr]

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
    nworkers=16, maxworkertasks=1000, skipepd=False):
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
    """

    print('%sZ: %s files to compute ACF statistics for' %
          (datetime.utcnow().isoformat(), len(tfafiles)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, outdir, eval_times_hr, skipepd) for x in tfafiles]

    # fire up the pool of workers
    results = pool.map(compute_acf_statistics_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {result for result in results}


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
    fluxes = pd.DataFrame(fluxes)[bright_limit]
    catmags = catmags[bright_limit]

    correctioncoeffs = []
    for ap_col in fluxes:
        try:
            c2 = np.polyfit(np.zeros_like(catmags),
                            catmags + 2.5*np.log10(fluxes[ap_col]), 0)[0]
            correctioncoeffs.append([1.0, c2])
        except:
            print('ERR! %sZ: could not compute correction coefficients' %
                  datetime.utcnow().isoformat())
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
                          fitslcnottxt=None):
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
    elif lcext == '.epdlc' or lcext == '.tfalc':
        # Check for tfalc
        if tfalcrequired and lcext == '.epdlc':
            tfalcfile = lcfile.replace('.epdlc', '.tfalc')
            if not os.path.exists(tfalcfile):
                print('%sZ: no TFA mags available for %s and '
                      'TFALC is required, skipping...' %
                      (datetime.utcnow().isoformat(), lcfile))
                raise FileNotFoundError
            lcfile = tfalcfile
        elif lcext == '.tfalc':
            tfalcrequired = True

        # Check that file is not empty
        if os.stat(lcfile).st_size == 0:
            print('%sZ: no mags available for %s, skipping...' %
                  (datetime.utcnow().isoformat(), lcfile))
            return None

        # RM and EP
        try:
            mags['rm'] = np.genfromtxt(lcfile, usecols=rmcols, unpack=True)
            mags['ep'] = np.genfromtxt(lcfile, usecols=epcols, unpack=True)
        except:
            if not epdlcrequired: 
                mags['ep'] = mags['rm']
            else:
                print('%sZ: no EPD mags available for %s!' %
                      (datetime.utcnow().isoformat(), lcfile))
        # TF
        if tfalcrequired:
            try:
                mags['tf'] = np.genfromtxt(lcfile, usecols=tfcols, unpack=True)
            except:
                print('%sZ: no TFA mags available for %s!' %
                      (datetime.utcnow().isoformat(), lcfile))
        # RF
        if rfcols:
            mags['rf'] = np.genfromtxt(lcfile, usecols=rfcols, unpack=True)

    else:
        raise NotImplementedError

    ##################################
    # get statistics for each column #
    ##################################

    # Create results dict
    results = {'lcfile': lcfile,
               'lcobj': os.path.splitext(os.path.basename(lcfile))[0]}

    # Now, get statistics for each column
    keys = ['rf', 'rm', 'ep', 'tf']
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
    print('%sZ: done with statistics for %s' %
          (datetime.utcnow().isoformat(), lcfile))
    return results


def lc_statistics_worker(task):
    """
    Wrapper for compute_lc_statistics function.
    """
    try:
        return compute_lc_statistics(task[0], **task[1])
    except Exception as e:
        print('SOMETHING WENT WRONG! task was %s' % task)
        return None


def parallel_compute_lc_statistics(lcdir, lcglob,
                                   fovcatalog, fovcathasgaiaids=False,
                                   tfalcrequired=False, num_aps=3,
                                   fovcatcols=(0,5), fovcatmaglabel='G',
                                   outfile=None,
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

    tasks = [[f, {'rmcols':rmcols,
                  'epcols':epcols,
                  'tfcols':tfcols,
                  'rfcols':rfcols,
                  'sigclip':sigclip,
                  'num_aps':num_aps,
                  'epdlcrequired':epdlcrequired,
                  'tfalcrequired':tfalcrequired}] for f in lcfiles]

    pool = mp.Pool(nworkers,maxtasksperchild=workerntasks)
    results = pool.map(lc_statistics_worker, tasks)
    pool.close()
    pool.join()

    print('%sZ: done. %s lightcurves processed.' %
          (datetime.utcnow().isoformat(), len(lcfiles)))

    if not outfile:
        outfile = os.path.join(lcdir, 'lightcurve-statistics.txt')

    outf = open(outfile, 'wb')

    # Write header for stats.
    default_header = '# total objects: %s, sigmaclip used: %s\n' % (
        len(lcfiles), sigclip)
    outf.write(default_header.encode('utf-8'))
    if outheader is not None:
        outf.write(outheader.encode('utf-8'))

    apkeys = ['rm', 'ep', 'tf', 'rf']
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
        elif not np.isnan(stat['median_tf3']):
            print('no catalog mag for %s, using median TF3 mag' % stat['lcobj'])
            stat['catmag'] = stat['median_tf3']
        elif not np.isnan(stat['median_ep3']):
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

    from matplotlib.ticker import FormatStrFormatter

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


def plot_lightcurve_and_ACF(
    rawlags, rawacf, epdlags, epdacf, tfalags, tfaacf,
    rawtime, rawflux, epdtime, epdflux, tfatime, tfaflux,
    ap, savpath=None,
    xlabeltime='BJD - 2457000',xlabelacf='ACF lag [days]',
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
            (np.mean(mag)  - 3*np.std(mag),
             np.mean(mag)  + 3*np.std(mag))
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
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    for ax in (a0,a2,a1,a3):
        ax.set_xticklabels([])

    fig.text(0.45,0.05, xlabeltime, ha='center')
    fig.text(0.88,0.05, xlabelacf, ha='center')

    # make the y label
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('instrument mag (left), ACF (right)',
                         fontsize='small', labelpad=3)

    if not savpath:
        savpath = 'temp_{:s}.png'.format(apstr)

    fig.tight_layout(h_pad=0., w_pad=0)
    fig.savefig(savpath, dpi=250, bbox_inches='tight')
    print('%sZ: made plot: %s' % (datetime.utcnow().isoformat(), savpath))




# -*- coding: utf-8 -*-
from __future__ import division, print_function
"""
functions for assessing statistics of lightcurves you've made.

percentiles_MAD_stats_and_plots:
    make csv files and (optionally) percentiles plots of MAD vs magnitude.

TODO:
move items from aperturephot.py here, and update all respective calls.

FIXME: break these out into their own module
1. RMS vs. MAG plot for EPD and TFA lightcurves
2. MAD vs. MAG plot for EPD and TFA lightcurves
3. ratios of RMS and MAD vs. MAG for CCD 6,7,8 to that of CCD 5
4. binned LC versions of these plots, using 10, 30, and 60 minute binning


"""

import os, pickle
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import aperturephot as ap
from glob import glob
from astropy import units as units, constants as constants
from datetime import datetime

from numpy import array as nparr

from astrobase.periodbase import macf
from astrobase.varbase.autocorr import autocorr_magseries

from datetime import datetime
import multiprocessing as mp


################################
# LIGHTCURVE READING FUNCTIONS #
################################

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

    30 ep1    TFA magnitude for aperture 1
    31 ep2    TFA magnitude for aperture 2
    32 ep3    TFA magnitude for aperture 3
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


#####################################
# FUNCTIONS TO CALCULATE STATISTICS #
#####################################

def compute_acf_statistics_worker(task, n_apertures=3, timename='btjd',
                                  filterwindow=7, istessffi=True):

    #NOTE : might want to check a couple smoothing values ("filterwindows")...

    if not istessffi:
        raise NotImplementedError(
        'this function assumes an ffi cadence of 30 minutes to evaluate ACFs. '
        'to generalize this you just need to interpolate, but I am lazy.'
        )

    tfafile, outdir, eval_times_hr = task

    lcdata = read_tfa_lc(tfafile)

    outpickle = os.path.join(
        outdir,
        os.path.basename(tfafile).replace('.tfalc','_acf_stats.pickle')
    )
    outcsv = os.path.join(
        outdir,
        os.path.basename(tfafile).replace('.tfalc','_acf_stats.csv')
    )

    d_pkl, outdf = {}, pd.DataFrame({})
    for ap in range(1,n_apertures+1):

        rawap = 'RM{:d}'.format(ap)
        epdap = 'EP{:d}'.format(ap)
        tfaap = 'TF{:d}'.format(ap)
        errap = 'RMERR{:d}'.format(ap)

        time = lcdata[timename]
        flux_raw = lcdata[rawap]
        flux_epd = lcdata[epdap]
        flux_tfa = lcdata[tfaap]
        err = lcdata[errap]

        acf_raw = autocorr_magseries(time, flux_raw, err, maxlags=None,
                                     fillgaps='noiselevel', sigclip=5,
                                     magsarefluxes=True,
                                     filterwindow=filterwindow)
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
        for acf, dtr in zip([acf_raw, acf_epd, acf_tfa],['raw','epd','tfa']):

            d_pkl[apstr]['acf_{}'.format(dtr)] = acf['acf']
            d_pkl[apstr]['lag_time_{}'.format(dtr)] = acf['lags']*acf['cadence']
            d_pkl[apstr]['itimes_{}'.format(dtr)] = acf['itimes']
            d_pkl[apstr]['ifluxs_{}'.format(dtr)] = acf['imags']
            d_pkl[apstr]['ierrs_{}'.format(dtr)] = acf['ierrs']
            d_pkl[apstr]['cadence'] = acf['cadence']

        # assuming FFI cadence of 30 minutes, evalute ACF at desired lags.
        eval_times_hr = nparr(eval_times_hr)
        outdf['LAG_TIME_HR'] = eval_times_hr
        for acf, colstr in zip(
            [acf_raw, acf_epd, acf_tfa],
            ['RAW{:d}'.format(ap),'EPD{:d}'.format(ap),'TFA{:d}'.format(ap)]
        ):
            # wonky indexing scheme. 30 minute cadence -> e.g., 1hr lag is
            # index number 2.
            these_acf_vals = acf['acf'][2*eval_times_hr]

            outdf[colstr+"_ACF"] = these_acf_vals

    with open(outpickle, 'w') as f:
        pickle.dump(d_pkl, f, pickle.HIGHEST_PROTOCOL)
    print('wrote {}'.format(outpickle))

    outdf.to_csv(outcsv, index=False)
    print('wrote {}'.format(outcsv))

    return 1


def parallel_compute_acf_statistics(
    tfafiles, outdir, eval_times_hr=[1,2,6,12,24,48,96,120,144,192],
    nworkers=16, maxworkertasks=1000
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
    """

    print('%sZ: %s files to compute ACF statistics for' %
          (datetime.utcnow().isoformat(), len(tfafiles)))

    pool = mp.Pool(nworkers,maxtasksperchild=maxworkertasks)

    tasks = [(x, outdir, eval_times_hr) for x in tfafiles]

    # fire up the pool of workers
    results = pool.map(compute_acf_statistics_worker, tasks)

    # wait for the processes to complete work
    pool.close()
    pool.join()

    return {result for result in results}



######################
# PLOTTING FUNCTIONS #
######################

def percentiles_MAD_stats_and_plots(statdir, outprefix, binned=False,
                                    make_percentiles_plot=True,
                                    percentiles_xlim=[4,17],
                                    percentiles_ylim=[1e-5,1e-1],
                                    percentiles=[2,25,50,75,98]):
    """
    make csv files and (optionally) plots of MAD vs magnitude statistics,
    evaluated at percentiles.
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

    # MAD vs magnitude statistics, at 1-mag differences, summarized as percentiles
    # at 2%, 25%, 50%, 75%, and 98% percentiles.
    for apstr in ['tf1','tf2','tf3',
                  'rm1','rm2','rm3',
                  'ep1','ep2','ep3']:

        try:
            medstr = 'med_'+apstr
            madstr = 'mad_'+apstr

            minmag = np.floor(np.nanmin(stats[medstr])).astype(int)
            maxmag = np.ceil(np.nanmax(stats[medstr])).astype(int)
            mag_bins = [
                (me, me+1) for me in np.arange(minmag, maxmag, 1)
            ]

            percentile_dict = {}
            for mag_bin in mag_bins:

                thismagmean = np.round(np.mean(mag_bin),1)
                percentile_dict[thismagmean] = {}

                thismin, thismax = min(mag_bin), max(mag_bin)
                sel = (stats[medstr] > thismin) & (stats[medstr] <= thismax)

                for percentile in percentiles:
                    val = np.nanpercentile(stats[sel][madstr], percentile)
                    percentile_dict[thismagmean][percentile] = np.round(val,7)

            pctile_df = pd.DataFrame(percentile_dict)

            if make_percentiles_plot:
                import itertools

                plt.close('all')
                fig, ax = plt.subplots(figsize=(4,3))

                # this was originally a bonafide percentiles plot. however getting
                # it to look right was surprisingly annoying.
                # this plot is therefore just lines for each percentile in
                # [2,25,50,75,98], binned for every magnitude interval.

                markers = itertools.cycle(('o', 'v', '>', 'D', 's', 'P'))

                for ix, row in pctile_df.iterrows():
                    pctile = row.name
                    label = '{}%'.format(str(pctile))

                    midbins = nparr(row.index)
                    vals = nparr(row)

                    ax.plot(midbins, vals, label=label, marker=next(markers))

                ax.legend(loc='best', fontsize='xx-small')

                ax.set_yscale('log')
                ax.set_xlabel('{:s} median instrument magnitude'.format(apstr.upper()))
                ax.set_ylabel('{:s} median abs. dev.'.format(apstr.upper()))
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
                    statdir,'percentiles_MAD_vs_med_mag_{:s}.png'.format(apstr.upper())
                ))
                fig.tight_layout()
                fig.savefig(savname, dpi=250)
                print('%sZ: made %s plot: %s' %
                      (datetime.utcnow().isoformat(), titlestr, savname))

            csvname = ( os.path.join(
                statdir,'percentiles_MAD_vs_med_mag_{:s}.csv'.format(apstr.upper())
            ))
            pctile_df.to_csv(csvname, index=False)
            print('%sZ: wrote %s' %
                  (datetime.utcnow().isoformat(), csvname))
        except Exception as e:
            print('%sZ: failed to make percentiles for %s, err was %s' %
                  (datetime.utcnow().isoformat(), apstr, e))


def plot_raw_epd_tfa(time, rawmag, epdmag, tfamag, ap_index, savpath=None,
                     xlabel='BTJD = BJD - 2457000'):
    """
    Plot a 3 row, 1 column plot with rows of:
        * raw mags vs time
        * EPD mags vs time
        * TFA mags vs time.

    args:
        time, rawmag, epdmag, tfamag (np.ndarray)

        ap_index (int): integer, e.g., "2" for aperture #2.
    """

    from matplotlib.ticker import FormatStrFormatter

    plt.close('all')
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6,4))

    axs = axs.flatten()

    apstr = 'AP{:d}'.format(ap_index)
    stagestrs = ['RM{:d}'.format(ap_index),
                 'EP{:d}'.format(ap_index),
                 'TF{:d}'.format(ap_index)]

    for ax, mag, txt in zip(axs, [rawmag,epdmag,tfamag], stagestrs):
        ax.plot(time, mag, c='k', linestyle='-', marker='o',
                markerfacecolor='k', markeredgecolor='k', ms=3, lw=0.,
                zorder=1, rasterized=True)
        ax.plot(time, mag, c='orange', linestyle='-', marker='o',
                markerfacecolor='orange', markeredgecolor='orange', ms=1.5,
                lw=0., zorder=2, rasterized=True)
        txt_x, txt_y = 0.99, 0.98
        t = ax.text(txt_x, txt_y, txt, horizontalalignment='right',
                verticalalignment='top', fontsize='small', zorder=3,
                transform=ax.transAxes)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    axs[2].set_xlabel(xlabel, fontsize='small')

    # make the y label
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('instrument mag', fontsize='small', labelpad=3)

    if not savpath:
        savpath = 'temp_{:s}.png'.format(apstr)

    fig.tight_layout(h_pad=-0.5)
    fig.savefig(savpath, dpi=250, bbox_inches='tight')
    print('%sZ: made plot: %s' % (datetime.utcnow().isoformat(), savpath))


def plot_lightcurve_and_ACF(
    rawlags, rawacf, epdlags, epdacf, tfalags, tfaacf,
    rawtime, rawflux, epdtime, epdflux, tfatime, tfaflux,
    ap, savpath=None,
    xlabeltime='BJD - 2457000',xlabelacf='ACF lag [days]'):
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
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(7,4),
                            gridspec_kw= {'width_ratios':[4,1],
                                          'height_ratios':[1,1,1]})
    a0,a1,a2,a3,a4,a5 = axs.flatten()

    apstr = 'AP{:d}'.format(ap)
    stagestrs = ['RM{:d}'.format(ap),
                 'EP{:d}'.format(ap),
                 'TF{:d}'.format(ap)]

    for ax, mag, time, txt in zip(
        (a0,a2,a4), [rawflux,epdflux,tfaflux], [rawtime,epdtime,tfatime],
        stagestrs
    ):
        ax.plot(time, mag, c='k', linestyle='-', marker='o',
                markerfacecolor='k', markeredgecolor='k', ms=3, lw=0.,
                zorder=1, rasterized=True)
        ax.plot(time, mag, c='orange', linestyle='-', marker='o',
                markerfacecolor='orange', markeredgecolor='orange', ms=1.5,
                lw=0., zorder=2, rasterized=True)
        txt_x, txt_y = 0.99, 0.98
        t = ax.text(txt_x, txt_y, txt, horizontalalignment='right',
                verticalalignment='top', fontsize='small', zorder=3,
                transform=ax.transAxes)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='x-small')

    for ax, acf, lag, txt in zip(
        (a1,a3,a5),
        [rawacf, epdacf, tfaacf],
        [rawlags, epdlags, tfalags],
        stagestrs
    ):

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

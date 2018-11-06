# -*- coding: utf-8 -*-
from __future__ import division, print_function
'''
functions for assessing statistics of lightcurves you've made.

whisker_MAD_stats_and_plots:
    make csv files and (optionally) whisker plots of MAD vs magnitude.

TODO:
move items from aperturephot.py here, and update all respective calls.

FIXME: break these out into their own module
1. RMS vs. MAG plot for EPD and TFA lightcurves
2. MAD vs. MAG plot for EPD and TFA lightcurves
3. ratios of RMS and MAD vs. MAG for CCD 6,7,8 to that of CCD 5
4. binned LC versions of these plots, using 10, 30, and 60 minute binning


'''

import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import aperturephot as ap
from glob import glob
from astropy import units as units, constants as constants
from datetime import datetime

################################
# LIGHTCURVE READING FUNCTIONS #
################################

def read_tfa_lc(tfafile,
                jdcol=0,
                timename='btjd',
                rmcols=[14,19,24],
                epcols=[27,28,29],
                tfcols=[30,31,32]):
    '''
    Get the times and mags from a .tfalc file. For speed, don't read any of the
    rest. The full contents of a .tfalc file are written below.

    Args:
        column ids for RAW, EPD, and TFA columns.

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
    '''

    lcmagcols = [pre+str(ix) for pre in ['RM','EP','TF'] for ix in range(1,4)]

    lcdata = np.genfromtxt(tfafile,
                           usecols=tuple([jdcol] + rmcols + epcols + tfcols),
                           names=[timename] + lcmagcols)

    return lcdata


####################
# HELPER FUNCTIONS #
####################

def _customized_box_plot(percentiles, axes, redraw = True, *args, **kwargs):
    '''
    Generates a customized boxplot based on the given percentile values. Stolen
    from

    `https://stackoverflow.com/questions/27214537/is-it-possible-to-draw-a-matplotlib-boxplot-given-the-percentile-values-instead`
    '''

    # Create len(percentiles) no of box plots
    n_box = len(percentiles)
    box_plot = axes.boxplot([[-9, -4, 2, 4, 9],]*n_box, *args, **kwargs)

    min_y, max_y = float('inf'), -float('inf')

    for box_no, (q1_start,
                 q2_start,
                 q3_start,
                 q4_start,
                 q4_end
                 ) in enumerate(percentiles):

        # Lower cap
        box_plot['caps'][2*box_no].set_ydata([q1_start, q1_start])
        # xdata is determined by the width of the box plot

        # Lower whiskers
        box_plot['whiskers'][2*box_no].set_ydata([q1_start, q2_start])

        # Higher cap
        box_plot['caps'][2*box_no + 1].set_ydata([q4_end, q4_end])

        # Higher whiskers
        box_plot['whiskers'][2*box_no + 1].set_ydata([q4_start, q4_end])

        # Box
        box_plot['boxes'][box_no].set_ydata([q2_start,
                                             q2_start,
                                             q4_start,
                                             q4_start,
                                             q2_start])

        # Median
        box_plot['medians'][box_no].set_ydata([q3_start, q3_start])

        min_y = min(q1_start, min_y)
        max_y = max(q4_end, max_y)

        # The y axis is rescaled to fit the new box plot completely with 10% 
        # of the maximum value at both ends
        axes.set_ylim([min_y*1.2, max_y*1.1])

    # If redraw is set to true, the canvas is updated.
    if redraw:
        axes.figure.canvas.draw()

    return box_plot


######################
# PLOTTING FUNCTIONS #
######################

def whisker_MAD_stats_and_plots(statdir, outprefix, binned=False,
                                make_whisker_plot=True, whisker_xlim=[4,17],
                                whisker_ylim=[1e-5,1e-1],
                                percentiles=[2,25,50,75,98]):
    '''
    make csv files and (optionally) plots of MAD vs magnitude statistics,
    evaluated at percentiles.
    '''

    # get and read the most complete stats file (the TFA one!)
    tfastatfile = glob(statdir+'*.tfastats')

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

    # MAD vs magnitude statistics, at 1-mag differences, summarized as whiskers
    # at 2%, 25%, 50%, 75%, and 98% percentiles.
    for apstr in ['tf1','tf2','tf3',
                  'rm1','rm2','rm3',
                  'ep1','ep2','ep3']:

        medstr = 'med_'+apstr
        madstr = 'mad_'+apstr

        minmag = np.floor(stats[medstr].min()).astype(int)
        maxmag = np.ceil(stats[medstr].max()).astype(int)
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
                val = np.percentile(stats[sel][madstr], percentile)
                percentile_dict[thismagmean][percentile] = np.round(val,7)

        pctile_df = pd.DataFrame(percentile_dict)

        if make_whisker_plot:

            plt.close('all')
            fig, ax = plt.subplots(figsize=(4,3))

            # whisker plot isn't aware of absolute x values, it just sees them
            # as different columns. thus we need to pre-select the magnitude
            # columns according to whisker_xlims that are passed.
            if whisker_xlim:
                sel_cols = (
                    (pctile_df.columns < max(whisker_xlim)) &
                    (pctile_df.columns > min(whisker_xlim))
                )
                # surprising that i couldn't figure out a smarter way ._.
                plotdf = pctile_df.transpose()[sel_cols].transpose()

            boxdata = np.array(plotdf).T

            b = _customized_box_plot(
                boxdata, ax, redraw=True, notch=0, vert=1, whis=1.5
            )

            ax.set_yscale('log')
            ax.set_xlabel('{:s} median instrument magnitude'.format(apstr.upper()))
            ax.set_ylabel('{:s} median abs. dev.'.format(apstr.upper()))
            xticklabels = []
            for ix, m in enumerate(np.mean(mag_bins, axis=1)):
                if ix%2==0:
                    xticklabels.append(m)
                else:
                    xticklabels.append('')
            ax.set_xticklabels(xticklabels)
            if whisker_ylim:
                ax.set_ylim(whisker_ylim)

            titlestr = '{:s} - {:d} - {:s}'.format(
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

            savname = (
                statdir+'whisker_MAD_vs_med_mag_{:s}.png'.format(apstr.upper())
            )
            fig.tight_layout()
            fig.savefig(savname, dpi=250)
            print('%sZ: made %s plot: %s' %
                  (datetime.utcnow().isoformat(), titlestr, savname))

        csvname = (
            statdir+'whisker_MAD_vs_med_mag_{:s}.csv'.format(apstr.upper())
        )
        pctile_df.to_csv(csvname, index=False)
        print('%sZ: wrote %s' %
              (datetime.utcnow().isoformat(), csvname))


def plot_raw_epd_tfa(time, rawmag, epdmag, tfamag, ap_index, savpath=None,
                     xlabel='BTJD = BJD - 2457000'):
    '''
    Plot a 3 row, 1 column plot with rows of:
        * raw mags vs time
        * EPD mags vs time
        * TFA mags vs time.

    args:
        time, rawmag, epdmag, tfamag (np.ndarray)

        ap_index (int): integer, e.g., "2" for aperture #2.
    '''

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

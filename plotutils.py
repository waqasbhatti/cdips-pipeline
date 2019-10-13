'''
plotutils.py

Various useful plotting functions.
'''

import re
import os
from glob import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from matplotlib import colors, patches

import imageutils as iu
import astropy.io.fits as pyfits


#################
#  RMS and MAD  #
#################
def get_plot_label(name, type='legend'):
    '''
    Utility function to get a properly formatted plot label.
    '''
    label = ''
    if name == 'catmag':
        return 'Catalog mag'
    if 'median' in name:
        label = 'Median '
    elif 'mean' in name:
        label = 'Mean '

    if ('rm' in name or 'ep' in name or 'tf' in name) \
            and type =='general':
        label += name[-3:].upper() + ' '

    if type == 'legend':
        if 'mad' in name:
            label += 'MAD '
        elif 'stdev' in name:
            label += 'RMS '
    elif type == 'ylabel':
        if 'mad' in name:
            label += 'Median absolute deviation'
        elif 'stdev' in name:
            label += 'Standard deviation'

    if 'sigclip' in name:
        label += '(sigclip)'

    if type == 'xlabel':
        label += 'Magnitude'

    return label

def make_rms_plot(stats, xcols, ycols, logy=True,
                  xlim=(6,14), ylim=None, figsize=(8,6),
                  xlabel=None, ylegends=None, ylabel=None,
                  title=None, newfig=True, savefig=False):
    '''
    Makes rms or MAD plots.

    Args:
        stats: np.recarray or DataFrame-like object.
        xcol (str): Column containing x variable.
        ycols (str, callable, or list): y variable.
            Can be a column name, a callable to be applied to stats,
            or a list. If a list is provided, will overlay plots.
        xlabel: Label to use for x axis.
        ylegends (str or list): Labels to use in legend.
        newfig: Creates a new figure
        savefig: Filename to save the figure.
    '''
    if newfig:
        fig = plt.figure(figsize=figsize)
    else:
        fig = plt.gcf()
    if type(ycols) != list:
        ycols = [ycols]
    if type(xcols) != list:
        xcols = [xcols] * len(ycols)
    elif len(xcols) == 1:
        xcols = [xcols] * len(ycols)
    if len(xcols) != len(ycols):
        print("Must provide either one xcol or as many xcols as ycols.")
        raise ValueError

    ax = plt.gca()
    if ylegends is None:
        ylegends = []
        for ycol in ycols:
            if callable(ycol):
                ylegends.append('')
            else:
                ylegends.append(get_plot_label(ycol))
    elif type(ylegends) != list:
        ylegends = [ylegends]

    for xcol, ycol, ylegend in zip(xcols, ycols, ylegends):
        if callable(ycol):
            ydata = ycol(stats)
            ax.scatter(stats[xcol], ydata, s=2, alpha=0.5, label=ylegend)
        else:
            ax.scatter(stats[xcol], stats[ycol], s=2, alpha=0.5, label=ylegend)

    if ylegends is not False:
        leg = ax.legend(markerscale=3)
        for lh in leg.legendHandles:
            lh.set_alpha(1.0)
    if logy:
        ax.set_yscale('log')
        if ylim is None:
            ylim = (10**-3, 1)
    else:
        if ylim is None:
            ylim = ax.get_ylim()

    # Get default x labels
    if xlabel is None:
        xlabel = get_plot_label(xcols[0], type='xlabel')
    ax.set_xlabel(xlabel)
    if ylabel is None:
        ylabel = get_plot_label(ycols[0], type='ylabel')
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid('on', linestyle=':', which='both', axis='both')

    if title is not None:
        plt.title(title, fontsize=14)

    plt.tight_layout()

    if savefig:
        plt.savefig(savefig)

    return fig, ax

def make_percentile_plot(stats, xcols, ycols, logy=True,
                         percentiles=[.02, .25, .50, .75, .98],
                         xlim=(6,14), ylim=None, figsize=(8,6),
                         xlabel=None, ylegends=None, ylabel=None,
                         title=None, newfig=True, savefig=False):
    '''
    Bins the results by magnitude and creates line plots for each percentile.
    '''
    stats = pd.DataFrame(stats)
    if newfig:
        fig = plt.figure(figsize=figsize)
    else:
        fig = plt.gcf()
    ax = plt.gca()
    linestyles = ['-', '--', ':', '-.']

    # Get y column labels
    if type(ycols) != list:
        ycols = [ycols]
    if ylegends is None:
        ylegends = []
        for ycol in ycols:
            if callable(ycol):
                ylegends.append('')
            else:
                ylegends.append(get_plot_label(ycol))
    elif type(ylegends) != list:
        ylegends = [ylegends]
    if type(xcols) != list:
        xcols = [xcols] * len(ycols)
    elif len(xcols) == 1:
        xcols = [xcols] * len(ycols)
    if len(xcols) != len(ycols):
        print("Must provide either one xcol or as many xcols as ycols.")
        raise ValueError

    for idx, (xcol, ycol, ylegend) in enumerate(zip(xcols, ycols, ylegends)):
        if callable(ycol):
            ydata = ycol(stats)
        else:
            ydata = stats[ycol]
        # Compute percentiles
        grouped_data = ydata.groupby(stats[xcol].round())
        for pidx, percentile in enumerate(percentiles):
            plabel = ylegend if int(pidx - len(percentiles)/2) == 0 else '_'
            ax.plot(grouped_data.quantile(percentile), 'o-',
                    ms=2, color='C{:d}'.format(idx % 10), label=plabel,
                    linestyle=linestyles[abs(int(pidx - len(percentiles)/2))])

    # Set scale and legend
    ax.legend()
    if logy:
        ax.set_yscale('log')
        if ylim is None:
            ylim = (10**-3, 1)
    else:
        if ylim is None:
            ylim = ax.get_ylim()

    # Get default x labels
    if xlabel is None:
        xlabel = get_plot_label(xcols[0], type='xlabel')
    ax.set_xlabel(xlabel)
    if ylabel is None:
        ylabel = get_plot_label(ycols[0], type='ylabel')
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid('on', linestyle=':', which='both', axis='both')

    if title is not None:
        plt.title(title, fontsize=14)

    plt.tight_layout()

    if savefig:
        plt.savefig(savefig)
        plt.close()
        return savefig

    return fig, ax


#################
#   Utilities   #
#################
def get_magbin_subset(stats, magcol, num_stars=50,
                      varcol=None, varmin=0.0, varmax=0.1,
                      permin=0.0, permax=None,
                      magmin=5.5, magmax=14.5):
    '''
    Get a subset of stars from each magnitude bin.

    Args:
        stats (pd.DataFrame): statistics table
        magcol (str): Column with magnitudes to bin on
        num_stars (int): Number of stars in each bin
        varcol (str): If not None, column to query on
        varmin, varmax: Query bounds.
    Returns: list of (mag, pd.DataFrame) tuples
    '''
    if type(stats) != pd.DataFrame:
        stats = pd.DataFrame(stats)
    if magmin is not None and magmax is not None:
        stats = stats.query('{:f} <= {:s} <= {:f}'.format(
            magmin, magcol, magmax))
    grouped_data = stats.groupby(stats[magcol].round())

    stars = []
    for mag, group in grouped_data:
        if varcol is not None:
            # Use percentile based bounds
            if permin is not None and permax is not None:
                bounds = group[varcol].quantile([permin, permax])
                bounds = [bounds.iloc[0], bounds.iloc[1]]
            else:
                bounds = [varmin, varmax]
            subset = group.query('{:f} <= {:s} < {:f}'.format(
                bounds[0], varcol, bounds[1]))

            if num_stars == 'all':
                stars.append((mag, subset))
            else:
                stars.append((mag, subset.iloc[0:num_stars]))
        else:
            if num_stars == 'all':
                stars.append((mag, group))
            else:
                stars.append((mag, group.iloc[0:num_stars]))

    return stars

def read_lc_mags(lcfile, jdcol=0, timename='rjd',
            framecol=None, framename='rstfc',
            rmcols=[14,19,24], epcols=[27,28,29],
            tfcols=[30,31,32], rfcols=[12,17,22],
            recols=[15,20,25]):
    '''
    Read only magnitude (and error) columns for LC.
    '''
    # Automatically detect type of file
    with open(lcfile, 'r') as f:
        line = f.readline()
        numcols = len(line.split())
    if epcols is not None and numcols <= epcols[-1]:
        epcols = None
    if tfcols is not None and numcols <= tfcols[-1]:
        tfcols = None

    num_aps = len(rmcols)
    colnames = [timename,]
    cols = [jdcol,]
    if framecol is not None:
        colnames += [framename]
        cols += [framecol]
    if rmcols is not None:
        colnames += ['rm%d' % (i + 1) for i in range(num_aps)]
        cols += rmcols
    if epcols is not None:
        colnames += ['ep%d' % (i + 1) for i in range(num_aps)]
        cols += epcols
    if tfcols is not None:
        colnames += ['tf%d' % (i + 1) for i in range(num_aps)]
        cols += rfcols
    if rfcols is not None:
        colnames += ['rf%d' % (i + 1) for i in range(num_aps)]
        cols += rfcols
    if recols is not None:
        colnames += ['re%d' % (i + 1) for i in range(num_aps)]
        cols += recols

    lcdata = np.genfromtxt(lcfile, usecols=cols, names=colnames)

    return lcdata


def mad(mags):
    return np.nanmedian(np.fabs(mags - np.nanmedian(mags)))

def sigmaclip(lc, magcol, sigclip=5.0):
    '''
    Iterative sigma-clipping
    '''
    std = np.nanstd(lc[magcol])
    mn = np.nanmean(lc[magcol])
    sigmask = np.logical_or(lc[magcol] < (-sigclip * std + mn),
                            lc[magcol] > (sigclip * std + mn))
    while np.sum(sigmask) > 0:
        lc = lc[~sigmask]
        std = np.nanstd(lc[magcol])
        mn = np.nanmean(lc[magcol])
        sigmask = np.logical_or(lc[magcol] < (-sigclip * std + mn),
                                lc[magcol] > (sigclip * std + mn))
    return lc

def finiteclip(lc, magcol):
    '''
    Clip lightcurve down to rows with finite magnitudes.
    '''
    finitemask = np.isfinite(lc[magcol])
    return lc[finitemask]


#################
#   LC Plots    #
#################
def best_magticks(ymin, ymax):
    '''
    Get best magticks for a certain magrange
    '''
    if (ymax - ymin) < 0.5:
        ystep = 0.1
    elif (ymax - ymin) < 2.5:
        ystep = 0.25
    else:
        ystep=0.5
    return (np.arange(ymin, ymax + 5e-2, ystep))


def plot_lc(lc, aps='all', apnums='all', tcol='rjd', tlim=None, sigclip=5.0,
            title=None, figheight=3, figwidth=8, savefig=False):
    '''
    Plot a lightcurve.

    Args:
        lc (np.recarray or pd.DataFrame): Lightcurve
        aps (list or 'all'): List of apertures (e.g. ['rm', 'ep'])
        apnums (list or 'all'): List of apnums (e.g. [1, 2, 3])
        figheight: Height of single lightcurve
        figwidth: Width of figure
    '''
    # Get number of plots to make
    lccols = list(lc.columns) if type(lc) == pd.DataFrame \
                              else list(lc.dtype.names)
    if type(aps) == str and aps == 'all':
        aps = [c for c in ['rm', 'ep', 'tf'] if c + str(1) in lccols]
    elif type(aps) == str:
        aps = [aps]
    if type(apnums) == int:
        apnums = list(range(1, apnums+1))
    elif type(apnums) == str and apnums == 'all':
        apnums = [int(c[-1:]) for c in lccols if aps[0] in c]

    lccols = [ap + str(i) for ap in aps for i in apnums]
    num_plots = len(aps) * len(apnums)

    fig, axs = plt.subplots(nrows=num_plots, ncols=1, sharex=True,
                            figsize=(figwidth, figheight*num_plots))

    if tlim is not None:
        tmask = (lc[tcol] > tlim[0]) & (lc[tcol] < tlim[1])
        lc = lc[tmask]

    for (ax, ap) in zip(axs, lccols):
        lc_ap = finiteclip(lc, ap)
        if sigclip:
            lc_ap = sigmaclip(lc, ap, sigclip)

        ax.scatter(lc_ap[tcol], lc_ap[ap], c='k', alpha=0.9, s=3, linewidths=0)
        ax.text(0.99, 0.95, ap.upper(),
                fontsize='medium',transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='top')
        ax.text(0.99, 0.87, 'MAD = %f' % mad(lc_ap[ap]),
                fontsize='medium', transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='top')
        ax.text(0.99, 0.79, 'RMS = %f' % np.nanstd(lc_ap[ap]),
                fontsize='medium', transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='top')
        yticks = ax.get_yticks()
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small')
        ax.set_yticks(best_magticks(yticks[1], yticks[-2]))

        ax.set_ylim(reversed(ax.get_ylim()))
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small')

    axs[-1].set_xlabel(tcol.upper(), fontsize='large')
    ax_hidden = fig.add_subplot(111, frameon=False)
    ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False)
    ax_hidden.set_ylabel('instrument mag', fontsize='large', labelpad=5)
    if title is not None:
        axs[0].set_title(title, fontsize='large')
    fig.tight_layout(h_pad=-1.0)

    if savefig:
        if savefig.endswith('.png')
            fig.savefig(savfig, bbox_inches='tight', dpi=250)
        elif savefig.endswith('.png')
            fig.savefig(savfig, bbox_inches='tight')
        else:
            fig.savefig(savfig, bbox_inches='tight')
        plt.close()

    return fig, ax


#################
#   ACF Plots   #
#################
def plot_dacf(dacf, overlayplots=False, showerrs=False, title=None,
              ylegends=None, figheight=4, figwidth=8, savefig=False):
    '''
    Plot a DACF object.

    Args:
        dacf (dict): Result from lcs.compute_dacf()
        overlayplots (bool): If true, overlays all ACFs; otherwise makes
            subplots.
        showerrs (bool): Whether to show error bars
    '''
    udcfarr = dacf['udcf']
    eudcfarr = dacf['eudcf']
    if udcfarr.ndim == 1:
        udcfarr = [udcfarr]
        eudcfarr = [eudcfarr]

    if overlayplots or len(udcfarr) == 1:
        fig = plt.figure(figsize=(figwidth, figheight))
        ax = plt.gca()
        subplots = False
    else:
        num_plots = len(udcfarr)
        fig, axs = plt.subplots(nrows=num_plots, ncols=1, sharex=True,
                                figsize=(figwidth, figheight*num_plots))
        subplots = True

    # Plot each dcf
    for (i, udcf) in enumerate(udcfarr):
        yerr = eudcfarr[i] if showerrs else None
        ylegend = ylegends[i] if ylegends is not None else None
        if subplots:
            ax = axs[i]
            color='k'
        else:
            color='C%d' % (i%10)
        ax.errorbar(dacf['lags'], udcf, yerr=yerr, ecolor='gray', 
                    marker='o', mec='none', ms=3, color=color,
                    label=ylegend)
        ax.grid('on', linestyle='dashed')
        ax.set_xlim(0, dacf['lags'][-1])
        ax.axhline(0, color='0.5', zorder=0)
        if ylegends is not None:
            ax.legend()

    if subplots:
        if title is not None:
            axs[0].set_title(title, fontsize='large')
        axs[-1].set_xlabel('Lags (days)', fontsize='medium')
        ax_hidden = fig.add_subplot(111, frameon=False)
        ax_hidden.tick_params(labelcolor='none', top=False, bottom=False,
                              left=False, right=False)
        ax_hidden.set_ylabel('Autocorrelation', fontsize='medium', labelpad=5)
        fig.tight_layout(h_pad=-1.0)
    else:
        if title is not None:
            ax.set_title(title, fontsize='large')
        ax.set_xlabel('Lags (days)', fontsize='medium')
        ax.set_ylabel('Autocorrelation')

    if savefig:
        plt.savefig(savefig)
        plt.close()

    return fig, ax


#################
#     Stamps    #
#################
def get_obj_coords(iphot, objlist, coord_cols=[1,2]):
    obj_coords = {}
    if type(objlist) is not list:
        objlist = [objlist]
    for line in open(iphot):
        splline = line.split()
        if splline[0] in objlist:
            obj_coords[splline[0]] = tuple(float(splline[i])
                                           for i in coord_cols)

    for obj in objlist:
        if obj not in obj_coords:
            obj_coords[obj] = None

    return obj_coords


def get_stamp_bounds(coords, size, limits):
    '''
    Gets the coordinates of a stamp around coords, with provided size.

    Args:
        coords (tuple): Central coords
        size (int): Size of square stamp
        limits (tuple): Dimensions of full fits object array
    '''
    xmin = int(round(coords[0] - size))
    xmin = xmin if xmin > 0 else 0
    xmax = int(round(coords[0] + size))
    xmax = xmax if (xmax < limits[0]) else limits[0]
    ymin = int(round(coords[1] - size))
    ymin = ymin if ymin > 0 else 0
    ymax = int(round(coords[1] + size))
    ymax = ymax if ymax < limits[1] else limits[1]

    return xmin, xmax, ymin, ymax


def shift_coords(coords, x, y):
    '''
    Shift coordinate system to one with origin at (x, y)
    '''
    return (coords[0] - x - 0.5, coords[1] - y - 0.5)


def draw_aperture(aperturestr, center, ax, colors='gray'):
    '''
    Aperture should be in same format as call to fiphot
    '''
    if aperturestr is None:
        return
    apertures = aperturestr.split(',')

    if colors =='gray':
        colors = ['r', 'b', 'b']
    elif colors == 'rgb':
        colors = ['k', 'w', 'w']

    for aperture in apertures:
        ap_radii = np.asarray(aperture.split(':'), dtype=float)
        ap_radii[2] = ap_radii[1] + ap_radii[2]

        for r, c in zip(ap_radii, colors):
            circ = patches.Circle(center, r, facecolor='none', edgecolor=c)
            ax.add_patch(circ)


def dump_stamps(objs, fits, iphot=None, fitstype='rsub', coord_cols=None,
                ofits=None, oiphot=None, includeoriginal=False,
                size=15, aperturestr=None, showcolorbar=True,
                scale_func=None, scale_func_params={},
                figsize=6, header='', savefig=False, computeresiduals=True,
                noplot=False):
    '''
    Dumps stamps around a list of objects from a given fits file.

    Args:
        objs (list): List of objects (and optionally, their coordinates)
        fits (str): Path to fits file.
        iphot (str): Path to iphot file.
            Not used if a list of coordinates is already provided.
            If None, searches for an appropriate iphot according to fitstype.
        fitstype (str): If rsub or nsub, helps in searching for iphot.
        ofits (str): Original (non image-subtracted) xtrns fits file. 
        oiphot (str): Original iphot.
        size (int): Size of stamp.
        aperturestr (str): aperture str (same format as fiphot)
    '''
    if fits is None:
        return None
    fitsframename = os.path.splitext(os.path.basename(fits))[0]
    if fitstype == 'rsub' or fitstype == 'nsub':
        fitsframename = re.sub('[r|n]sub-.*?-', '', fitsframename)\
                            .replace('-xtrns', '')
    # Get coordinates if they are not provided.
    coords_provided = True
    if not type(objs) == dict:
        coords_provided = False
    else:
        for obj in objs:
            if type(objs[obj]) != tuple:
                coords_provided = False
            break
    if not coords_provided:
        # Find iphot if not provided.
        if iphot is None:
            if fitstype == 'rsub':
                iphotglob = re.sub('rsub-.*?-', 'rsub-*-', fits)\
                                .replace('-xtrns.fits', '.iphot')
            elif fitstype == 'nsub':
                iphotglob = re.sub('nsub-.*?-', 'nsub-*-', fits)\
                                .replace('-xtrns.fits', '.iphot')
            else:
                iphotglob = fits.replace('.fits', '.fiphot')

            try:
                iphot = glob(iphotglob)[0]
            except:
                print('Could not find iphot file %s.' % iphotglob)
                return None
        # Read iphot to get coords
        # Coordinate columns are [1, 2] in rsub iphot and [2, 3] in fiphot.
        if coord_cols is None and fitstype == 'rsub' \
                               or fitstype == 'nsub':
            coord_cols = [1, 2]
        elif coord_cols is None:
            coord_cols = [2, 3]
        objs = get_obj_coords(iphot, objs, coord_cols)

    # Get original (non subtracted) fits if requested
    ofits = None
    oiphot = None
    if (fitstype == 'rsub' or fitstype == 'nsub') and includeoriginal:
        if includeoriginal == 'xtrns':
            ofitsglob = re.sub('[r|n]sub-.*?-', '', fits)
            ofits_objs = objs
        else:
            ofitsglob = re.sub('[r|n]sub-.*?-', '', fits)\
                                .replace('-xtrns.fits', '.fits')
            # Also need to find oiphot
            oiphotglob = ofitsglob.replace('.fits', '.fiphot')
            oiphotglob = glob(oiphotglob)
            oiphot = oiphotglob[0] if oiphotglob is not None else None

            # Read in coordinates
            if oiphot is not None:
                objlist = [obj for obj in objs]
                ofits_objs = get_obj_coords(oiphot, objlist, [2, 3])
            else:
                ofits_objs is None

        oglob = glob(ofitsglob)
        ofits = oglob[0] if oglob is not None else None
        if ofits is None:
            print(('Could not find original fits file for %s, ' % fits) +
                  'continuing without... ')
        elif ofits_objs is None:
            ofits = None
            print(('Could not find original fiphot file for %s, ' % fits) +
                  'continuing without...')
        else:
            ofitsarr = pyfits.getdata(ofits)

    # Read fits
    fitsarr = pyfits.getdata(fits)

    # Plot for each object
    fluxresiduals = []
    for obj in objs:
        coords = objs[obj]
        if coords is None:
            continue

        xmin, xmax, ymin, ymax = get_stamp_bounds(coords, size, fitsarr.shape)
        corrected_center = shift_coords(coords, xmin, ymin)

        if computeresiduals:
            resid_mean = np.mean(fitsarr[ymin:ymax, xmin:xmax])
            resid_rms = np.std(fitsarr[ymin:ymax, xmin:xmax])
            fluxresiduals.append([fitsframename, obj, resid_mean, resid_rms])

        if noplot:
            continue

        # Create figure
        if ofits is None or ofits_objs[obj] is None:
            fig = plt.figure(figsize=(figsize, figsize))
            ax = plt.gca()
        else:
            fig, axes = plt.subplots(nrows=1, ncols=2,
                                     figsize=(figsize*2,figsize))
            ax = axes[0]

        if fitstype == 'rsub':
            cscale = np.max(np.abs(fitsarr[ymin:ymax, xmin:xmax]))
            cscale = 10**(np.floor(np.log10(cscale)) + 1)
            diffnorm = colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                         vmin=-cscale, vmax=cscale)
            pimg = ax.imshow(fitsarr[ymin:ymax, xmin:xmax], cmap='RdBu_r',
                             vmin=-cscale, vmax=cscale, norm=diffnorm)
            ax.set_aspect('equal')
            ax.set_title(fitstype + ' frame')
            plt.tight_layout()
            if showcolorbar:
                if includeoriginal:
                    cb_ax = fig.add_axes([0.1, 0.1, 0.02, 0.8], frameon=False)
                    leftpad = 0.15
                else:
                    cb_ax = fig.add_axes([0.1, 0.1, 0.04, 0.8], frameon=False)
                    leftpad = 0.2
                cb = fig.colorbar(pimg, cax=cb_ax, extend='both')
                cb_ax.yaxis.set_ticks_position('left')
                cb.set_ticks([-1e3,-1e2,-1e1,0,1e1,1e2,1e3])
                cb.set_ticklabels(['-$10^3$','-$10^2$','-$10^1$','0',
                                    '$10^1$','$10^2$','$10^3$'])
                fig.subplots_adjust(left=leftpad)
        else:
            if scale_func is None:
                ax.imshow(fitsarr[ymin:ymax, xmin:xmax],
                          cmap='gray', vmin=0, vmax=255)
            else:
                ax.imshow(scale_func(fitsarr[ymin:ymax, xmin:xmax],
                                     **scale_func_params),
                          cmap='gray', vmin=0, vmax=255)
            ax.set_title(fitstype + ' frame')

        draw_aperture(aperturestr, corrected_center, ax, colors='rgb')

        # Plot original fits
        if ofits is not None and ofits_objs[obj] is not None:
            ax = axes[1]

            ofits_coords = ofits_objs[obj]
            xmin, xmax, ymin, ymax = get_stamp_bounds(ofits_coords, size,
                                                      ofitsarr.shape)
            corrected_center = shift_coords(ofits_coords, xmin, ymin)

            ax.imshow(iu.clipped_linscale_img(ofitsarr[ymin:ymax, xmin:xmax]),
                      vmin=0, vmax=255, cmap='gray')
            ax.set_aspect('equal')
            ax.set_title('original frame', fontsize='large')

            draw_aperture(aperturestr, corrected_center, ax, colors='gray')

        # Make title
        plt.suptitle(obj + '\n' + fitsframename + header,
                     fontsize='x-large', fontweight='bold',
                     ha='center', va='center')

        if savefig:
            outdir = os.path.join(savefig, obj + '/')
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            outplotfile = os.path.join(outdir, fitsframename + '.png')
            plt.savefig(outplotfile)
            plt.close()
            print('Saved figure to %s' % outplotfile)

    if computeresiduals and fitstype == 'rsub':
        return fluxresiduals
    else:
        return None


###################
#     Dump JPG    #
###################
def fits_to_full_jpeg(fits, fitstype='rsub', outfile=None, outprefix='',
                      flip=None, annotate=True, fits_jdsrc=None,
                      frame_time=None, colorscheme=None, norm=None,
                      scale_func=None, scale_func_params={},
                      fieldinfo=None, showfig=False, overwrite=False):
    '''
    Converts a FITS image to a full frame JPEG.

    Args:
        flip: Axis to flip along
    '''
    # Get outfilename
    if fits is None:
        return None
    fitsframename = os.path.splitext(os.path.basename(fits))[0]
    if outfile is None:
        outdir = os.path.dirname(fits)
        outfile = os.path.join(outdir, outprefix + fitsframename + '.jpg')
        if os.path.exists(outfile) and not overwrite:
            print('%s already exists , continuing...' % outfile)

    img, hdr = pyfits.getdata(fits, header=True)
    if fieldinfo is not None:
        projid = fieldinfo['projectid']
        field = fieldinfo['field']
    else:
        projid = hdr['PROJID'] if 'PROJID' in hdr else 'unknown'
        field = hdr['OBJECT'] if 'OBJECT' in hdr else 'objunknown'
    exptime = hdr['EXPTIME'] if 'EXPTIME' in hdr else 'expunknown'

    # Default values
    if fitstype == 'rsub' or fitstype == 'nsub':
        if colorscheme is None:
            colorscheme = 'RdBu_r'
        if norm is None:
            norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                     vmin=-10000, vmax=10000)
    else:
        if colorscheme is None:
            colorscheme = 'gray'
        if scale_func is None:
            scale_func = iu.clipped_linscale_img

    # Scale image values
    if scale_func is not None:
        img = scale_func(img, **scale_func_params)

    # Flip image
    if flip is not None:
        img = np.flip(img, flip)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0., 0., 1., 1.,], frameon=False)
    ax.set_axis_off()
    ax.imshow(img, cmap=colorscheme, norm=norm)

    # Annotate
    if annotate:
        annotation = "%s; PROJ%s - %s; Exp: %s s" % (
            fitsframename, projid, field, exptime)
        txtcolor = 'k' if colorscheme is not 'gray' else 'w'
        ax.text(15, 15, annotation, fontsize='x-large',
                va='top', color=txtcolor)

    plt.savefig(outfile)
    if not showfig:
        plt.close()



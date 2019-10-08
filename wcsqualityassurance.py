"""
Tools to tell you how good your WCS solutions are. In particular, how precisely
can you convert between sky and pixel coordinates? (If worse than your aperture
size, this is a big no-no).
"""

import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, shutil
from numpy import array as nparr

from astropy.io import fits
from astropy import wcs

from aperturephot import (
    make_frameprojected_catalog,
    extract_frame_sources,
    astrometrydotnet_solve_frame
)

def does_wcs_pass_quality_check(
    wcsfile,
    fitsfile,
    reformedfovcatalog,
    isspocwcs=True,
    N_bright=1000,
    N_faint=9000,
    fistarpath=None,
    matchedoutpath=None,
    make_plots=True,
    qualitycondition={'median_px':0.2, '90th_px':0.4, 'std_px': 1.7},
    ccdextent={'x':[0.,2048.],'y':[0.,2048.]},
    pixborders=0.0,
    gain=5.35,
    fluxthreshold=500,
    zeropoint=11.82,
    exptime=1800
):
    """
    You have a WCS for some image.  You want to be sure that the WCS
    transformation is good to a required precision, e.g., to do forced aperture
    photometry.

    This function projects Gaia (ra,dec) to pixel coordinates using the passed
    WCS.  The (x,y) positions are also measured with fistar (and either passed
    via fistarpath, or else found by this function using fistar).

    The function then matches the projected positions to the measured positions
    with the `grmatch` symmetric point-matching algorithm.  The sources chosen
    for this are in the top N_bright to N_faint stars in both lists.

    It can optionally then make three plots to assess the quality of the WCS
    solution:
        * a histogram of the separation distribution
        * a scatter plot with the separations colored
        * a vector plot with the position differences as arrows and colored

    Args:

        wcsfile : path to .wcs file

        fitsfile : path to .fits file the WCS was derived for

        reformedfovcatalog : made by `ap.reform_gaia_fov_catalog`

        isspocwcs : bool

        N_bright : integer brightest number of stars for grmatch. (For TESS,
        1000 avoids saturated stars).

        N_faint: grmatch down to how faint?

        fistarpath : if None, this function makes the fistar file. If a path is
        given, it is assumed to contain measured star positions.

        matchedoutpath : if None, the file with the matched fistar and
        projected coordinates is written with a "*.matched" extension. Else the
        given path is used.

        qualitycondition : dict with lower acceptable limits on astrometric
        residual separation distribution's median, 90th percentile, and
        standard deviation.

        From
            ccdextent={'x':[0.,2048.],'y':[0.,2048.]},
            pixborders=0.0,
            gain=5.35,
            fluxthreshold=500,
            zeropoint=11.82,
            exptime=1800
        first two are used for frame projection, last 4 are used for fistar
        source extraction if not already done.


    Returns:

        True if all the qualitycondition conditions are ok, else False.
    """

    outprojcat = fitsfile.replace('.fits','.projcat')

    if matchedoutpath is None:
        matchedoutpath = fitsfile.replace('.fits','.matched')

    if not os.path.exists(wcsfile):

        _wcs = os.path.join(
            os.path.dirname(wcsfile),
            'badframes',
            os.path.basename(wcsfile)
        )

        if os.path.exists(_wcs):
            return True
        else:
            raise AssertionError('did not find wcsfile')

    #
    # given WCS and reformed celestial coordinate catalog, make the frame
    # projected catalog
    #
    projcatfile = make_frameprojected_catalog(fitsfile,
                                              reformedfovcatalog,
                                              ccdextent=ccdextent,
                                              out=outprojcat,
                                              wcs=wcsfile,
                                              removetemp=True,
                                              pixborders=pixborders)
    #
    # if source extraction not done, then do it
    #
    if not os.path.exists(fistarpath):
        extract_frame_sources(fitsfile,
                              fistarpath,
                              fistarexec='fistar',
                              ccdextent='0:2048,0:2048',
                              ccdgain=gain,
                              fluxthreshold=fluxthreshold,
                              zeropoint=zeropoint,
                              exptime=exptime)


    quality_check_result = impose_wcs_quality_check(
        fitsfile,
        fistarpath,
        matchedoutpath,
        isspocwcs=isspocwcs,
        projcatfile=projcatfile,
        N_bright=N_bright,
        N_faint=N_faint,
        make_plots=make_plots,
        qualitycondition=qualitycondition
    )

    return quality_check_result


def impose_wcs_quality_check(
    fitsfile,
    fistarpath,
    matchedoutpath,
    isspocwcs=False,
    projcatfile=None,
    N_bright=1000,
    N_faint=9000,
    make_plots=True,
    qualitycondition=None
):
    """
    impose check as described in does_wcs_pass_quality_check docstring
    """

    #
    # see what stars are projected on silicon; sort by Gaia Rp magnitude; take
    # stars from N_bright to N_faint indices; write to a "cut" version of
    # projcatfile
    #
    df = pd.read_csv(
        projcatfile,
        delim_whitespace=True,
        names='id,ra,dec,xi,eta,G,Rp,Bp,plx,pmra,pmdec,varflag,x,y'.split(',')
    )

    outdf = df.sort_values(by='Rp').iloc[N_bright:N_faint]
    cutprojcat = os.path.join(
        os.path.dirname(projcatfile),
        'cut_{}'.format(os.path.basename(projcatfile))
    )
    outdf.to_csv(cutprojcat, index=False, sep=' ', header=False)

    #
    # same cut, but for stars MEASURED to be on silicon 
    #
    with open(fistarpath, 'r') as f:
        fistarlines = f.readlines()

    fistaroutlines = fistarlines[6+N_bright:6+N_faint] # 6 comment lines
    cutfistar = os.path.join(
        os.path.dirname(fistarpath),
        'cut_{}'.format(os.path.basename(fistarpath))
    )

    with open(cutfistar, 'w') as f:
        f.writelines(fistaroutlines)

    #
    # run `grmatch` to match the xy positions of the two (x,y) lists. Andras'
    # thesis describes the triangles that accomplish this.  read in the merged
    # grmatch output.
    #
    matchcmd = (
        'grmatch --match-points -i {cutfistar} --col-inp 2,3 '
        '-r {cutprojcat} --col-ref 13,14 --output {matched}'.
        format(cutfistar=cutfistar, cutprojcat=cutprojcat, matched=matchedoutpath)
    )

    rc = os.system(matchcmd)
    print(matchcmd)

    if rc != 0:
      raise RuntimeError('grmatch failed somehow')

    df = pd.read_csv(
        matchedoutpath,
        delim_whitespace=True,
        names=('id,ra,dec,xi,eta,G,Rp,Bp,plx,pmra,pmdec,varflag,x_proj,y_proj,'+
               'Ident,x_meas,y_meas,Bg,Amp,S,D,K,Flux,S/N').split(',')
    )

    #
    # projected and measured positions seem to systematically differ by 0.5
    # pixels in row and column position. NOTE: ensure this is OK.
    #
    if isspocwcs:

        df['x_proj'] -= 0.5
        df['y_proj'] -= 0.5

        #
        # SPOC rows/cols were trimmed, which further shifts WCS
        #
        # FIXME FIXME this needs to be generalized elsewhere, to ensure that the
        # projected positions used for forced aperture photometry are
        # correctedly shifted (and correctly span 0:2048, instead of like
        # 0:2048-SCIROWS
        hdulist = fits.open(fitsfile)
        hdr = hdulist[0].header

        df['x_proj'] -= (hdr['SCCSA']-1)
        df['y_proj'] -= (hdr['SCIROWS']-1)

        hdulist.close()

    else:
        df['x_proj'] -= 0.5
        df['y_proj'] -= 0.5

    df['sep'] = sep(
        df['x_meas'],
        df['y_meas'],
        df['x_proj'],
        df['y_proj']
    )

    #
    # make plots, and return whether or not WCS is precise enough
    #
    if make_plots:

        plt.close('all')

        pre = os.path.splitext(os.path.basename(fistarpath))[0]
        if isspocwcs:
            pre = pre+'_spocwcs'
        outpath = os.path.join(os.path.dirname(fistarpath),
                               '{}_sep_hist.png'.format(pre))
        plot_sep_hist(df, outpath)

        outpath = os.path.join(os.path.dirname(fistarpath),
                               '{}_quiver_meas_proj_sep.png'.format(pre))
        plot_quiver_meas_proj_sep(df, outpath)

        # not very useful plot given quiver, but leave it in
        outpath = os.path.join(os.path.dirname(fistarpath),
                               '{}_scatter_x_y_sep.png'.format(pre))
        plot_scatter_x_y_sep(df, outpath)


    if (
        float(df['sep'].median()) < qualitycondition['median_px']
        and
        float(df['sep'].std()) < qualitycondition['std_px']
        and
        float(df['sep'].quantile(q=0.9)) < qualitycondition['90th_px']
    ):

        quality_check_result = True

    else:

        quality_check_result = False

    return quality_check_result


def write_wcs_from_spoc(infilename, outfilename=None, observatory='tess',
                        assert_ctype_intact=True):
    """
    For TESS, the SPOC CAL FFIs contain a WCS solution that is typically VERY
    good.  (Like, <0.1 pixel median astrometric residual).

    Given a *.fits `infilename`, with the WCS in the header, this function
    creates a matching *.wcs file with the full header re-written, so that
    almost all WCS-reading utilities can use it.

    If `assert_ctype_intact` is True, it will require that the WCS header has
    the "CTYPE1" key. If this is found to be false, the frame and wcs are moved
    to "badframes".
    """

    if observatory != 'tess':
        raise ValueError('write_wcs_from_spoc assumes TESS FFI headers exist')

    if outfilename is None:
        outfilename = infilename

    wcspath = outfilename.replace('.fits','.wcs')

    if os.path.exists(wcspath):
        print('found {}'.format(wcspath))
        return

    hdulist = fits.open(infilename)

    header = hdulist[0].header

    # Use header to create a new PrimaryHDU and write it to a file.
    hdu = fits.PrimaryHDU(header=header)

    hdu.writeto(wcspath)
    print('made {} from SPOC header'.format(wcspath))

    hdulist.close()

    if assert_ctype_intact:

        hdulist = fits.open(wcspath)
        hdr = hdulist[0].header
        hdulist.close()

        if not ('CTYPE1' in hdr):

            fits_src = outfilename
            fits_dst = os.path.join(
                os.path.dirname(fits_src),
                'badframes',
                os.path.basename(fits_src)
            )

            wcs_src = wcspath
            wcs_dst = os.path.join(
                os.path.dirname(wcs_src),
                'badframes',
                os.path.basename(wcs_src)
            )

            print('WRN: CTYPE1 NOT FOUND IN WCS HEADER. MOVE WCS,IMG,FISTAR\n'
                  '{} -> {}'.format(fits_src,fits_dst))

            shutil.move(fits_src, fits_dst)
            shutil.move(wcs_src, wcs_dst)

            fistar_src = wcs_src.replace('.wcs','.fistar')
            fistar_dst = os.path.join(
                os.path.dirname(fistar_src),
                'badframes',
                os.path.basename(fistar_src)
            )

            if os.path.exists(fistar_src):
                shutil.move(fistar_src, fistar_dst)

        else:
            pass


def sep(x0,y0,x1,y1):
    return np.sqrt((x1-x0)**2 + (y1-y0)**2)


def plot_sep_hist(df, outpath):

    f, ax = plt.subplots(figsize=(4,3))

    weights = np.ones_like(nparr(df['sep']))/float(len(df))

    ax.hist(df['sep'], bins=np.logspace(-2, 1, 19), weights=weights)

    ax.text(
        0.98, 0.98,
        'mean={:.3f}px\nstd={:.3f}px\nmed={:.3f}px\n90th={:.3f}px'.
        format(df['sep'].mean(),df['sep'].std(),
               df['sep'].median(),df['sep'].quantile(q=0.9)),
        va='top', ha='right',
        transform=ax.transAxes
    )
    ax.set_xscale('log')
    ax.set_xlabel('separation [pixels]')
    ax.set_ylabel('relative fraction')
    f.savefig(outpath, bbox_inches='tight', dpi=300)


def plot_scatter_x_y_sep(df, outpath):

    f, ax = plt.subplots(figsize=(4,4))

    norm = mpl.colors.Normalize(vmin=-2.,vmax=1.)

    cax = ax.scatter(df['x_meas'], df['y_meas'], c=np.log10(df['sep']),
                     cmap='viridis', s=1, rasterized=True, linewidths=0,
                     zorder=1, norm=norm)

    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label('log10(measured vs proj separation [px])')

    f.savefig(outpath, bbox_inches='tight', dpi=300)


def plot_quiver_meas_proj_sep(df, outpath):

    f, ax = plt.subplots(figsize=(4,4))

    norm = mpl.colors.Normalize(vmin=-2.,vmax=1.)

    cax = ax.quiver(df['x_meas'], df['y_meas'],
                    df['x_proj']-df['x_meas'],
                    df['y_proj']-df['y_meas'],
                    np.log10(df['sep']),
                    cmap='viridis',
                    norm=norm)

    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')

    cbar = f.colorbar(cax, extend='both')
    cbar.set_label('log10(separation [px])')

    f.savefig(outpath, bbox_inches='tight', dpi=400)


def debug(
    qualitycondition = {'median_px':0.2, '90th_px':0.4, 'std_px': 1.0}
    ):
    ##########################################
    # change these
    ##########

    # # cam1 ccd4
    # reformedfovcatalog = '/home/lbouma/astrometry_tests/GAIADR2-RA86.6153477457052-DEC11.8674772184025-SIZE24.reformed_catalog'
    # fitsfile = '/home/lbouma/astrometry_tests/tess2018350135939-s0006-1-4-0126_cal_img_bkgdsub.fits' #cam1ccd4
    # ra, dec = 86.6153, 11.8674

    # # cam3 ccd3
    # reformedfovcatalog = '/home/lbouma/astrometry_tests/GAIADR2-RA81.7485401006594-DEC-48.4508505192356-SIZE24.reformed_catalog'
    # fitsfile = '/home/lbouma/astrometry_tests/proj1510-s0006-cam3-ccd3-combinedphotref-onenight.fits'
    # ra, dec = 81.7485, -48.4509

    # # cam 2 ccd 1
    # reformedfovcatalog = '/home/lbouma/astrometry_tests/GAIADR2-RA85.1053734708905-DEC-24.6352467680916-SIZE24.reformed_catalog'
    # fitsfile = '/home/lbouma/astrometry_tests/tess2018353175939-s0006-2-1-0126_cal_img_bkgdsub.fits'
    # ra, dec = 85.1054, -24.6352

    # the photometric references are supposedly all shifted to the
    # astrometric reference... so if you want the best photometric WCS, find
    # and use the SPOC astrometric frame WCS.
    # # cam4 ccd 2, photometric reference frame. nb. without tuning pixelorder,
    # your WCS method gets a wonky result!!
    reformedfovcatalog = '/home/lbouma/astrometry_tests/GAIADR2-RA78.554087901003-DEC-59.386568726486-SIZE24.reformed_catalog'
    fitsfile = '/home/lbouma/astrometry_tests/proj1513-s0006-cam4-ccd2-combinedphotref-onenight.fits'
    astromreffile = '/home/lbouma/astrometry_tests/proj1513-camera4-ccd2-astromref-tess2018364145938-s0006-4-2-0126_cal_img_bkgdsub.fits'
    ra, dec = 78.554, -59.387

    use_spoc_wcs = 0

    ##########################################

    ##########################################
    # wont change the following
    ccdextent = {'x':[0.,2048.],'y':[0.,2048.]}
    pixborders = 0.0

    gain = 5.35
    fluxthreshold = 500
    zeropoint = 11.82
    exptime = 1800

    outprojcat = fitsfile.replace('.fits','.projcat')
    fistarpath = fitsfile.replace('.fits','.fistar')
    matchedoutpath = fitsfile.replace('.fits','.matched')
    wcsfile = fitsfile.replace('.fits','.wcs')

    for f in [outprojcat, fistarpath, matchedoutpath, wcsfile]:
        if os.path.exists(f):
            os.remove(f)

    if not os.path.exists(wcsfile) and not use_spoc_wcs:
        astrometrydotnet_solve_frame(fitsfile,
                                     wcsfile,
                                     ra,
                                     dec,
                                     radius=20,
                                     scalelow=20,
                                     scalehigh=22,
                                     tweakorder=6,
                                     useimagenotfistar=True,
                                     pixelerror=0.5,
                                     downsample=2,
                                     uniformize=10)
    else:
        if use_spoc_wcs and 'photref' not in fitsfile:
            write_wcs_from_spoc(fitsfile)
        elif use_spoc_wcs and 'photref' in fitsfile:
            write_wcs_from_spoc(astromreffile, outfilename=fitsfile)


    # given WCS and reformed celestial coord catalog, make the frame projected
    # catalog...
    if not os.path.exists(outprojcat):
        projcatfile = make_frameprojected_catalog(fitsfile,
                                                  reformedfovcatalog,
                                                  ccdextent=ccdextent,
                                                  out=outprojcat,
                                                  wcs=wcsfile,
                                                  removetemp=True,
                                                  pixborders=pixborders)
    # do source extraction...
    if not os.path.exists(fistarpath):
        extract_frame_sources(fitsfile,
                              fistarpath,
                              fistarexec='fistar',
                              ccdextent='0:2048,0:2048',
                              ccdgain=gain,
                              fluxthreshold=fluxthreshold,
                              zeropoint=zeropoint,
                              exptime=exptime)

    quality_check_result = impose_wcs_quality_check(
        fitsfile,
        fistarpath,
        matchedoutpath,
        isspocwcs=use_spoc_wcs,
        projcatfile=projcatfile,
        N_bright=1000,
        N_faint=9000,
        make_plots=True,
        qualitycondition=qualitycondition
    )


if __name__ == "__main__":
    debug()

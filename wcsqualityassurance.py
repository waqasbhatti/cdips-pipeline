"""
Tools to tell you how good your WCS solutions are. In particular, how precisely
can you convert between sky and pixel coordinates? (If worse than your aperture
size, this is a big no-no).
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from numpy import array as nparr

from astropy.io import fits
from astropy import wcs

from aperturephot import (
    make_frameprojected_catalog,
    extract_frame_sources,
    astrometrydotnet_solve_frame
)

def main():

    ##########################################
    # change these
    ##########

    # # cam1 ccd4
    # reformedfovcatalog = '/home/lbouma/GAIADR2-RA86.6153477457052-DEC11.8674772184025-SIZE24.reformed_catalog'
    # fitsfile = '/home/lbouma/tess2018350135939-s0006-1-4-0126_cal_img_bkgdsub.fits' #cam1ccd4
    # ra, dec = 86.6153, 11.8674

    # # cam3 ccd3
    # reformedfovcatalog = '/home/lbouma/GAIADR2-RA81.7485401006594-DEC-48.4508505192356-SIZE24.reformed_catalog'
    # fitsfile = '/home/lbouma/proj1510-s0006-cam3-ccd3-combinedphotref-onenight.fits'
    # ra, dec = 81.7485, -48.4509

    # cam 2 ccd 1
    reformedfovcatalog = '/home/lbouma/GAIADR2-RA85.1053734708905-DEC-24.6352467680916-SIZE24.reformed_catalog'
    fitsfile = '/home/lbouma/tess2018353175939-s0006-2-1-0126_cal_img_bkgdsub.fits'
    ra, dec = 85.1054, -24.6352

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
    matchedpath = fitsfile.replace('.fits','.matched')
    wcsfile = fitsfile.replace('.fits','.wcs')

    for f in [outprojcat, fistarpath, matchedpath, wcsfile]:
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
                                     pixelerror=0.15,
                                     downsample=2,
                                     uniformize=10)
    else:
        write_wcs_from_spoc(fitsfile)


    # given WCS and celestial coord catalog, make the frame projected catalog...
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
    match_and_assess(
        fitsfile,
        use_spoc_wcs,
        incatpath = projcatfile,
        N_bright = 1000,
        N_faint = 9000,
        fistarpath = fistarpath,
        matchedpath = matchedpath
    )


def write_wcs_from_spoc(filename):

    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)

    header = hdulist[0].header

    # # Parse the WCS keywords in the primary HDU
    # w = wcs.WCS(hdulist[0].header)

    # # Now, write out the WCS object as a FITS header
    # header = w.to_header()

    # Use header to create a new PrimaryHDU and write it to a file.
    hdu = fits.PrimaryHDU(header=header)

    # Save to FITS file
    wcspath = filename.replace('.fits','.wcs')
    hdu.writeto(wcspath)
    print('made {} from SPOC header'.format(wcspath))

    hdulist.close()


def match_and_assess(
    fitsfile,
    use_spoc_wcs,
    incatpath = '/home/lbouma/foo.projcat',
    N_bright = 1000, # say top 1000 are too bright
    N_faint = 9000, # go down to get total of 8000 stars
    fistarpath = '/home/lbouma/foo.fistar',
    matchedpath = '/home/lbouma/foo.matched',
):
    """
    You have a WCS for some image.  You use it to project a background catalog
    (e.g, Gaia) to your pixel coordinates.  This gives you the "projected
    catalog", incatpath.

    You then measure a bunch of star positions with fistar, to get fistarpath.

    This function matches the projected positions to the measured positions
    with the symmetric point-matching algorithm implemented in `grmatch`.

    It then makes three plots to assess the quality of the WCS solution:
        * a histogram of the separation distribution
        * a scatter plot with the separations colored
        * a vector plot with the position differences as arrows and colored
    """

    df = pd.read_csv(
        incatpath,
        delim_whitespace=True,
        names='id,ra,dec,xi,eta,G,Rp,Bp,plx,pmra,pmdec,varflag,x,y'.split(',')
    )

    outdf = df.sort_values(by='Rp').iloc[N_bright:N_faint]
    cutprojcat = os.path.join(
        os.path.dirname(incatpath),
        'cut_{}'.format(os.path.basename(incatpath))
    )
    outdf.to_csv(cutprojcat, index=False, sep=' ', header=False)

    with open(fistarpath, 'r') as f:
        fistarlines = f.readlines()

    fistaroutlines = fistarlines[6+N_bright:6+N_faint] # 6 comment lines
    cutfistar = os.path.join(
        os.path.dirname(fistarpath),
        'cut_{}'.format(os.path.basename(fistarpath))
    )

    with open(cutfistar, 'w') as f:
        f.writelines(fistaroutlines)

    matchcmd = (
        'grmatch --match-points -i {cutfistar} --col-inp 2,3 '
        '-r {cutprojcat} --col-ref 13,14 --output {matched}'.
        format(cutfistar=cutfistar, cutprojcat=cutprojcat, matched=matchedpath)
    )

    rc = os.system(matchcmd)
    print(matchcmd)

    if rc != 0:
      raise RuntimeError('grmatch failed somehow')

    df = pd.read_csv(
        matchedpath,
        delim_whitespace=True,
        names=('id,ra,dec,xi,eta,G,Rp,Bp,plx,pmra,pmdec,varflag,x_proj,y_proj,'+
               'Ident,x_meas,y_meas,Bg,Amp,S,D,K,Flux,S/N').split(',')
    )

    if not use_spoc_wcs:
        # NOTE: means I somehow messed up the coordinates somewhere?!
        df['x_proj'] -= 0.5
        df['y_proj'] -= 0.5

    else:
        # SPOC rows/cols were trimmed, which shifts WCS
        hdulist = fits.open(fitsfile)
        hdr = hdulist[0].header

        df['x_proj'] -= (hdr['SCCSA']-1)
        df['y_proj'] -= (hdr['SCIROWS']-1)
        df['x_proj'] -= 0.5
        df['y_proj'] -= 0.5

    df['sep'] = sep(
        df['x_meas'],
        df['y_meas'],
        df['x_proj'],
        df['y_proj']
    )

    pre = os.path.splitext(os.path.basename(fistarpath))[0]
    if use_spoc_wcs:
        pre = pre+'_spocwcs'
    outpath = os.path.join(os.path.dirname(fistarpath),
                           '{}_sep_hist.png'.format(pre))
    plot_sep_hist(df, outpath)

    outpath = os.path.join(os.path.dirname(fistarpath),
                           '{}_scatter_x_y_sep.png'.format(pre))
    plot_scatter_x_y_sep(df, outpath)

    outpath = os.path.join(os.path.dirname(fistarpath),
                           '{}_quiver_meas_proj_sep.png'.format(pre))
    plot_quiver_meas_proj_sep(df, outpath)


def sep(x0,y0,x1,y1):
    return np.sqrt((x1-x0)**2 + (y1-y0)**2)

def plot_sep_hist(df, outpath):

    f, ax = plt.subplots(figsize=(4,3))

    weights = np.ones_like(nparr(df['sep']))/float(len(df))

    ax.hist(df['sep'], bins=np.logspace(-2, 1, 19), weights=weights)

    ax.text(
        0.02, 0.98,
        'mean={:.3f}px\nstd={:.3f}px\nmed={:.3f}px'.
        format(df['sep'].mean(),df['sep'].std(),df['sep'].median()),
        va='top', ha='left',
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


if __name__ == "__main__":
    main()

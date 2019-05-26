"""
Tools to tell you how good your WCS solutions are. In particular, how precisely
can you convert between sky and pixel coordinates? (If worse than your aperture
size, this is a big no-no).
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import os

from aperturephot import make_frameprojected_catalog

def main():

    #FIXME FIXME: automate the process of generating the .projcat and the
    #.fistar files. b/c we'll need to tweak the WCS solution, most likely!!

    # and/or 

    if overwrite or not os.path.exists(outprojcat):
        #FIXME: you can pass the wcs here, to make the projcat files from the
        # TESS headers!!
        projcatfile = make_frameprojected_catalog(fits,
                                                  reformedfovcatalog,
                                                  ccdextent=ccdextent,
                                                  out=outprojcat,
                                                  removetemp=removesourcetemp,
                                                  pixborders=pixborders)

    match_and_assess(
        incatpath = '/home/lbouma/foo.projcat',
        N_bright = 1000,
        N_faint = 9000,
        fistarpath = '/home/lbouma/foo.fistar',
        matchedpath = '/home/lbouma/foo.matched'
    )


def match_and_assess(
    incatpath = '/home/lbouma/foo.projcat',
    N_bright = 1000, # say top 1000 are too bright
    N_faint = 9000, # go down to get total of 8000 stars
    fistarpath = '/home/lbouma/foo.fistar',
    matchedpath = '/home/lbouma/foo.matched'
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

    if rc != 0:
      raise RuntimeError('grmatch failed somehow')

    df = pd.read_csv(
        matchedpath,
        delim_whitespace=True,
        names=('id,ra,dec,xi,eta,G,Rp,Bp,plx,pmra,pmdec,varflag,x_proj,y_proj,'+
               'Ident,x_meas,y_meas,Bg,Amp,S,D,K,Flux,S/N').split(',')
    )

    df['sep'] = sep(
        df['x_meas'],
        df['y_meas'],
        df['x_proj'],
        df['y_proj']
    )

    outpath = os.path.join(os.path.dirname(fistarpath), 'foo_sep_hist.png')
    plot_sep_hist(df, outpath)

    outpath = os.path.join(os.path.dirname(fistarpath),
                           'foo_scatter_x_y_sep.png')
    plot_scatter_x_y_sep(df, outpath)

    outpath = os.path.join(os.path.dirname(fistarpath),
                           'foo_quiver_meas_proj_sep.png')
    plot_quiver_meas_proj_sep(df, outpath)


def sep(x0,y0,x1,y1):
    return np.sqrt((x1-x0)**2 + (y1-y0)**2)

def plot_sep_hist(df, outpath):

    f, ax = plt.subplots(figsize=(4,3))
    ax.hist(df['sep'], bins=np.logspace(-2, 1, 10))
    ax.text(
        0.02, 0.98,
        '$\mu$={:.3f}px\n$\sigma$={:.3f}px'.
        format(df['sep'].mean(),df['sep'].std()),
        va='top', ha='left',
        transform=ax.transAxes
    )
    ax.set_xscale('log')
    ax.set_xlabel('separation [pixels]')
    ax.set_ylabel('relative fraction')
    f.savefig(outpath, bbox_inches='tight', dpi=300)

def plot_scatter_x_y_sep(df, outpath):

    f, ax = plt.subplots(figsize=(4,4))

    cax = ax.scatter(df['x_meas'], df['y_meas'], c=np.log10(df['sep']),
                     cmap='viridis', s=1, rasterized=True, linewidths=0,
                     zorder=1)

    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')

    cbar = f.colorbar(cax)
    cbar.set_label('log10(measured vs proj separation [px])')

    f.savefig(outpath, bbox_inches='tight', dpi=300)


def plot_quiver_meas_proj_sep(df, outpath):

    f, ax = plt.subplots(figsize=(4,4))

    cax = ax.quiver(df['x_meas'], df['y_meas'],
                    df['x_proj']-df['x_meas'],
                    df['y_proj']-df['y_meas'],
                    np.log10(df['sep']),
                    cmap='viridis')

    ax.set_xlabel('x [px]')
    ax.set_ylabel('y [px]')

    cbar = f.colorbar(cax)
    cbar.set_label('log10(separation [px])')

    f.savefig(outpath, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()

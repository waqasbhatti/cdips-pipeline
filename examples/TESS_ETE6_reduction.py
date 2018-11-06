# -*- coding: utf-8 -*-
'''
usage: TESS_ETE6_reduction.py [-h] [--fitsdir FITSDIR] [--fitsglob FITSGLOB]
                              [--projectid PROJECTID] [--field FIELD]
                              [--outdir OUTDIR] [--lcdirectory LCDIRECTORY]
                              [--aperturelist APERTURELIST]
                              [--kernelspec KERNELSPEC]
                              [--photdisjointradius PHOTDISJOINTRADIUS]
                              [--anetfluxthreshold ANETFLUXTHRESHOLD]
                              [--anettweak ANETTWEAK]
                              [--anetradius ANETRADIUS]
                              [--initccdextent INITCCDEXTENT]
                              [--epdsmooth EPDSMOOTH]
                              [--epdsigclip EPDSIGCLIP] [--nworkers NWORKERS]
                              [--convert_to_fitsh_compatible]
                              [--no-convert_to_fitsh_compatible]
                              [--tuneparameters TUNEPARAMETERS]

Given images, make lightcurves.

optional arguments:
  -h, --help            show this help message and exit
  --fitsdir FITSDIR     e.g., /foo/FFI/RED
  --fitsglob FITSGLOB   e.g., *.fits.gz. Chosen frames are fitsdir+fitsglob.
  --projectid PROJECTID
                        unique integer identifying this reduction run
  --field FIELD         field identifier string, e.g., "sector1"
  --outdir OUTDIR       e.g., /foo/FFI/RED
  --lcdirectory LCDIRECTORY
                        e.g., /foo/RED/LC
  --aperturelist APERTURELIST
                        formatted as "inner radius: outer radius: bkgd annulus
                        width"
  --kernelspec KERNELSPEC
                        https://fitsh.net/wiki/man/ficonv explains the
                        formatting.i/<spatial order> identity kernel (a.k.a.
                        flux term) with thespecified order of polynomial
                        spatial variations.b/<spatial order> constant offset
                        kernel (a.k.a. backgroundterm) with the specified
                        order of polynomial spatialvariationsd=<size>/<spatial
                        order> discrete kernel with the half-size of<size> and
                        the specified order of polynomial spatial
                        variationsg=<size>,<sigma>,<order>/<spatial order>
                        Gaussian kernel withthe half-size of <size>, standard
                        deviation of <sigma> andHermite basis order of
                        <order>, with the specified order ofpolynomial spatial
                        variations
  --photdisjointradius PHOTDISJOINTRADIUS
                        https://fitsh.net/wiki/man/fiphot gives details.
                        During the bacground determination on the aperture
                        annuli, omit the pixels which are closer to the other
                        centroids than the specified radius.
  --anetfluxthreshold ANETFLUXTHRESHOLD
                        flux threshold used to identify bright stars in anet
                        initial wcs solution attempt of frames.
  --anettweak ANETTWEAK
                        SIP (Simple Imaging Polynomial) order. Higher = more
                        flexibility for frame distortion. Max of 10.
  --anetradius ANETRADIUS
                        Radius (in deg) over which to search for index files
                        in astrometry.net or net. For TESS, you want to go big
                        -- e.g., 30.
  --initccdextent INITCCDEXTENT
                        section <x1>:<x2>,<y1>:<y2> for initial source
                        extraction and astrometry. (not for initial
                        photometry).
  --epdsmooth EPDSMOOTH
                        number of cadences used when passing a median filter
                        over themagnitude time series before EPD (see
                        imagesubphot.epd_magseries_imagesub). For HATPI
                        default is 21, or a ~10 minute median filter. For TESS
                        FFIs, try 11 = 5.5hrs.
  --epdsigclip EPDSIGCLIP
                        sigma clipping aplpied to epd lightcurves. assumed
                        symmetric.
  --nworkers NWORKERS   number of workers to thread over
  --convert_to_fitsh_compatible
                        converts calibrated frames to single-fits-header fitsh
                        readableobjects. if already done, will find them and
                        skip.
  --no-convert_to_fitsh_compatible
                        don't do fitsh-compatibility conversion.
  --tuneparameters TUNEPARAMETERS
                        TUNING: iterate through different img subtraction
                        parameters if true. Gain speed by selecting a small
                        set of frames. Otherwise, run in FULL reduction mode.
'''
from __future__ import division, print_function

import os, time
import matplotlib as mpl
mpl.use('AGG')
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import aperturephot as ap, shared_variables as sv, autoimagesub as ais, \
       imagesubphot as ism, tessutils as tu, lcstatistics as lcs, \
       imageutils as iu
from glob import glob
from tqdm import tqdm
from astropy.io import fits
from astropy import units as units, constants as constants
from datetime import datetime

import argparse

np.random.seed(42)


def get_files_needed_before_image_subtraction(
        fitsdir, fitsglob, outdir, initccdextent, ccdgain, zeropoint, exptime,
        ra_nom, dec_nom,
        catra, catdec, ccd_fov,
        catalog, catalog_file, reformed_cat_file,
        fnamestr='*-1-1-0016_cal_img.fits', anetfluxthreshold=20000,
        fistarglob='*.fistar',
        width=13, anettweak=6, anetradius=30, xpix=2048, ypix=2048, cols=(2,3),
        brightrmag=6.0, faintrmag=13.0,
        fistarfluxthreshold=1000,
        aperturelist='1.45:7.0:6.0,1.95:7.0:6.0,2.45:7.0:6.0',
        nworkers=20,
        useastrometrydotnet=True,
        useimagenotfistar=True,
        extractsources=True
    ):
    '''
    get .fistar, .fiphot, and .wcs files needed before image subtraction

    1. run parallel_extract_sources on all frames with threshold ~ 10000 to get
    bright stars for astrometry.
    2. run parallel_anet to get precise WCS headers for all frames.
    3. run make_fov_catalog to get a FOV source catalog for the field.
    4. run reform_fov_catalog to cut this down to the columns needed for magfit
    only.
    5. run parallel_fitsdir_photometry for photometry on all frames (via
    fistar)
    '''

    ap.parallel_extract_sources(fitsdir, outdir, ccdextent=initccdextent,
                                ccdgain=ccdgain,
                                fluxthreshold=anetfluxthreshold,
                                zeropoint=zeropoint, exptime=exptime,
                                tailstr='.fits',
                                fnamestr=fnamestr)

    if useastrometrydotnet:

        ap.fistardir_to_xy(fitsdir, fistarglob=fistarglob)

        ap.parallel_astrometrydotnet(
            fitsdir, outdir, ra_nom, dec_nom,
            fistarfitsxyglob=fistarglob.replace('.fistar','.fistar-fits-xy'),
            tweakorder=anettweak, radius=anetradius, xpix=xpix, ypix=ypix,
            nworkers=nworkers, scalelow=10, scalehigh=30,
            scaleunits='arcsecperpix', nobjs=200, xcolname='ximage',
            ycolname='yimage', useimagenotfistar=useimagenotfistar,
            downsample=4
        )

    else:
        ap.parallel_anet(fitsdir, outdir, ra_nom, dec_nom,
                         fistarglob=fistarglob,
                         infofromframe=False, width=width, tweak=anettweak,
                         radius=anetradius,
                         xpix=xpix, ypix=ypix,
                         cols=cols # columns with x,y in fistar file.
        )

    # This function gets all the sources in the field of view of the frame, given
    # its central pointing coordinates and plate-scale from 2MASS. This catalog
    # file is then be used as input to make_source_list below.
    _ = ap.make_fov_catalog(ra=catra, dec=catdec, size=ccd_fov,
                            brightrmag=brightrmag, faintrmag=faintrmag,
                            fits=None, outfile=None, outdir=outdir,
                            catalog=catalog, catalogpath=None,
                            columns=None, observatory='tess')

    # This converts the full output catalog from 2massread, etc. to the format
    # required for magfit. Also useful for general reforming of the columns.
    ap.reform_fov_catalog(catalog_file, reformed_cat_file)

    if extractsources==True:
        fiphot_xycols = '7,8'
    else:
        fiphot_xycols = '13,14'
    ap.parallel_fitsdir_photometry(fitsdir, outdir, reformed_cat_file,
                                   fluxthreshold=fistarfluxthreshold,
                                   ccdextent={'x':[0.,2048.],'y':[0.,2048.]},
                                   pixborders=0.0,
                                   aperturelist=aperturelist,
                                   removesourcetemp=True,
                                   removesourcelist=False, binaryoutput=False,
                                   nworkers=nworkers, maxtasksperworker=1000,
                                   saveresults=True, rejectbadframes=True,
                                   minsrcbgv=200.0, maxmadbgv=150.0,
                                   maxframebgv=2000.0, minnstars=500,
                                   fitsglob=fitsglob, ccdgain=ccdgain,
                                   ccdexptime=exptime, zeropoint=zeropoint,
                                   extractsources=extractsources,
                                   fovcat_xycols=(12,13),
                                   projcat_xycols=(24,25),
                                   fiphot_xycols=fiphot_xycols,
                                   observatory='tess'
                                  )


def initial_wcs_worked_well_enough(outdir, fitsglob):
    # check whether initial anet converged on over half of frames.

    N_fitsfiles = len(glob(outdir+fitsglob))
    N_wcsfiles = len(glob(outdir+fitsglob.replace('.fits','.wcs')))

    if N_wcsfiles < N_fitsfiles/2:
        return False
    else:
        return True


def is_presubtraction_complete(outdir, fitsglob, lcdir, percentage_required=95,
                               extractsources=False):
    '''
    require at least e.g., 95% of the initial astrometry, photometry, etc to
    exist to return True. in that case, or if any stats_files products are
    found, move on to image subtraction.  else, returns False.
    '''

    N_fitsfiles = len(glob(outdir+fitsglob))
    N_fistarfiles = len(glob(outdir+fitsglob.replace('.fits','.fistar')))
    N_wcsfiles = len(glob(outdir+fitsglob.replace('.fits','.wcs')))
    N_projcatalog = len(glob(outdir+fitsglob.replace('.fits','.projcatalog')))
    N_sourcelistfiles = len(glob(outdir+fitsglob.replace('.fits','.sourcelist')))
    N_fiphotfiles = len(glob(outdir+fitsglob.replace('.fits','.fiphot')))

    if extractsources:
        N_files = [N_fitsfiles, N_fistarfiles, N_wcsfiles, N_fiphotfiles,
                   N_sourcelistfiles, N_projcatalog]
    else:
        N_files = [N_fitsfiles, N_fistarfiles, N_wcsfiles, N_fiphotfiles,
                   N_projcatalog]

    statsdir = lcdir+'stats_files'
    N_statsfiles_products = len(glob(statsdir+"*"))
    if N_statsfiles_products >= 1:
        print('found stats_files product. skipping to detrending.')
        return True

    elif (np.any( np.abs(np.diff(N_files)/N_fitsfiles)>(1-percentage_required/100) )
        or
        np.any( np.array(N_files)==0 )
    ):
        print('did not find {:d}% completed source extraction, astrometry, '.
              format(percentage_required)+
              'and initial photometry')
        return False

    else:
        print('found {:d}% completed source extraction, astrometry, '.
              format(percentage_required)+
              'and initial photometry. skipping presubtraction steps.')
        return True


def is_imagesubtraction_complete(fitsdir, fitsglob, lcdir):

    N_subfitslist = len(glob(
        fitsdir+'rsub-????????-'+fitsglob.replace('.fits','-xtrns.fits')))
    N_iphotlist = len(glob(
        fitsdir+'rsub-????????-'+fitsglob.replace('.fits','.iphot')))
    N_kernellist = len(glob(
        fitsdir+'rsub-????????-'+
        fitsglob.replace('.fits','-xtrns.fits-kernel')))
    N_subconvjpglist = len(glob(
        fitsdir+'JPEG-SUBTRACTEDCONV-rsub-*-tess*.jpg'))
    N_lcs = len(glob(
        lcdir+'*.grcollectilc'))

    N_files = [N_subfitslist, N_iphotlist, N_kernellist, N_subconvjpglist]
    full_N_files = N_files.append(N_lcs)

    statsdir = lcdir+'stats_files'
    N_statsfiles_products = len(glob(statsdir+"*"))
    if N_statsfiles_products >= 1:
        print('found stats_files product. skipping to detrending.')
        return True

    elif np.any( np.diff(N_files) ) or np.any( np.array(full_N_files)==0 ):
        print('did not find completed image-subtracted photometry products '
              'including photometry, images, and lightcurves.')
        return False

    else:
        return True


def run_imagesubtraction(fitsdir, fitsglob, fieldinfo, photparams, fits_list,
                         photreftype, dbtype, reformed_cat_file, xtrnsglob,
                         iphotpattern, lcdirectory, kernelspec='b/4;i/4;d=4/4',
                         refdir=sv.REFBASEDIR, nworkers=1,
                         aperturelist='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                         photdisjointradius=2, colorscheme='bwr',
                         photreffluxthreshold=30000, extractsources=True):

    ccdgain = photparams['ccdgain']
    exptime = photparams['ccdexptime']
    zeropoint = photparams['zeropoint']

    # Step ISP0.
    _ = ais.parallel_frames_to_database(fitsdir, 'calibratedframes',
                                        observatory='tess', fitsglob=fitsglob,
                                        overwrite=False,
                                        badframetag='badframes',
                                        nonwcsframes_are_ok=False,
                                        nworkers=nworkers, maxworkertasks=1000)

    # Step ISP1.
    arefinfo = ais.dbgen_get_astromref(fieldinfo, makeactive=True,
                                       observatory='tess', overwrite=False,
                                       refdir=sv.REFBASEDIR, database=None)
    if arefinfo == None:
        raise AssertionError('you need an astrometric reference!!')

    # Step ISP2. Takes ~1 sec per 20-30 frames.
    _ = ism.get_smoothed_xysdk_coeffs(fitsdir, fistarglob='*.fistar',
                                      nworkers=nworkers, maxworkertasks=1000)

    # Step ISP3. ~600 in 5 minutes --> 2 frames per second, running over 20 workers.
    # If called with warpcheck=True and warpthreshold=2000, many things will be
    # moved to badframes that are not in fact badframes. The warpthreshold is a bad
    # empirical thing that should be avoided.
    _ = ais.framelist_make_xtrnsfits(fits_list, fitsdir, fitsglob, outdir=None,
                                     refinfo='foobar', warpcheck=False,
                                     warpthreshold=15000.0, warpmargins=100,
                                     nworkers=nworkers, observatory='tess',
                                     maxworkertasks=1000, fieldinfo=fieldinfo)

    # # Optional Step ISP3.5: parallelized move of astrometry ref shifted frames to database
    # out = ais.parallel_frames_to_database(fitsdir, 'arefshifted_frames',
    #                                       fitsglob='1-???????_?-xtrns.fits',
    #                                       network='HP', overwrite=False,
    #                                       badframetag='badframes',
    #                                       nonwcsframes_are_ok=False, nworkers=nworkers,
    #                                       maxworkertasks=1000)
    # 

    # Step ISP4: the next thing to do is to select a bunch of frames that can serve
    # as photometric reference frames (photrefs).

    xtrnsfiles = glob(fitsdir+xtrnsglob)
    photrefinfo = ais.generate_photref_candidates_from_xtrns(
        xtrnsfiles, minframes=50, observatory='tess',
        maxbackgroundstdevpctile=100., maxbackgroundmedianpctile=70.,
        minngoodobjectpctile=70., forcecollectinfo=False, nworkers=nworkers,
        maxworkertasks=1000)

    # Optional Step ISP5: amend the list, if needed.
    # photrefinfo = ais.amend_candidate_photrefs(photrefinfo)

    # Step ISP6: make a photometric reference frame
    _ = ais.generate_combined_photref(
        photrefinfo, photreftype, dbtype,
        photref_reformedfovcat=reformed_cat_file, makeactive=True, field=None,
        ccd=None, projectid=None, combinemethod='median',
        kernelspec=kernelspec, ccdgain=ccdgain, zeropoint=zeropoint,
        ccdexptime=exptime, extractsources=extractsources,
        apertures=aperturelist, framewidth=None, searchradius=8.0,
        nworkers=nworkers, maxworkertasks=1000, observatory='tess',
        fieldinfo=fieldinfo, photreffluxthreshold=photreffluxthreshold)

    # Step ISP7: convolve and subtract all FITS files in the xtrnsfits list from the
    # photometric reference.  With 30 workers, at best process ~few frames per
    # second.

    _ = ais.parallel_xtrnsfits_convsub(
        xtrnsfiles, photreftype, fitsdir=fitsdir, fitsglob=fitsglob,
        outdir=None, observatory='tess', fieldinfo=fieldinfo,
        reversesubtract=True, kernelspec=kernelspec, nworkers=nworkers,
        maxworkertasks=1000, colorscheme=colorscheme)

    # Step ISP8: do photometry on your subtracted frames to produce .iphot files.
    # With 30 workers, at best process ~few frames per second.

    subfitslist = glob(fitsdir+'rsub-????????-'+
                       fitsglob.replace('.fits','-xtrns.fits'))
    _ = ais.parallel_convsubfits_staticphot(
        subfitslist, fitsdir=fitsdir, fitsglob=fitsglob,
        photreftype=photreftype, kernelspec=kernelspec,
        lcapertures=aperturelist, photdisjointradius=photdisjointradius,
        outdir=None, fieldinfo=fieldinfo, observatory='tess',
        nworkers=nworkers, maxworkertasks=1000, photparams=photparams)

    # Step ISP9 + 10 : dump lightcurves.
    ism.dump_lightcurves_with_grcollect(
        iphotpattern, lcdirectory, '4g', lcextension='grcollectilc',
        objectidcol=3, observatory='tess')

    # # Alternative Step 9 + 10: add the .iphot files to postgres database. Collect
    # # LCs from postgres, and then make difference image lightcurve (*.ilc) files in
    # # lcdirectory. AKA "light curve dump", or "the transposition problem".
    # # Surprisingly computationally expensive.  An alternative approach is
    # # fitsh's `grcollect`, which skips the database architecture entirely.
    # 
    # print('beginning insert_phots_into_database')
    # ais.insert_phots_into_database(sv.REDPATH, frameglob='rsub-*-xtrns.fits',
    #                                photdir=None, photglob='rsub-*-%s.iphot',
    #                                maxframes=None, overwrite=False, database=None)
    # 
    # hatidlist = ais.get_hatidlist_from_cmrawphot(projectid, field, ccd, photreftype)
    # 
    # print('beginning lightcurve dump')
    # ais.parallel_dbphot_lightcurves_hatidlist(hatidlist, lcdirectory)


def run_detrending(epdstatfile, tfastatfile, lcdirectory, epdlcglob,
                   reformed_cat_file, statsdir, field, epdsmooth=11,
                   epdsigclip=10, nworkers=10, binlightcurves=False):
    '''
    Step ISP11: do EPD on all the LCs, and collect stats on the results.
    for ISP LCs, use lcmagcols=([27,28,29],[30,],[30,],[30,])

    Step ISP12: do TFA on all the LCs. First, choose TFA template stars using the
    .epdlc stats. Then run TFA, to get .tfalc.TF{1,2,3} files. Turn them into
    single .tfalc files. Then collect statistics.
    '''

    if not os.path.exists(epdstatfile):

        _ = ism.parallel_run_epd_imagesub(lcdirectory,
                                          ilcglob='*.grcollectilc',
                                          outdir=None, smooth=epdsmooth,
                                          sigmaclip=epdsigclip, nworkers=nworkers,
                                          maxworkertasks=1000, minndet=100)

        ap.parallel_lc_statistics(lcdirectory, epdlcglob,
                                  reformed_cat_file, tfalcrequired=False,
                                  fovcatcols=(0,9), # objectid, magcol to use
                                  fovcatmaglabel='r', outfile=epdstatfile,
                                  nworkers=nworkers,
                                  workerntasks=500, rmcols=[14,19,24],
                                  epcols=[27,28,29], tfcols=[30,31,32],
                                  rfcols=None, correctioncoeffs=None,
                                  sigclip=5.0)
    else:
        print('already made EPD LC stats file')

    epdmadplot = glob(os.path.join(statsdir, '*median-EP1-vs-mad-*png'))
    if not epdmadplot:
        ap.plot_stats_file(epdstatfile, statsdir, field, binned=False,
                           logy=True, logx=False, correctmagsafter=None,
                           rangex=(5.9,16), observatory='tess')
    else:
        print('already made EPD LC plots')

    if not os.path.exists(lcdirectory+'aperture-1-tfa-template.list'):
        _ = ap.choose_tfa_template(epdstatfile, reformed_cat_file, lcdirectory,
                                   ignoretfamin=False, fovcat_idcol=0,
                                   fovcat_xicol=3, fovcat_etacol=4,
                                   fovcat_magcol=9, min_ndet=100,
                                   min_nstars=50, max_nstars=1000,
                                   brightest_mag=8.5, faintest_mag=13.0,
                                   max_rms=0.1, max_sigma_above_rmscurve=5.0,
                                   outprefix=statsdir, tfastage1=True)
    if not os.path.exists(tfastatfile):
        templatefiles = glob(lcdirectory+'aperture-?-tfa-template.list')
        ism.parallel_run_tfa(lcdirectory, templatefiles, epdlc_glob='*.epdlc',
                            epdlc_jdcol=0, epdlc_magcol=(27,28,29),
                            template_sigclip=5.0, epdlc_sigclip=5.0, nworkers=nworkers,
                            workerntasks=1000)
        ap.parallel_lc_statistics(lcdirectory, '*.epdlc', reformed_cat_file,
                                  tfalcrequired=True,
                                  fovcatcols=(0,9), # objectid, magcol from fovcat
                                  fovcatmaglabel='r',
                                  outfile=tfastatfile,
                                  nworkers=nworkers,
                                  workerntasks=500,
                                  rmcols=[14,19,24],
                                  epcols=[27,28,29],
                                  tfcols=[30,31,32],
                                  rfcols=None,
                                  correctioncoeffs=None,
                                  sigclip=5.0)
    else:
        print('already made TFA LC stats file')

    tfaboolstatusfile = os.path.join(statsdir,'are_tfa_plots_done.txt')
    if not os.path.exists(tfaboolstatusfile):
        ap.plot_stats_file(tfastatfile, statsdir, field, binned=False,
                           logy=True, logx=False, correctmagsafter=None,
                           rangex=(5.9,16), observatory='tess')
        with open(tfaboolstatusfile+'','w') as f:
            f.write('1\n')
    else:
        print('found done TFA plots (unbinned) through {:s}. continuing.'.
              format(tfaboolstatusfile))

    if binlightcurves:

        binsizes = [3600,21600]

        binnedlcfiles = glob(lcdirectory+'*.binned-*sec-lc.pkl')
        if len(binnedlcfiles) == 0:
            ap.parallel_bin_lightcurves(
                lcdirectory, '*epdlc', binsizes=binsizes,
                lcexts=('epdlc', 'tfalc'),
                lcmagcols=([27,28,29],[30,31,32]),
                jdcol=0, nworkers=nworkers, workerntasks=1000)
        else:
            print('found >=1 binned lightcurve. continuing.')

        onehr_binstatfile = (
            os.path.join(statsdir,'onehr-binned-lightcurve-statistics.txt'))
        onehrglob = '*.binned-3600sec-lc.pkl'
        sixhr_binstatfile = (
            os.path.join(statsdir,'sixhr-binned-lightcurve-statistics.txt'))
        sixhrglob = '*.binned-21600sec-lc.pkl'

        # make statfiles and MAD vs. mag plots for all binned LCs.
        for binstatfile, binglob, cadence in zip(
            [onehr_binstatfile, sixhr_binstatfile], [onehrglob, sixhrglob],
            binsizes):

            if not os.path.exists(binstatfile):
                ap.parallel_binnedlc_statistics(
                    lcdirectory, binglob, reformed_cat_file, fovcatcols=(0,9),
                    fovcatmaglabel='r', corrmagsource=None, corrmag_idcol=0,
                    corrmag_magcols=[122,123,124], outfile=binstatfile,
                    nworkers=nworkers, workerntasks=500, sigclip=5)
            else:
                print('found {:s}, continue'.format(binstatfile))

            outprefix = field+'-'+str(cadence)
            ap.plot_stats_file(binstatfile, statsdir, outprefix,
                               binned=cadence, logy=True, logx=False,
                               correctmagsafter=None, rangex=(5.9,16),
                               observatory='tess')
    else:
        print('will not bin lightcurves or make assoiated statplots')


def run_detrending_on_raw_photometry():
    '''
    run detrending on raw photometry, to compare vs. image subtracted

    Steps [1-5] of raw photometry (source extraction, astrometry, creation of a
    source catalog, reforming the source catalog, and raw photometry) have been
    performed.  To detrend, we first "magnitude fit" (see Sec 5.5 of Zhang,
    Bakos et al 2016).  This procedure empirically fits magnitude differences
    (vs. a reference) for each star in each frame, using on-camera position,
    catalog magnitude, catalog color, and subpixel position. Then we run EPD
    and TFA. Procedurally:

    6. run get_magfit_frames to select a single magfit photometry reference and
    set up per-CCD work directories, symlinks, etc. for the next steps.

    7. run make_magfit_config to generate magfit config files for
    MagnitudeFitting.py

    8. run make_fiphot_list to make lists of fiphot files for each CCD.

    9. run MagnitudeFitting.py in single reference mode.

    10. run do_masterphotref.py to get the master mag fit reference.

    11. run MagnitudeFitting.py in master reference mode.

    12. run parallel_collect_lightcurves to collect all lightcurves into .rlc
    files.

    13. run serial_run_epd or parallel_run_epd to do EPD on all LCs.

    14. run parallel_lc_statistics to collect stats on .epdlc files.

    15. run choose_tfa_template to choose TFA template stars using the .epdlc
    stats.

    16. run parallel_run_tfa for TFA to get .tfalc files

    17. run parallel_lc_statistics to collect stats on .tfalc files.

    20. run plot_stats_file to make MAD vs. mag plots for all unbinned and
    binned LCs.
    '''

    raise NotImplementedError

    mfdict = ap.get_magfit_frames(fitsdir, sv.LOCAL_GLOBPATTERN, fitsdir,
                                  selectreference=True, linkfiles=False,
                                  framestats=False, observatory='tess')

    # note: magnitude fitting is probably not needed for space data.
    # what should our comparison case be here?

    # ap.make_magfit_config()

    # ap.make_fiphot_list()

    # ap.run_magfit()

    # ap.get_master_photref()

    # ap.run_magfit()

    # ap.parallel_collect_aperturephot_lightcurves(fitsdir, lcdirectory,
    #                                              photext='fiphot',
    #                                              skipcollectedlcs=False,
    #                                              nworkers=nworkers)

    # ap.parallel_run_epd()

    # ap.parallel_lc_statistics()

    # ap.choose_tfa_template()

    # ap.parallel_run_tfa()

    # ap.parallel_lc_statistics()

    # ap.plot_stats_file()


def assess_run(statsdir, lcdirectory, starttime, outprefix, fitsdir,
               projectid, field, camera, ccd,
               binned=False, make_whisker_plot=True, whisker_xlim=[4,17],
               whisker_ylim=[1e-5,1e-1]):
    '''
    write files with summary statistics of run.

    args:
        statsdir (str): e.g.,
        '/nfs/phtess1/ar1/TESS/SIMFFI/LC/TUNE/orbit-10/ISP_1-2/stats_files/'

        lcdirectory (str): one level up from statsdir.

        starttime (datetime obj)

    kwargs:
        is statfile binned?

        make_whisker_plot (bool)
    '''

    # how long did the pipeline take?
    endtime = datetime.utcnow()
    howlong = (endtime - starttime).total_seconds()*units.s

    SUMMARY={
        'starttime': starttime.isoformat(),
        'endttime': endtime.isoformat(),
        'howlong_day':howlong.to(units.day).value,
        'howlong_hr':howlong.to(units.hr).value,
        'howlong_min':howlong.to(units.minute).value,
        'howlong_sec':howlong.to(units.second).value
    }
    summarydf = pd.DataFrame(SUMMARY, index=[0])
    summarypath = statsdir+'timing_summary.csv'
    summarydf.to_csv(summarypath, index=False)

    whiskerfiles = glob(statsdir+'whisker_*png')
    if not whiskerfiles:
        lcs.whisker_MAD_stats_and_plots(statsdir, outprefix, binned=binned,
                                        make_whisker_plot=make_whisker_plot,
                                        whisker_xlim=whisker_xlim,
                                        whisker_ylim=whisker_ylim)
    else:
        print('found whisker plots')

    is_image_noise_gaussian(fitsdir, projectid, field, camera, ccd)


def _plot_normalized_subtractedimg_histogram(
    subimg_normalized, rsubimgfile, zerooutnans=True):
    '''subimg_normalized: subtracted, normalized image.'''

    if zerooutnans:
        subimg_normalized[np.isnan(subimg_normalized)] = 0
        print('ERR! you maybe should not have nans in your image...')
    else:
        raise NotImplementedError

    plt.close('all')
    f, ax = plt.subplots(figsize=(4,3))

    bins = np.arange(-100,100,0.1)
    bins, edges = np.histogram(subimg_normalized, bins=bins)
    left, right = edges[:-1], edges[1:]

    # data
    X = np.array([left,right]).T.flatten()
    Y = np.array([bins,bins]).T.flatten()

    ax.plot(X, Y, label='data')

    # theory: a gaussian
    # 1/sqrt(2pi*sigma^2) * exp( -(x-mu)^2 / (2 sigma^2)  )
    mu, sigma = 0, 1
    x_theory = np.linspace(X.min(), X.max(), num=int(5e3))
    y_theory = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(
        - (x_theory - mu)**2 / (2*sigma**2)
    )

    n_pixels = len(subimg_normalized.flatten())
    ax.plot(x_theory, y_theory*n_pixels, label='theory')

    ax.set_xlim([-4*sigma, 4*sigma])

    ax.legend(loc='best')
    ax.set_xlabel('normalized difference pixel value')
    ax.set_ylabel('number of pixels')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))

    savdir = os.path.dirname(rsubimgfile)
    savname = (
        os.path.basename(rsubimgfile).replace(
            '.fits','-normalized_histogram.png')
    )
    savpath = os.path.join(savdir, savname)

    f.tight_layout()
    f.savefig(savpath, dpi=200, bbox_inches='tight')
    print('%sZ: wrote %s' % (datetime.utcnow().isoformat(), savpath))


def is_image_noise_gaussian(
    fitsdir, projectid, field, camera, ccd,
    photrefdir='/nfs/phtess1/ar1/TESS/SIMFFI/BASE/reference-frames/'):
    '''
    The noise in the differenced image should be gaussian.  Oelkers & Stassun
    (2018) suggest the following approach to check whether it is.

    For each subtracted frame, normalize it by the combination of the noise
    from the science frame  and the photref frame:

        expected noise = sqrt( science frame  + photref frame ).

    The normalized pixels should be a gaussian centered at zero, with a std
    devn of 1.
    '''

    rsubglob = 'rsub-*-tess*-xtrns.fits'
    rsubfiles = np.sort(glob(fitsdir+rsubglob))

    imgglob = 'tess*-xtrns.fits'
    sciimgfiles = np.sort(glob(fitsdir+imgglob))

    # e.g., proj43-orbit-10-cam1-ccd2-combinedphotref-onenight.fits
    photrefglob = ('proj{:d}-{:s}-cam{:s}-ccd{:s}-*.fits'.format(
                   projectid,field,str(camera),str(ccd)))
    photreffile = np.sort(glob(photrefdir+photrefglob))

    if len(photreffile) != 1:
        raise AssertionError('expected a single photometric reference')

    photreffile = photreffile[0]
    photrefimg, _ = iu.read_fits(photreffile, ext=0)

    for rsubimgfile, sciimgfile in zip(rsubfiles, sciimgfiles):

        subimg, subhdr = iu.read_fits(rsubimgfile, ext=0)
        sciimg, scihdr = iu.read_fits(sciimgfile, ext=0)

        expected_noise = np.sqrt(photrefimg + sciimg)

        subimg_normalized = subimg / expected_noise

        _plot_normalized_subtractedimg_histogram(
            subimg_normalized, rsubimgfile)


def main(fitsdir, fitsglob, projectid, field, outdir=sv.REDPATH,
         lcdirectory=None, nworkers=1,
         aperturelist='1.45:7.0:6.0,1.95:7.0:6.0,2.45:7.0:6.0',
         kernelspec='b/4;i/4;d=4/4', convert_to_fitsh_compatible=True,
         anetfluxthreshold=20000, anettweak=6, initccdextent='0:2048,0:2048',
         anetradius=30,
         zeropoint=11.82,
         epdsmooth=21, epdsigclip=10, photdisjointradius=2,
         tuneparameters='true', is_ete6=True,
         catalog_faintrmag=13, fistarfluxthreshold=1000,
         photreffluxthreshold=1000, extractsources=True, binlightcurves=False
         ):
    '''
    args:

        projectid (int): ...

        field (str): ...

    kwargs:

        initccdextent (string): <x1>:<x2>,<y1>:<y2> for initial source
        extraction and astrometry.

        anetradius (float): how far apart the index files used for anet can be.
        make this a good deal bigger than the FoV of the image you're using,
        especially if it's wide-field (> 10deg x 10deg).

        anetfluxthreshold (float/int): this value sets the minimum flux for
        bright stars used in the astrometric solution. if it's too low, then
        you get too many stars, and astrometry.net gets confused.

        zeropoint (float): 11.82, from "get_zero_point.py" in ete6 repo.
        since we're doing relative photometry, and the same stars are always in
        the same cameras, the exact value isn't very important.

        extractsources (bool): if True, uses fistar for source extraction
        (with the fistarfluxthreshold). if False, uses the background catalog
        (e.g., 2MASS), projected onto the frame with WCS, and the
        catalog_faintrmag cutoff to make the list of sources.
    '''

    starttime = datetime.utcnow()

    ###########################################################################
    # get list of ete6 reduced images. (different format from images fitsh can
    # work with). trim images, and save to a single-extension fits file.

    if is_ete6:
        ete6_reddir = '/nfs/phtess1/ar1/TESS/SIMFFI/RED/'
        ete6_list = np.sort(
            glob(ete6_reddir+fitsglob.replace('_cal_img.fits','-s_ffic.fits'))
        )
    else:
        raise NotImplementedError('need path for real images')

    if tuneparameters=='true':
        # select 150 sequential images for pipeline tuning
        ete6_list = ete6_list[400:550]
    else:
        pass

    if convert_to_fitsh_compatible:
        tu.from_ete6_to_fitsh_compatible(ete6_list, outdir)

    ###########################################################################

    # TESS specific variables.
    ccd_fov = 12 # degrees. 24/2.

    fits_list = np.sort(glob(fitsdir + fitsglob))

    # get gain, ccdextent, zeropoint, exposure time from header.
    rand_fits = fits_list[np.random.randint(0, high=len(fits_list))]
    hdu_list = fits.open(rand_fits)
    hdr = hdu_list[0].header

    ccdgain = hdr['GAINA'] # electrons/count, from CCD output A. (there are four!)
    exptime = int(np.round(hdr['TELAPSE']*24*60*60)) # in seconds, 1800
    ra_nom = hdr['CRVAL1']  # RA at CRPIX1, CRPIX2. Not "camera boresight".
    dec_nom = hdr['CRVAL2'] # DEC at CRPIX1, CRPIX2

    catalog, catra, catdec, catbox = '2MASS', ra_nom, dec_nom, ccd_fov
    catalog_file = (
        '%s-RA%s-DEC%s-SIZE%s.catalog' % (catalog, catra, catdec, catbox)
    )
    reformed_cat_file = catalog_file.replace('.catalog', '.reformed_catalog')

    catalog_file = outdir + catalog_file
    reformed_cat_file = outdir + reformed_cat_file

    ###########################################################################

    if not is_presubtraction_complete(outdir, fitsglob, lcdirectory,
                                      extractsources=extractsources):
        get_files_needed_before_image_subtraction(
            fitsdir, fitsglob, outdir, initccdextent, ccdgain, zeropoint, exptime,
            ra_nom, dec_nom, catra, catdec, ccd_fov, catalog, catalog_file,
            reformed_cat_file, fnamestr=fitsglob,
            anetfluxthreshold=anetfluxthreshold, anetradius=anetradius,
            fistarglob='*.fistar', width=13, anettweak=anettweak, xpix=2048,
            ypix=2048, cols=(2,3), brightrmag=6.0, faintrmag=catalog_faintrmag,
            fistarfluxthreshold=fistarfluxthreshold, aperturelist=aperturelist,
            nworkers=nworkers, extractsources=extractsources)

    else:
        print('found fistar, fiphot, and wcs files. proceeding to image '
              'subtraction.')

    if not initial_wcs_worked_well_enough(outdir, fitsglob):
        raise AssertionError('if initial wcs failed on majority of frames, '
                             'you need to fix that.')

    #############################################################
    # run image subtraction convolution, then do the photometry #
    #############################################################

    camera = int(fitsglob.split('-')[1])
    ccd = int(fitsglob.split('-')[2])

    fieldinfo = {}
    fieldinfo['camera'] = camera
    fieldinfo['ccd'] = ccd
    fieldinfo['projectid'] = projectid
    fieldinfo['field'] = field

    photparams = {
        'ccdgain': ccdgain, 'ccdexptime': exptime, 'zeropoint': zeropoint
    }

    photreftype, dbtype = 'onenight', 'postgres'

    epdlcglob, tfalcglob = '*.epdlc', '*.tfalc'
    statsdir = os.path.dirname(lcdirectory) + '/stats_files/'
    for dirname in [lcdirectory, statsdir]:
        if not os.path.exists(dirname):
            os.mkdir(dirname)
    epdstatfile = statsdir + 'camera' + str(camera) + '_ccd' + str(ccd) + '.epdstats'
    tfastatfile = statsdir + 'camera' + str(camera) + '_ccd' + str(ccd) + '.tfastats'

    xtrnsglob = fitsglob.replace('.fits','-xtrns.fits')
    iphotpattern = fitsdir+'rsub-????????-'+fitsglob.replace('.fits','.iphot')

    if not is_imagesubtraction_complete(fitsdir, fitsglob, lcdirectory):

        run_imagesubtraction(fitsdir, fitsglob, fieldinfo, photparams,
                             fits_list, photreftype, dbtype, reformed_cat_file,
                             xtrnsglob, iphotpattern, lcdirectory,
                             kernelspec=kernelspec, refdir=sv.REFBASEDIR,
                             nworkers=nworkers, aperturelist=aperturelist,
                             photdisjointradius=photdisjointradius,
                             colorscheme='bwr',
                             photreffluxthreshold=photreffluxthreshold,
                             extractsources=extractsources)
    else:
        print('found that image subtraction is complete.')

    run_detrending(epdstatfile, tfastatfile, lcdirectory, epdlcglob,
                   reformed_cat_file, statsdir, field, epdsmooth=epdsmooth,
                   epdsigclip=epdsigclip, nworkers=nworkers,
                   binlightcurves=binlightcurves)

    statsdir = os.path.dirname(epdstatfile)+'/'
    assess_run(statsdir, lcdirectory, starttime, field, fitsdir, projectid,
               field, camera, ccd, binned=False, make_whisker_plot=True)

    # TODO: maybe change the statsfile format, and include CDPP? or some
    # duration-aware RMS measure


def check_args(args):
    if not args.fitsdir:
        print('fitsdir must not be null')
    if not args.fitsglob:
        print('fitsglob must not be null')
    if not args.field:
        print('field must not be null')
    if not (args.tuneparameters=='true') or (args.tuneparameters=='false'):
        raise AssertionError('boolean tuneparameters must be set as string.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=('Given images, make lightcurves.')
    )

    parser.add_argument('--fitsdir', type=str, default=None,
        help=('e.g., /foo/FFI/RED'))

    parser.add_argument('--fitsglob', type=str, default=None,
        help=('e.g., *.fits.gz. Chosen frames are fitsdir+fitsglob.'))

    parser.add_argument('--projectid', type=int, default=1,
        help=('unique integer identifying this reduction run'))

    parser.add_argument('--field', type=str, default=None,
        help=('field identifier string, e.g., "sector1"'))

    parser.add_argument('--outdir', type=str, default=None,
        help=('e.g., /foo/FFI/RED'))

    parser.add_argument('--lcdirectory', type=str, default=None,
        help=('e.g., /foo/RED/LC'))

    parser.add_argument(
        '--aperturelist', type=str,
        default="1.45:7.0:6.0,1.95:7.0:6.0,2.45:7.0:6.0",
        help=('formatted as "inner radius: outer radius: bkgd annulus width"')
    )

    parser.add_argument(
        '--kernelspec', type=str, default='b/4;i/4;d=4/4',
        help=('https://fitsh.net/wiki/man/ficonv explains the formatting.'
              'i/<spatial order> identity kernel (a.k.a. flux term) with the'
              'specified order of polynomial spatial variations.'
              'b/<spatial order> constant offset kernel (a.k.a. background'
              'term) with the specified order of polynomial spatial'
              'variations'
              'd=<size>/<spatial order> discrete kernel with the half-size of'
              '<size> and the specified order of polynomial spatial variations'
              'g=<size>,<sigma>,<order>/<spatial order> Gaussian kernel with'
              'the half-size of <size>, standard deviation of <sigma> and'
              'Hermite basis order of <order>, with the specified order of'
              'polynomial spatial variations'
             )
    )

    parser.add_argument(
        '--photdisjointradius', type=int, default=2,
        help=('https://fitsh.net/wiki/man/fiphot gives details. '
              'During the bacground determination on the aperture annuli, '
              'omit the pixels which are closer to the other centroids than '
              'the specified radius.'
             )
    )

    parser.add_argument('--anetfluxthreshold', type=int, default=20000,
        help=(
            'flux threshold used to identify bright stars in anet '
            'initial wcs solution attempt of frames.'
        )
    )
    parser.add_argument('--anettweak', type=int, default=6,
        help=(
            'SIP (Simple Imaging Polynomial) order. Higher = more flexibility '
            'for frame distortion. Max of 10.'
        )
    )
    parser.add_argument('--anetradius', type=float, default=30.,
        help=(
            'Radius (in deg) over which to search for index files in '
            'astrometry.net or net. For TESS, you want to go big -- e.g., 30.'
        )
    )
    parser.add_argument(
        '--initccdextent', type=str,
        default="0:2048,0:2048",
        help=(' section <x1>:<x2>,<y1>:<y2> for initial source extraction and '
              'astrometry. (not for initial photometry).')
    )

    parser.add_argument('--epdsmooth', type=int, default=21,
        help=(
            'number of cadences used when passing a median filter over the'
            'magnitude time series before EPD (see '
            'imagesubphot.epd_magseries_imagesub). For HATPI default is 21, or'
            ' a ~10 minute median filter. For TESS FFIs, try 11 = 5.5hrs.'
        )
    )
    parser.add_argument('--epdsigclip', type=float, default=10.0,
        help=(
            'sigma clipping aplpied to epd lightcurves. assumed symmetric.'
        )
    )

    parser.add_argument('--nworkers', type=int, default=1,
        help=('number of workers to thread over'))

    parser.add_argument(
        '--convert_to_fitsh_compatible', dest='cfc', action='store_true',
        help=('converts calibrated frames to single-fits-header fitsh readable'
              'objects. if already done, will find them and skip.')
    )
    parser.add_argument(
        '--no-convert_to_fitsh_compatible', dest='cfc', action='store_false',
        help=('don\'t do fitsh-compatibility conversion.')
    )
    parser.set_defaults(cfc=True)

    parser.add_argument(
        '--binlightcurves', dest='binlcs', action='store_true',
        help=('will bin lightcurves to 1hr and 6hr cadence. only purpose: '
              'making plots to understand the red noise')
    )
    parser.add_argument(
        '--no-binlightcurves', dest='binlcs', action='store_false',
        help=('don\'t bin lightcurves (it takes a while).')
    )
    parser.set_defaults(binlcs=True)

    parser.add_argument(
        '--tuneparameters', type=str,
        default="true",
        help=('TUNING: iterate through different img subtraction parameters '
              'if true. Gain speed by selecting a small set of frames. '
              'Otherwise, run in FULL reduction mode.')
    )

    parser.add_argument('--catalog_faintrmag', type=float, default=13.0,
        help=('faint 2MASS catalog r magnitude used to make fovcat, '
              'which sets the stars that get photometered (TODO: is this '
              'true?)'))
    parser.add_argument('--fistarfluxthreshold', type=int, default=1000,
        help=('faint flux threshold used to do initial raw photometry on '
              'frame, in pre-img subtraction.'
             ))
    parser.add_argument('--photreffluxthreshold', type=int, default=1000,
        help=('faint flux threshold used to extract sources on photometric '
              'reference frame, in image subtraction. (these sources will '
              'then be subtracted against).'
             ))
    parser.add_argument('--extractsources', type=int, default=1,
        help=(
            'extractsources (int): if 1, uses fistar for source extraction '
            '(with the fistarfluxthreshold). if 0, use the background catalog '
            '(e.g., 2MASS), projected onto the frame with WCS, and the '
            'catalog_faintrmag cutoff to make the list of sources. ')
         )

    args = parser.parse_args()

    check_args(args)

    if args.extractsources not in [0,1]:
        raise AssertionError('extractsources must be 0 or 1')
    if args.extractsources == 1:
        extractsources = True
    else:
        extractsources = False

    main(args.fitsdir, args.fitsglob, args.projectid, args.field,
         outdir=args.outdir, lcdirectory=args.lcdirectory,
         nworkers=args.nworkers, aperturelist=args.aperturelist,
         convert_to_fitsh_compatible=args.cfc, kernelspec=args.kernelspec,
         epdsmooth=args.epdsmooth, epdsigclip=args.epdsigclip,
         photdisjointradius=args.photdisjointradius,
         tuneparameters=args.tuneparameters,
         anetfluxthreshold=args.anetfluxthreshold,
         anettweak=args.anettweak, initccdextent=args.initccdextent,
         anetradius=args.anetradius,
         catalog_faintrmag=args.catalog_faintrmag,
         fistarfluxthreshold=args.fistarfluxthreshold,
         photreffluxthreshold=args.photreffluxthreshold,
         extractsources=extractsources,
         binlightcurves=args.binlcs
    )

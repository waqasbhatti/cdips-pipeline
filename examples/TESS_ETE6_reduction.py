'''
Start with: calibrated ETE6 frames optained from MAST.

Then do:
    - initial calibration steps
    - image subtraction photometry
    - detrend the resulting lightcurves.

TODO:
* currently reducing only channel 1 from camera 1.
* what is the difference between CCDA, CCDB, CCDC, CCDD?

SOMEDAY:
* manage the whole thing with airflow

====================

USAGE:

from within `examples/` subdir, in the "trex_27" environment,
$ python TESS_ETE6_reduction.py

If you want to do this iteratively, you'll probably need to hard-delete:
    * all rows of `astromrefs`, `photrefs`, and `calibratedframes` SQL tables,
    * "/LC/stats_file*.epdstats"
    * "/LC/stats_file*.tfastats"
    * all epdlcs and tfalcs, if you want to remake them (e.g., you have new
      data points)

because smart overwriting has not yet been implemented.

If you want to convert tess to fitsh frames, set
    convert_to_fitsh_compatible = True

====================
'''

import os
import numpy as np
import aperturephot as ap, shared_variables as sv, autoimagesub as ais, \
       imagesubphot as ism, tessutils as tu
from glob import glob
from tqdm import tqdm
from astropy.io import fits

# insane "globals" to prevent rewriting all the code
camera = 1
ccd = 1
projectid = 42 # an integer, to identify the "project" (mostly avoids rewriting code)
field = 'ete6_field0'

if __name__ == '__main__':

    np.random.seed(42)

    ###########################################################################
    ###########################################################################
    # Job specific variables
    nworkers = 20

    # TESS specific variables.
    convert_to_fitsh_compatible = False
    tess_fov = 24 # degrees
    zeropoint = 11.82 # via "get_zero_point.py" in ete6 repo

    fitsdir = sv.REDPATH
    fitsglob = sv.LOCAL_GLOBPATTERN
    fits_list = np.sort(glob(fitsdir + fitsglob))

    # get gain, ccdextent, zeropoint, exposure time.
    rand_fits = fits_list[np.random.randint(0, high=len(fits_list))]
    hdu_list = fits.open(rand_fits)
    hdr = hdu_list[0].header

    ccdgain = hdr['GAINA'] # electrons/count, from CCD output A. (there are four!)
    exptime = int(np.round(hdr['TELAPSE']*24*60*60)) # in seconds, 1800
    ra_nom = hdr['CRVAL1']  # RA at CRPIX1, CRPIX2. Not "camera boresight".
    dec_nom = hdr['CRVAL2'] # DEC at CRPIX1, CRPIX2

    ccdextent = '0:2048,0:2048'

    # get list of ete6 reduced images.
    ete6_reddir = '/nfs/phtess1/ar1/TESS/SIMFFI/RED/'
    outdir = sv.REDPATH
    ete6_list = np.sort(glob(
                    ete6_reddir+'tess?????????????-1-1-0016-s_ffic.fits'))
    ete6_list = ete6_list[100:400] #FIXME temporary, to run on 300 frames w ISP

    catalog = '2MASS'
    catra, catdec = ra_nom, dec_nom
    catbox=tess_fov
    catalog_file = '%s-RA%s-DEC%s-SIZE%s.catalog' % (catalog, catra, catdec,
                                                     catbox)
    reformed_cat_file = catalog_file.replace('.catalog', '.reformed_catalog')

    catalog_file = sv.REDPATH + catalog_file
    reformed_cat_file = sv.REDPATH + reformed_cat_file

    # trim images, and save to a single-extension fits file.
    if convert_to_fitsh_compatible:
        tu.from_ete6_to_fitsh_compatible(ete6_list, outdir)

    ### ########################################################################
    ### # get .fistar, .fiphot, and .wcs files needed before image subtraction #
    ### ########################################################################
    ### # 1. run parallel_extract_sources on all frames with threshold ~ 10000 to get
    ### #    bright stars for astrometry.
    ### # 2. run parallel_anet to get precise WCS headers for all frames.
    ### # 3. run make_fov_catalog to get a FOV source catalog for the field.
    ### # 4. run reform_fov_catalog to cut this down to the columns needed for magfit
    ### #    only.
    ### # 5. run parallel_fitsdir_photometry for photometry on all frames

    ### ap.parallel_extract_sources(fitsdir, outdir, ccdextent=ccdextent,
    ###                             ccdgain=ccdgain, fluxthreshold=10000,
    ###                             zeropoint=zeropoint, exptime=exptime,
    ###                             tailstr='.fits',
    ###                             fnamestr='*-1-1-0016_cal_img.fits')

    ### # Really width=24 deg for TESS frames. But this makes anet fail. Instead,
    ### # get wcs for the inside of the frame. [and stuff on the edge less
    ### # accurate?]
    ### ap.parallel_anet(fitsdir, outdir, ra_nom, dec_nom, fistarglob='*.fistar',
    ###                  infofromframe=False, width=13, tweak=6, radius=13,
    ###                  xpix=2048, ypix=2048,
    ###                  cols=(2,3) # columns with x,y in fistar file.
    ###                 )

    ### # This function gets all the sources in the field of view of the frame, given
    ### # its central pointing coordinates and plate-scale from 2MASS. This catalog
    ### # file is then be used as input to make_source_list below.
    ### _ = ap.make_fov_catalog(ra=catra, dec=catdec, size=tess_fov,
    ###                         brightrmag=6.0, faintrmag=13.0, fits=None,
    ###                         outfile=None, outdir=outdir, catalog=catalog,
    ###                         catalogpath=None, columns=None,
    ###                         observatory='tess')

    ### # This converts the full output catalog from 2massread, etc. to the format
    ### # required for magfit. Also useful for general reforming of the columns.
    ### ap.reform_fov_catalog(catalog_file, reformed_cat_file)

    ### ap.parallel_fitsdir_photometry(fitsdir, outdir, reformed_cat_file,
    ###                                fluxthreshold=1000.0,
    ###                                ccdextent={'x':[0.,2048.],'y':[0.,2048.]},
    ###                                pixborders=0.0,
    ###                                aperturelist='1.45:7.0:6.0,1.95:7.0:6.0,2.45:7.0:6.0',
    ###                                removesourcetemp=True,
    ###                                removesourcelist=False, binaryoutput=False,
    ###                                nworkers=nworkers, maxtasksperworker=1000,
    ###                                saveresults=True, rejectbadframes=True,
    ###                                minsrcbgv=200.0, maxmadbgv=150.0,
    ###                                maxframebgv=2000.0, minnstars=500,
    ###                                fitsglob=fitsglob, ccdgain=ccdgain,
    ###                                ccdexptime=exptime, zeropoint=zeropoint)


    ##########################################################
    # run image subtraction convolution, then the photometry #
    ##########################################################
    fovcatalog = catalog_file

    lcdirectory = sv.LCPATH + "ISP_1-1/"

    fieldinfo = {}
    fieldinfo['camera'] = camera
    fieldinfo['ccd'] = ccd
    fieldinfo['projectid'] = projectid
    fieldinfo['field'] = field

    photparams = {'ccdgain': ccdgain, 'ccdexptime': exptime,
                  'zeropoint': zeropoint}

    photreftype = 'onenight'
    dbtype = 'postgres'

    epdlcglob = '*.epdlc'
    tfalcglob = '*.tfalc'
    statspath = sv.LCPATH + 'stats_files/'
    epdstatfile = statspath + 'camera' + str(camera) + '_ccd' + str(ccd) + '.epdstats'
    tfastatfile = statspath + 'camera' + str(camera) + '_ccd' + str(ccd) + '.tfastats'

    xtrnsglob = sv.LOCAL_GLOBPATTERN.replace('.fits','-xtrns.fits')
    iphotpattern = sv.REDPATH+'rsub-????????-'+sv.LOCAL_GLOBPATTERN.replace('.fits','.iphot')

    ### # FIXME: commented if done, but overwrites have not been fixed to work
    ### # Step ISP0.
    ### _ = ais.parallel_frames_to_database(fitsdir, 'calibratedframes',
    ###                                     observatory='tess', fitsglob=fitsglob,
    ###                                     overwrite=False,
    ###                                     badframetag='badframes',
    ###                                     nonwcsframes_are_ok=False,
    ###                                     nworkers=nworkers, maxworkertasks=1000)

    ### # Step ISP1.
    ### _ = ais.dbgen_get_astromref(fieldinfo, makeactive=True, observatory='tess',
    ###                             overwrite=False, refdir=sv.REFBASEDIR,
    ###                             database=None)

    ### # Step ISP2. Takes ~1 sec per 20-30 frames.
    ### _ = ism.get_smoothed_xysdk_coeffs(fitsdir, fistarglob='*.fistar',
    ###                                   nworkers=nworkers, maxworkertasks=1000)

    ### # Step ISP3. ~600 in 5 minutes --> 2 frames per second, running over 20 workers.
    ### # If called with warpcheck=True and warpthreshold=2000, many things will be
    ### # moved to badframes that are not in fact badframes. The warpthreshold is a bad
    ### # empirical thing that should be avoided.
    ### _ = ais.framelist_make_xtrnsfits(fits_list, outdir=None, refinfo='foobar',
    ###                                  warpcheck=False, warpthreshold=15000.0,
    ###                                  warpmargins=100, nworkers=nworkers,
    ###                                  observatory='tess', maxworkertasks=1000,
    ###                                  fieldinfo=fieldinfo)

    ### # # Optional Step ISP3.5: parallelized move of astrometry ref shifted frames to database
    ### # out = ais.parallel_frames_to_database(fitsdir, 'arefshifted_frames',
    ### #                                       fitsglob='1-???????_?-xtrns.fits',
    ### #                                       network='HP', overwrite=False,
    ### #                                       badframetag='badframes',
    ### #                                       nonwcsframes_are_ok=False, nworkers=nworkers,
    ### #                                       maxworkertasks=1000)
    ### # 

    ### # Step ISP4: the next thing to do is to select a bunch of frames that can serve
    ### # as photometric reference frames (photrefs).

    ### xtrnsfiles = glob(fitsdir+xtrnsglob)
    ### photrefinfo = ais.generate_photref_candidates_from_xtrns(xtrnsfiles,
    ###                                                          minframes=50,
    ###                                                          observatory='tess',
    ###                                                          maxbackgroundstdevpctile=100.,
    ###                                                          maxbackgroundmedianpctile=70.,
    ###                                                          minngoodobjectpctile=70.,
    ###                                                          forcecollectinfo=False,
    ###                                                          nworkers=nworkers,
    ###                                                          maxworkertasks=1000)

    ### # Optional Step ISP5: amend the list, if needed.
    ### # photrefinfo = ais.amend_candidate_photrefs(photrefinfo)

    ### # Step ISP6: make a photometric reference frame
    ### _ = ais.generate_combined_photref(photrefinfo, photreftype, dbtype,
    ###                                   photref_reformedfovcat=reformed_cat_file,
    ###                                   makeactive=True, field=None, ccd=None,
    ###                                   projectid=None, combinemethod='median',
    ###                                   kernelspec='b/4;i/4;d=4/4',
    ###                                   ccdgain=ccdgain, zeropoint=zeropoint,
    ###                                   ccdexptime=exptime, extractsources=True,
    ###                                   astrometrysrcthreshold=25000,
    ###                                   apertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
    ###                                   framewidth=None, searchradius=8.0,
    ###                                   nworkers=nworkers, maxworkertasks=1000,
    ###                                   observatory='tess', fieldinfo=fieldinfo)

    ### # Step ISP7: convolve and subtract all FITS files in the xtrnsfits list from the
    ### # photometric reference.  With 30 workers, at best process ~few frames per
    ### # second.

    ### _ = ais.parallel_xtrnsfits_convsub(xtrnsfiles, photreftype, outdir=None,
    ###                                    observatory='tess', fieldinfo=fieldinfo,
    ###                                    reversesubtract=True,
    ###                                    kernelspec='b/4;i/4;d=4/4',
    ###                                    nworkers=nworkers, maxworkertasks=1000)

    ### # Step ISP8: do photometry on your subtracted frames to produce .iphot files.
    ### # With 30 workers, at best process ~few frames per second.

    ### subfitslist = glob(sv.REDPATH+'rsub-????????-'+
    ###                    sv.LOCAL_GLOBPATTERN.replace('.fits','-xtrns.fits'))
    ### _ = ais.parallel_convsubfits_staticphot(subfitslist,
    ###                                         photreftype=photreftype,
    ###                                         kernelspec='b/4;i/4;d=4/4',
    ###                                         lcapertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
    ###                                         photdisjointradius=2, outdir=None,
    ###                                         fieldinfo=fieldinfo,
    ###                                         observatory='tess', nworkers=30,
    ###                                         maxworkertasks=1000,
    ###                                         photparams=photparams)

    ### # Step ISP9 + 10 : dump lightcurves.
    ### ism.dump_lightcurves_with_grcollect(iphotpattern, lcdirectory, '4g',
    ###                                     lcextension='grcollectilc',
    ###                                     objectidcol=3,
    ###                                     observatory='tess')

    ### # # Alternative Step 9 + 10: add the .iphot files to postgres database. Collect
    ### # # LCs from postgres, and then make difference image lightcurve (*.ilc) files in
    ### # # lcdirectory. AKA "light curve dump", or "the transposition problem".
    ### # # Surprisingly computationally expensive.  An alternative approach is
    ### # # fitsh's `grcollect`, which skips the database architecture entirely.
    ### # 
    ### # print('beginning insert_phots_into_database')
    ### # ais.insert_phots_into_database(sv.REDPATH, frameglob='rsub-*-xtrns.fits',
    ### #                                photdir=None, photglob='rsub-*-%s.iphot',
    ### #                                maxframes=None, overwrite=False, database=None)
    ### # 
    ### # hatidlist = ais.get_hatidlist_from_cmrawphot(projectid, field, ccd, photreftype)
    ### # 
    ### # print('beginning lightcurve dump')
    ### # ais.parallel_dbphot_lightcurves_hatidlist(hatidlist, lcdirectory)


    # Step ISP11: do EPD on all the LCs, and collect stats on the results.
    # for ISP LCs, use lcmagcols=([27,28,29],[30,],[30,],[30,])
    if not os.path.exists(epdstatfile):
        _ = ism.parallel_run_epd_imagesub(lcdirectory,
                                          ilcglob='*.grcollectilc',
                                          outdir=None, smooth=21,
                                          sigmaclip=5.0, nworkers=nworkers,
                                          maxworkertasks=1000, minndet=100)
        ap.parallel_lc_statistics(lcdirectory, epdlcglob, reformed_cat_file, tfalcrequired=False,
                                  fovcatcols=(0,9), # objectid, magcol to use
                                  fovcatmaglabel='r', outfile=epdstatfile, nworkers=nworkers,
                                  workerntasks=500, rmcols=[14,19,24],
                                  epcols=[27,28,29], tfcols=[30,31,32], rfcols=None,
                                  correctioncoeffs=None, sigclip=5.0)
        ap.plot_stats_file(epdstatfile, statspath, field, binned=False,
                           logy=True, logx=False, correctmagsafter=None,
                           rangex=(5.9,13.1), observatory='tess')
    else:
        print('already made EPD LC stats and plots')


    # Step ISP12: do TFA on all the LCs. First, choose TFA template stars using the
    # .epdlc stats. Then run TFA, to get .tfalc.TF{1,2,3} files. Turn them into
    # single .tfalc files. Then collect statistics.
    if not os.path.exists(lcdirectory+'aperture-1-tfa-template.list'):
        _ = ap.choose_tfa_template(epdstatfile, reformed_cat_file, lcdirectory,
                                   ignoretfamin=False, fovcat_idcol=0,
                                   fovcat_xicol=3, fovcat_etacol=4,
                                   fovcat_magcol=9, min_ndet=100,
                                   min_nstars=50, max_nstars=1000,
                                   brightest_mag=8.5, faintest_mag=12.0,
                                   max_rms=0.1, max_sigma_above_rmscurve=5.0,
                                   outprefix=None, tfastage1=True)
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
                                  tfcols=[30,30,30], #NOTE weird call b/c .TF[1-3]
                                  rfcols=None,
                                  correctioncoeffs=None,
                                  sigclip=5.0)
        ap.plot_stats_file(tfastatfile, statspath, field, binned=False,
                           logy=True, logx=False, correctmagsafter=None,
                           rangex=(5.9,13.1), observatory='tess')
    else:
        print('already made TFA LC stats and plots')

    #####################################################################
    # run detrending on raw photometry, to compare vs. image subtracted #
    #####################################################################
    # Steps [1-5] of raw photometry (source extraction, astrometry, creation of
    # a source catalog, reforming the source catalog, and raw photometry) have
    # been performed.  To detrend, we first "magnitude fit" (see Sec 5.5 of
    # Zhang, Bakos et al 2016).  This procedure empirically fits magnitude
    # differences (vs. a reference) for each star in each frame, using
    # on-camera position, catalog magnitude, catalog color, and subpixel
    # position. Then we run EPD and TFA. Procedurally:
    # 6. run get_magfit_frames to select a single magfit photometry reference
    #    and set up per-CCD work directories, symlinks, etc. for the next 
    #    steps.
    # 7. run make_magfit_config to generate magfit config files for
    #    MagnitudeFitting.py
    # 8. run make_fiphot_list to make lists of fiphot files for each CCD.
    # 9. run MagnitudeFitting.py in single reference mode.
    # 10. run do_masterphotref.py to get the master mag fit reference.
    # 11. run MagnitudeFitting.py in master reference mode.
    # 12. run parallel_collect_lightcurves to collect all lightcurves into .rlc
    #     files.
    # 13. run serial_run_epd or parallel_run_epd to do EPD on all LCs.
    # 14. run parallel_lc_statistics to collect stats on .epdlc files.
    # 15. run choose_tfa_template to choose TFA template stars using the .epdlc
    #     stats.
    # 16. run parallel_run_tfa for TFA to get .tfalc.TF{1,2,3} files (FIXME: still
    #     need to collect into single .tfalc files for all apertures)
    # 17. run parallel_lc_statistics to collect stats on .tfalc files.
    # 20. run plot_stats_file to make RMS vs. mag plots for all unbinned and binned
    #     LCs.

    ### fitsdir = '/home/lbouma/proj/ete6/data/RED_1-1_SAP'
    ### lcdirectory = sv.LCPATH + "SAP_1-1/"

    ### mfdict = ap.get_magfit_frames(fitsdir, sv.LOCAL_GLOBPATTERN, fitsdir,
    ###                               selectreference=True, linkfiles=False,
    ###                               framestats=False, observatory='tess')

    ### #FIXME from here.
    ### # Joel notes: magnitude fitting is probably not needed for space data
    ### ap.make_magfit_config()

    ### # ap.make_fiphot_list()

    ### # ap.run_magfit()

    ### # ap.get_master_photref()

    ### # ap.run_magfit()

    ### ap.parallel_collect_aperturephot_lightcurves(fitsdir, lcdirectory,
    ###                                              photext='fiphot',
    ###                                              skipcollectedlcs=False,
    ###                                              nworkers=nworkers)

    ### # ap.parallel_run_epd()

    ### # ap.parallel_lc_statistics()

    ### # ap.choose_tfa_template()

    ### # ap.parallel_run_tfa()

    ### # ap.parallel_lc_statistics()

    ### # ap.plot_stats_file()

'''
Starting from calibrated frames, this script does image subtraction, and
detrends the resulting lightcurves.

Commented-out "steps" must be run; they are simply commented out because
this script is not yet runnable as a "push this button, and it will do
everything for you" utility.

The current usage involves running

$ python example_field_577_reduction.py

or something analogous, with whatever steps you _want_ to do commented in, and
the others commented out.

TODO: make it so that one could hypothetically push this button, and have a
script like this spit out good lightcurves.

'''

from glob import glob
from tqdm import tqdm
import aperturephot as ap, shared_variables as sv, autoimagesub as ais, \
       imagesubphot as ism
import os

fitsbasedir = sv.REDPATH

fovcatalog = sv.CATPATH + 'G1830-2230_577-gri-20.0.catalog' #NOTE could do better

projectid = 12
field = 'G1830-2230_577'
ccd = 8

fitsglob = '1-???????_?.fits'
fitsfiles = glob(fitsbasedir+fitsglob)

photreftype = 'onenight'
dbtype = 'postgres'

dirstr = 'lgb_projid{:s}_ccd{:s}_field{:s}/'.format(str(projectid), str(ccd), field)
lcdirectory = sv.LCPATH + dirstr

epdlcglob = '*.epdlc'
tfalcglob = '*.tfalc'
statspath = sv.LCPATH + 'stats_files/'
epdstatfile = statspath + field + '_' + str(projectid) + '_' + str(ccd) + '.epdstats'
tfastatfile = statspath + field + '_' + str(projectid) + '_' + str(ccd) + '.tfastats'


# # Step -1. (Before beginning "image subtraction" steps)
# ap.make_fov_catalog(faintrmag=15,
#                     fits=sv.REDPATH+'1-430036c_8.fits', #also selected as a photref
#                     outdir=sv.CATPATH)

# # Step -0.5
# incat = sv.CATPATH+'2MASS-RA277.5-DEC-22.5-SIZE14.0.catalog'
# outcat = sv.CATPATH+'G1830-2230_577-gri-20.0.catalog' # these names matter! used later
# ap.reform_fov_catalog(incat, outcat)

# Step 0.
# out = ais.parallel_frames_to_database(fitsbasedir, 'calframes',
#                                       fitsglob='1-???????_?.fits',
#                                       network='HP', overwrite=False,
#                                       badframetag='badframes',
#                                       nonwcsframes_are_ok=False, nworkers=40,
#                                       maxworkertasks=1000)

# # Step 1.
# arefinfo = ais.dbgen_astromref_projectidfieldccd(projectid, field, ccd,
#                                                  makeactive=True,
#                                                  overwrite=False,
#                                                  refdir=sv.REFBASEDIR,
#                                                  database=None)

# # Step 2. Takes ~1 sec per 20-30 frames.
# smoothedcoeffs = ism.get_smoothed_xysdk_coeffs(fitsbasedir,
#                                                fistarglob='*.fistar',
#                                                nworkers=40,
#                                                maxworkertasks=1000)


# Step 3. ~600 in 5 minutes --> 2 frames per second, running over 20 workers.
# If called with warpcheck=True and warpthreshold=2000, many things will be
# moved to badframes that are not in fact badframes. The warpthreshold is a bad
# empirical thing that should be avoided.
framelist_xtrnsfits = ais.framelist_make_xtrnsfits(fitsfiles, outdir=None,
                                                   refinfo='foobar',
                                                   warpcheck=False,
                                                   warpthreshold=15000.0,
                                                   warpmargins=100,
                                                   nworkers=20,
                                                   maxworkertasks=1000)

# # Step 3.5: parallelized move of astrometry ref shifted frames to database
# out = ais.parallel_frames_to_database(fitsbasedir, 'arefshifted_frames',
#                                       fitsglob='1-???????_?-xtrns.fits',
#                                       network='HP', overwrite=False,
#                                       badframetag='badframes',
#                                       nonwcsframes_are_ok=False, nworkers=40,
#                                       maxworkertasks=1000)
# 

# Step 4: the next thing to do is to select a bunch of frames that can serve
# as photometric reference frames (photrefs).

xtrnsglob = '1-???????_?-xtrns.fits'
xtrnsfiles = glob(fitsbasedir+xtrnsglob)

# photrefinfo = ais.generate_photref_candidates_from_xtrns(xtrnsfiles,
#                                                          minframes=50,
#                                                          maxhourangle=3.0,
#                                                          maxmoonphase=25.0,
#                                                          maxmoonelev=0.0,
#                                                          maxzenithdist=30.0,
#                                                          maxbackgroundstdevpctile=100.,
#                                                          maxbackgroundmedianpctile=20.,
#                                                          minngoodobjectpctile=85.,
#                                                          forcecollectinfo=False,
#                                                          nworkers=40,
#                                                          maxworkertasks=1000)

# Step 5: amend the list, if needed. NOTE: best to avoid this step if possible.
# It means more manual work.

#photrefinfo = ais.amend_candidate_photrefs(photrefinfo)

# Step 6: make a photometric reference frame

# photrefinfo = ais.generate_combined_photref(photrefinfo, photreftype, dbtype,
#                                             makeactive=True, field=None,
#                                             ccd=None, projectid=None,
#                                             combinemethod='median',
#                                             kernelspec='b/4;i/4;d=4/4',
#                                             ccdgain=None, zeropoint=None,
#                                             ccdexptime=None,
#                                             extractsources=True,
#                                             astrometrysrcthreshold=25000,
#                                             apertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
#                                             framewidth=None, searchradius=8.0,
#                                             nworkers=20, maxworkertasks=1000)
# 

# Step 7: convolve and subtract all FITS files in the xtrnsfits list from the
# photometric reference.  With 30 workers, at best process ~few frames per
# second.

results = ais.parallel_xtrnsfits_convsub(xtrnsfiles, photreftype, outdir=None,
                                         reversesubtract=True,
                                         kernelspec='b/4;i/4;d=4/4',
                                         nworkers=20, maxworkertasks=1000)

# Step 8: do photometry on your subtracted frames to produce .iphot files.
# With 30 workers, at best process ~few frames per second.

subfitslist = glob(sv.REDPATH+'rsub-????????-?-???????_?-xtrns.fits')

results = ais.parallel_convsubfits_staticphot(subfitslist,
                                              photreftype=photreftype,
                                              kernelspec='b/4;i/4;d=4/4',
                                              lcapertures='1.95:7.0:6.0,2.45:7.0:6.0,2.95:7.0:6.0',
                                              photdisjointradius=2,
                                              outdir=None, nworkers=30,
                                              maxworkertasks=1000)

# # Step 9 + 10 : dump lightcurves.
# 
# photfiles = glob(sv.REDPATH+'rsub-????????-?-???????_?.iphot')[:10]
# ism.dump_lightcurves_with_grcollect(photfiles, lcdirectory,
#                                     '4g',
#                                     lcextension='grcollectilc')

#    # Alternative Step 9 + 10: add the .iphot files to postgres database. Collect
#    # LCs from postgres, and then make difference image lightcurve (*.ilc) files in
#    # lcdirectory. AKA "light curve dump", or "the transposition problem".
#    # Surprisingly computationally expensive.  An alternative approach is
#    # fitsh's `grcollect`, which skips the database architecture entirely.
#    
#    print('beginning insert_phots_into_database')
#    ais.insert_phots_into_database(sv.REDPATH, frameglob='rsub-*-xtrns.fits',
#                                   photdir=None, photglob='rsub-*-%s.iphot',
#                                   maxframes=None, overwrite=False, database=None)
#    
#    hatidlist = ais.get_hatidlist_from_cmrawphot(projectid, field, ccd, photreftype)
#    
#    print('beginning lightcurve dump')
#    ais.parallel_dbphot_lightcurves_hatidlist(hatidlist, lcdirectory)


# Step 11: do EPD on all the LCs, and collect stats on the results.
results = ism.parallel_run_epd_imagesub(lcdirectory, ilcglob='*.ilc', outdir=None,
                                        smooth=21, sigmaclip=5.0, nworkers=20,
                                        maxworkertasks=1000)

# for ISM LCs, use lcmagcols=([27,28,29],[30,],[30,],[30,])
if not os.path.exists(epdstatfile):
    ap.parallel_lc_statistics(lcdirectory, epdlcglob, fovcatalog, tfalcrequired=False,
                              fovcatcols=(0,9), # objectid, magcol to use
                              fovcatmaglabel='r', outfile=epdstatfile, nworkers=20,
                              workerntasks=500, rmcols=[14,19,24],
                              epcols=[27,28,29], tfcols=[30,31,32], rfcols=None,
                              correctioncoeffs=None, sigclip=5.0)
else:
    print('already made EPD LC stats')

ap.plot_stats_file(epdstatfile, statspath, field, binned=False, logy=True,
                   logx=False, correctmagsafter=None, rangex=(5.9,13.1))

# Step 12: do TFA on all the LCs. First, choose TFA template stars using the
# .epdlc stats. Then run TFA, to get .tfalc.TF{1,2,3} files. Turn them into
# single .tfalc files. Then collect statistics.

if not os.path.exists(lcdirectory+'aperture-1-tfa-template.list'):
    tfadict = ap.choose_tfa_template(epdstatfile, fovcatalog, lcdirectory,
                                     ignoretfamin=False, fovcat_idcol=0,
                                     fovcat_xicol=3, fovcat_etacol=4,
                                     fovcat_magcol=9, min_ndet=100, min_nstars=50,
                                     max_nstars=1000, brightest_mag=8.5,
                                     faintest_mag=12.0, max_rms=0.1,
                                     max_sigma_above_rmscurve=5.0, outprefix=None,
                                     tfastage1=True)

templatefiles = glob(lcdirectory+'aperture-?-tfa-template.list')
ism.parallel_run_tfa(lcdirectory, templatefiles, epdlc_glob='*.epdlc',
                    epdlc_jdcol=0, epdlc_magcol=(27,28,29),
                    template_sigclip=5.0, epdlc_sigclip=5.0, nworkers=20,
                    workerntasks=1000)

#FIXME need to turn tfalc.TF1, .TF2, etc files to .tfalc
assert 0

if not os.path.exists(tfastatfile):
    ap.parallel_lc_statistics(lcdirectory, tfalcglob, fovcatalog,
                              tfalcrequired=True,
                              fovcatcols=(0,9), # objectid, magcol from fovcat
                              fovcatmaglabel='r',
                              outfile=tfastatfile,
                              nworkers=20,
                              workerntasks=500,
                              rmcols=[14,19,24],
                              epcols=[27,28,29],
                              tfcols=[30,31,32],
                              rfcols=None,
                              correctioncoeffs=None,
                              sigclip=5.0)

else:
    print('already made TFA LC stats')

ap.plot_stats_file(tfastatfile, statspath, field, binned=False, logy=True,
                   logx=False, correctmagsafter=None, rangex=(5.9,13.1))

#FIXME: still need to collect into single .tfalc files for all apertures

# Step 13: bin LCs to 10 minutes.
ap.parallel_bin_lightcurves(lcdirectory,
                            epdlc_glob=epdlcglob,
                            binsizes=[600,3600],
                            lcexts=('epdlc',
                                    'tfalc.TF1','tfalc.TF2','tfalc.TF3'),
                            jdcol=0,
                            lcmagcols=([27,28,29],[30,],[30,],[30,]),
                            nworkers=20,
                            workerntasks=1000)

# # Step 14: run plot_stats_file to make RMS vs. mag plots for all unbinned and
# # binned LCs.
# 
# ap.plot_stats_file(tfastatfile, statspath, field, binned=False, logy=True,
#                    logx=False, correctmagsafter=None, rangex=(5.9,13.1))


# run parallel_binnedlc_statistics to collect stats for the binned LCs.
assert 0 #what is appropriate glob for binned LCs?

ap.parallel_binnedlc_statistics(lcdirectory, lcglob, fovcatalog,
                                fovcatcols=(0,9), fovcatmaglabel='r',
                                corrmagsource=None, corrmag_idcol=0,
                                corrmag_magcols=[122,123,124], outfile=None,
                                nworkers=20, workerntasks=500, sigclip=5.0)

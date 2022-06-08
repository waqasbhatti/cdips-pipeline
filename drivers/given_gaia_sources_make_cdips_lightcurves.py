"""
Alternate mode for CDIPS reductions.

The default mode (in TESS_reduction.py) is to go down to G_RP<16 for the
cluster list, and G_RP<13 for field stars.  Often though, we have extra stars
for which we want to use the same machinery to get light curves from the
**already existing** difference images.
"""
###########
# imports #
###########
import os, shutil, uuid
from datetime import datetime
from glob import glob
import pandas as pd, numpy as np
from numpy import array as nparr
from copy import deepcopy
from astropy.io import fits

# non-standard: https://github.com/christopherburke/tess-point
from tess_stars2px import tess_stars2px_function_entry

from cdips.utils.gaiaqueries import gaia2read_given_df

# non-standard: cdips-pipeline imports
import aperturephot as ap
import autoimagesub as ais
import imagesubphot as ism
import lcutils as lcu
from TESS_reduction import main, run_detrending

######################
# variable arguments #
######################

# unique string identifying the reduction, to be used in directories, projid
# strings, etc.  lives at /cdips-pipeline/drivers/targetlists/{reduc_id}.csv
# needs at minimum the key "dr2_source_id".  optionally "ra" and "dec".
reduc_id = 'Meingast_2021_allstars'

######################
# validate arguments #
######################

# how far to search for stars on silicon.
MAX_SECTOR = 26

# comma-separated CSV file contanining at minimum a column named dr2_source_id,
# with the Gaia DR2 source identifiers.  If 'ra' and 'dec' columns are present,
# they are assumed to be correct.
targetlistcsv = os.path.join('targetlists', f'{reduc_id}.csv')
assert os.path.exists(targetlistcsv)
targetdf = pd.read_csv(targetlistcsv)
assert 'dr2_source_id' in targetdf

###############
# main script #
###############

#
# get ra/dec for the stars if needed, from the source_ids
#
srcpath = os.path.join('targetlists', f'{reduc_id}_sources_only.csv')
dstpath = os.path.join('targetlists', f'{reduc_id}_gaia2read.csv')

if not ( ('ra' in targetdf) and ('dec' in targetdf)):

    targetdf['dr2_source_id'].to_csv(srcpath, index=False, header=False)

    if not os.path.exists(dstpath):
        gaia2readcmd = f"gaia2read --header --extra --idfile {srcpath} --out {dstpath}"
        print(f'Beginning {gaia2readcmd}')
        returncode = os.system(gaia2readcmd)
        if returncode != 0: raise AssertionError('gaia2read cmd failed!!')
        print(f'Ran {gaia2readcmd}')
    else:
        print(f'Found {dstpath}')

    targetdf = pd.read_csv(dstpath, delim_whitespace=True)
    targetdf = targetdf.rename({
        '#Gaia-ID[1]':'dr2_source_id',
        'RA[deg][2]':'ra',
        'Dec[deg][3]':'dec',
        'phot_g_mean_mag[20]':'phot_g_mean_mag',
        'phot_bp_mean_mag[25]':'phot_bp_mean_mag',
        'phot_rp_mean_mag[30]':'phot_rp_mean_mag',
    }, axis='columns')

#
# figure out what data are *expected* for these stars using tess-point.
#
tesspointpath = os.path.join('targetlists', f'{reduc_id}_tesspoint.csv')

ra = nparr(targetdf.ra)
dec = nparr(targetdf.dec)
starids = nparr(targetdf.dr2_source_id).astype(str)

if not os.path.exists(tesspointpath):

    print(f'Beginning tess-point call...')
    (outID, outEclipLong, outEclipLat,
     sector, cam, ccd,
     colpix, rowpix, scinfo ) = (
         tess_stars2px_function_entry(starids, ra, dec)
     )

    tpdf = pd.DataFrame({
        'dr2_source_id': nparr(outID).astype(str),
        'sector': sector,
        'cam': cam,
        'ccd': ccd,
        'colpix': colpix,
        'rowpix': rowpix
    })

    tpdf = tpdf[tpdf.sector <= MAX_SECTOR]

    tpdf.to_csv(tesspointpath, index=False)

tpdf = pd.read_csv(tesspointpath)
tpdf['dr2_source_id'] = tpdf['dr2_source_id'].astype(str)

#
# figure out what data are *available* (as fully-formatted CDIPS light curves)
# for the passed stars.
#
from cdips.utils.lcutils import find_cdips_lc_paths

availpath = os.path.join('targetlists', f'{reduc_id}_tesspoint_existing.csv')

if not os.path.exists(availpath):

    lcpathdict = {}
    print('Beginning already existing CDIPS LC search...')
    for dr2_source_id in starids:
        lcpaths = find_cdips_lc_paths(dr2_source_id, raise_error=False)
        if lcpaths is None:
            lcpathdict[dr2_source_id] = [-1]
        else:
            lcpathdict[dr2_source_id] = lcpaths
    print('Completed already existing CDIPS LC search...')

    _lcpaths, hasmatch = [], []
    for ix, r in tpdf.iterrows():
        dr2_source_id = r['dr2_source_id']
        sector = r['sector']

        match = -1
        _hasmatch = 0
        for l in lcpathdict[dr2_source_id]:
            if isinstance(l, int):
                continue
            if f'sector-{sector}' in l:
                match = l
                _hasmatch = 1

        _lcpaths.append(match)
        hasmatch.append(_hasmatch)

    df = deepcopy(tpdf)
    df['hasmatch'] = hasmatch
    df['lcpath'] = _lcpaths

    df.to_csv(availpath, index=False)

df = pd.read_csv(availpath)

#
# re-calculate reduction/project identifier numbers
#
from cdips.utils.collect_cdips_lightcurves import (
    given_sector_cam_ccd_get_projid
)
df['projid'] = [
    given_sector_cam_ccd_get_projid(sec,cam,ccd)
    for sec,cam,ccd in zip(
        df.sector, df.cam, df.ccd
    )
]

#
# these are the stars that do not already have light curves, but need them.
#
sel = ~df.hasmatch.astype(bool)
sdf = df[sel]
N_needed = len(sdf)

#
# iterate over projids (i.e., sector/camera/ccd combinations)
#
uprojids = np.unique(sdf.projid)

returncodes = {}

outdirbase = f"/nfs/phtess2/ar0/TESS/REREDUC/{reduc_id}"
outlcdirbase = f"/nfs/phtess2/ar0/TESS/REREDUC/{reduc_id}_LC"
outcsvpath = os.path.join(outdirbase, 'reduc_status.log')

if not os.path.exists(outcsvpath):

    for ix, projid in enumerate(uprojids):

        # these are the stars we need light curves for on this sec/cam/ccd
        _sel = (sdf.projid == projid)
        _sdf = sdf[_sel]
        N_desired = len(_sdf)

        starttime = datetime.utcnow().isoformat()
        print(10*'-')
        print(f'{starttime}: {ix}/{len(uprojids)}: Begin projid {projid} with '
              f'{N_desired} stars.')
        print(_sdf)
        print(10*'-')

        sector = np.unique(_sdf.sector)
        assert len(sector) == 1
        sector = int(sector)
        cam = int(np.unique(_sdf.cam))
        ccd = int(np.unique(_sdf.ccd))

        ################
        # define paths #
        ################
        LOCAL_IMGBASE = "/nfs/phtess2/ar0/TESS/FFI/RED_IMGSUB/FULL"

        sectordir = LOCAL_IMGBASE + f"/s{str(sector).zfill(4)}/"
        fitsdir = sectordir + f"RED_{cam}-{ccd}-{projid}_ISP/"
        LOCAL_GLOBPATTERN = f'tess?????????????-s{str(sector).zfill(4)}-{cam}-{ccd}-*_cal_img_bkgdsub.fits'
        fitsglob = LOCAL_GLOBPATTERN
        lcbase = "/nfs/phtess2/ar0/TESS/FFI/LC/FULL"
        lcsector =  lcbase + f"/s{str(sector).zfill(4)}/"
        lcdir = lcsector + f"ISP_{cam}-{ccd}-{projid}/"

        #################################
        # reduction-specific parameters #
        #################################
        tuneparameters = False
        nworkers = 32
        aperturelist = "1.0:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0"
        epdsmooth = 11
        epdsigclip = 10000
        photdisjointradius = 0.5
        anetfluxthreshold = 1000
        anettweak = 6
        anetradius = 30
        initccdextent = "0:2048,0:2048"
        kernelspec="i/2;d=3/2"
        cluster_faint_Rp_mag = 17.5
        field_faint_Rp_mag = 13
        fiphotfluxthreshold = 300
        photreffluxthreshold = 300
        extractsources = 0
        binlightcurves = 0
        translateimages = 1
        reversesubtract = 1
        skipepd = 1
        zeropoint = 11.82
        useimagenotfistar = True

        np.random.seed(42)

        from TESS_reduction import given_fits_list_get_gain_exptime_ra_dec
        fits_list = np.sort(glob(os.path.join(fitsdir, fitsglob)))
        assert len(fits_list) > 0
        ccdgain, exptime, ra_nom, dec_nom = (
            given_fits_list_get_gain_exptime_ra_dec(fits_list)
        )

        #
        # ACTUAL STEPS IN THIS SUB-PIPELINE:
        # make_fov_catalog
        # reform_gaia_fov_catalog
        # parallel_convsubfits_staticphot
        # dump_lightcurves_with_grcollect
        # run_detrending : the full wrapper
        #

        #
        # make_fov_catalog
        #
        outdir = f"/nfs/phtess2/ar0/TESS/REREDUC/{reduc_id}/s{str(sector).zfill(4)}-{cam}-{ccd}-{projid}/"
        outlcdir = f"/nfs/phtess2/ar0/TESS/REREDUC/{reduc_id}_LC/s{str(sector).zfill(4)}-{cam}-{ccd}-{projid}/"
        for d in [outdirbase, outdir, outlcdirbase, outlcdir]:
            if not os.path.exists(d):
                os.mkdir(d)
                print(f"Made {d}")

        _lcnames = [os.path.join(outlcdir, f'{s}_llc.fits') for s in _sdf.dr2_source_id]
        newlcexists = [os.path.exists(f) for f in _lcnames]
        if np.all(newlcexists):
            print(f'Found LCs {_lcnames}, continue.')
            returncode = 0
            returncodes[projid] = (returncode, N_desired)
            continue
        if np.any(newlcexists):
            print(f'Found {np.sum(newlcexists)}/{len(newlcexists)} of the '
                  f'LCs for {projid}:\n{_lcnames}, continue.')
            returncode = 2
            returncodes[projid] = (returncode, N_desired)
            continue

        catalog_file = os.path.join(
            outdir, f"{reduc_id}-s{str(sector).zfill(4)}-{cam}-{ccd}.catalog"
        )
        if not os.path.exists(catalog_file):

            missedsrcpath = os.path.join(
                'targetlists', f'{reduc_id}_missed_sources_only.csv'
            )
            _sdf['dr2_source_id'].to_csv(
                missedsrcpath, index=False, header=False
            )

            gaia2readcmd = (
                f"gaia2read --header --extra --xieta-coords {ra_nom},{dec_nom} "
                f"--idfile {missedsrcpath} --out {catalog_file}"
            )
            print(f'Beginning {gaia2readcmd}')
            returncode = os.system(gaia2readcmd)
            if returncode != 0: raise AssertionError('gaia2read cmd failed!!')
            print(f'Ran {gaia2readcmd}')
        else:
            print(f'Found {catalog_file}')

        #
        # reform_gaia_fov_catalog
        #
        reformed_cat_file = catalog_file.replace('.catalog', '.reformed_catalog')
        ap.reform_gaia_fov_catalog(catalog_file, reformed_cat_file)

        #
        # parallel_convsubfits_staticphot
        # do photometry on your subtracted frames to produce .iphot files.
        #

        # first get metadata.  then run the convsubfits forced aperture photometry.
        field = 's'+str(sector).zfill(4)
        photreftype = 'onenight'
        fieldinfo = {}
        fieldinfo['camera'] = cam
        fieldinfo['ccd'] = ccd
        fieldinfo['projectid'] = int(projid)
        fieldinfo['field'] = field

        photparams = {
            'ccdgain': ccdgain, 'ccdexptime': exptime, 'zeropoint': zeropoint
        }

        if len(glob(os.path.join(outdir, '*.iphot'))) < int(1e9):

            # run photometry on the combinedphotref and generate a cmrawphot file. this
            # produces the base photometry values that we'll be diffing from those
            # found in the difference images to get difference magnitudes.
            # if extractsources==False, a placeholder fistar file with SDK values is
            # nonetheless generated, to be used for photometric reference frame
            # statistical bookkeeping.

            # the output combinedphotref path (previously made!)
            photreffname = ('proj{projid}-{field}-cam{cam}-ccd{ccd}'
                            '-combinedphotref-{photreftype}.{fileext}')
            photreffglob = (f'*proj{projid}-{field}-cam{cam}-ccd{ccd}'
                            f'-combinedphotref-{photreftype}*')
            astromreffglob = (f'*proj{projid}-camera{cam}-ccd{ccd}'
                              '-astromref*')

            refdir = f'/nfs/phtess2/ar0/TESS/FFI/BASE/rereduce-reference-frames/{reduc_id}'
            oldrefdir = '/nfs/phtess2/ar0/TESS/FFI/BASE/reference-frames'
            if not os.path.exists(refdir): os.mkdir(refdir)

            combinedphotrefpath = os.path.join(
                refdir,
                photreffname.format(
                    projid=projid,
                    field=field,
                    cam=cam,
                    ccd=ccd,
                    photreftype=photreftype,
                    fileext='fits'
                )
            )

            # copy all the old photometric and astrometric reference files into our
            # new "rereduce-reference-frames" directory. back up the .cmrawphot and
            # .rawphot files from the original reduction to be doubly-redundant
            # (and because they will be used to concatenate).
            photreffiles = glob(os.path.join(oldrefdir, photreffglob))
            astromreffiles = glob(os.path.join(oldrefdir, astromreffglob))
            filelist = photreffiles + astromreffiles

            for _file in filelist:
                src = _file
                dst = os.path.join(refdir, os.path.basename(_file))
                if os.path.exists(dst):
                    os.remove(dst)
                    print(f'Cleaning out {dst}')
                shutil.copy(src, dst)
                print(f'Copy {src} -> {dst}')

            photpaths = glob( os.path.join(
                refdir, f'*proj{projid}-{field}-cam{cam}-ccd{ccd}*.*phot'
            ))
            assert len(photpaths) == 3
            for photpath in photpaths:
                src = photpath
                ext = os.path.basename(photpath).split('.')[-1]
                dst = photpath.replace(ext,ext+'-backup')
                assert os.path.exists(src)
                os.rename(src, dst)
                print(f'Rename {src}->{dst}')

            photref_reformedfovcatpath = reformed_cat_file

            hdu_list = fits.open(combinedphotrefpath)
            hdr = hdu_list[0].header

            ccdgain = hdr['GAINA'] # electrons/count, from CCD output A. (there are four!)
            ccdexptime = int(np.round(hdr['TELAPSE']*24*60*60)) # in seconds, 1800
            ra_nom = hdr['CRVAL1']  # RA at CRPIX1, CRPIX2. Roughly "camera boresight".
            dec_nom = hdr['CRVAL2'] # DEC at CRPIX1, CRPIX2
            hdu_list.close()

            astrometrydownsample = 2
            pixelerror = 0.3
            uniformize = 10

            cphotref_photometry = ism.photometry_on_combined_photref(
                combinedphotrefpath,
                photref_reformedfovcatpath,
                ra_nom,
                dec_nom,
                ccd,
                ccdgain=ccdgain,
                zeropoint=zeropoint,
                ccdexptime=ccdexptime,
                extractsources=extractsources,
                apertures=aperturelist,
                framewidth=None,
                searchradius=8.0,
                photreffluxthreshold=photreffluxthreshold,
                observatory='tess',
                useimagenotfistar=useimagenotfistar,
                astrometrydownsample=astrometrydownsample,
                pixelerror=pixelerror,
                uniformize=uniformize,
                reformed_cat_file=reformed_cat_file,
                projid=projid
            )

            if cphotref_photometry is None:
                # means: the only star(s) in this projid were not found to be on
                # silicon
                returncode = 1
                returncodes[projid] = (returncode, N_desired)
                continue

            subfitslist = glob(os.path.join(
                fitsdir, '[r|n]sub-????????-'+
                fitsglob.replace('.fits','-xtrns.fits'))
            )

            out = ais.parallel_convsubfits_staticphot(
                subfitslist, fitsdir=fitsdir, fitsglob=fitsglob,
                photreftype=photreftype, kernelspec=kernelspec,
                lcapertures=aperturelist, photdisjointradius=photdisjointradius,
                outdir=outdir, fieldinfo=fieldinfo, observatory='tess',
                nworkers=nworkers, maxworkertasks=1000, photparams=photparams,
                dorereduce=reduc_id)

            if out==42:
                raise AssertionError('fatal error in convsubfits_staticphot')

        else:
            print('found .iphot files. skipping their production.')

        #
        # dump_lightcurves_with_grcollect
        #
        iphotpattern = os.path.join(
            outdir, 'rsub-????????-'+fitsglob.replace('.fits','.iphot')
        )
        if len(glob(os.path.join(outlcdir, '*.grcollectilc'))) < int(1e9):
            ism.dump_lightcurves_with_grcollect(
                iphotpattern, outlcdir, '4g', lcextension='grcollectilc',
                objectidcol=3, observatory='tess', fitsdir=fitsdir)
        else:
            print('found grcollectilc files. skipping lc dump.')

        #
        # run_detrending : the full wrapper.  kind of obnoxious, because most of
        # the steps are for TFA, which I don't care about.  nonetheless, we are
        # bound by format homogeneity.
        #
        epdlcglob, tfalcglob = '*_llc.fits', '*_llc.fits'
        statsdir = os.path.dirname(lcdir) + '/stats_files/'
        for dirname in [lcdir, statsdir]:
            assert os.path.exists(dirname)
        epdstatfile = statsdir + 'camera' + str(cam) + '_ccd' + str(ccd) + '.epdstats'
        tfastatfile = statsdir + 'camera' + str(cam) + '_ccd' + str(ccd) + '.tfastats'
        vartoolstfastatfile = os.path.join(statsdir, 'vartools_tfa_stats.txt')

        #
        # hacky call to:
        # parallel_convert_grcollect_to_fits_lc
        # parallel_apply_barycenter_time_correction
        #
        run_detrending(epdstatfile, tfastatfile, vartoolstfastatfile, outlcdir,
                       epdlcglob, reformed_cat_file, statsdir, field, fitsdir,
                       fitsglob, cam, ccd, projid, epdsmooth=epdsmooth,
                       epdsigclip=epdsigclip, nworkers=nworkers,
                       binlightcurves=binlightcurves,
                       tfa_template_sigclip=5.0, tfa_epdlc_sigclip=5.0,
                       skipepd=skipepd, fixedtfatemplate=None,
                       nmax_flow_logic=N_needed, escapeafterbarycenter=1,
                       barycenterparallel=False)

        #
        # call TFA runner here, instead of in run_detrending, since all the work
        # has already been done
        #
        lcfiles = glob(os.path.join(outlcdir, '*_llc.fits'))
        tfalclist_path = lcu._make_tfa_lc_list(lcfiles, outlcdir)
        trendlisttfa_paths = glob(os.path.join(statsdir, 'trendlist_tfa_ap?.txt'))
        assert len(trendlisttfa_paths) == 3
        datestfa_path = os.path.join(statsdir, 'dates_tfa.txt')
        assert os.path.exists(datestfa_path)

        npixexclude = 20
        tfafromirm = skipepd

        lcu.run_tfa(tfalclist_path, trendlisttfa_paths, datestfa_path,
                    outlcdir, statsdir, nworkers=nworkers,
                    do_bls_ls_killharm=False, npixexclude=npixexclude,
                    tfafromirm=tfafromirm)

        lcu.parallel_merge_tfa_lcs(outlcdir, nworkers=nworkers)

        # HYPE!
        endtime = datetime.utcnow().isoformat()
        print(f'{endtime}: {ix}/{len(uprojids)}: '
              f'Finished projid {projid} with {len(_sdf)} stars.')

        returncode = 0
        returncodes[projid] = (returncode, N_desired)

    print(f'finished {reduc_id}!')
    outdf = pd.DataFrame(returncodes, index=['returncode','n_desired'])
    outdf.to_csv(outcsvpath, index=False)
    print(f'wrote final status to {outcsvpath}')

#
# assess whether the light curves were actually made
#
print(f'Found {outcsvpath}')

newlcpaths = []
hasmatches = []
for ix, r in sdf.iterrows():
    source_id = str(r['dr2_source_id'])
    sector = int(r['sector'])
    cam = int((r['cam']))
    ccd = int((r['ccd']))
    projid = int((r['projid']))

    lcdir = f"/nfs/phtess2/ar0/TESS/REREDUC/{reduc_id}_LC"

    trylcpath = os.path.join(
        lcdir,
        f"s{str(sector).zfill(4)}-{cam}-{ccd}-{projid}/{source_id}_llc.fits"
    )
    if os.path.exists(trylcpath):
        newlcpaths.append(trylcpath)
        hasmatches.append(1)
    else:
        newlcpaths.append(-1)
        hasmatches.append(0)

odf = deepcopy(sdf)
odf['lcpath'] = newlcpaths
odf['hasmatch'] = hasmatches

gdf = gaia2read_given_df(odf, outdirbase)

assert len(odf) == len(gdf)

odf = odf.reset_index()
odf['phot_g_mean_mag'] = gdf['phot_g_mean_mag']

finalcsvpath = os.path.join(outdirbase, 'reduc_status_newlcs.log')
odf = odf.drop('index', axis='columns')
odf.to_csv(finalcsvpath, index=False)
print(f'Wrote {finalcsvpath}')

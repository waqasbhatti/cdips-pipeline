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
import os
import pandas as pd, numpy as np
from numpy import array as nparr
from copy import deepcopy

# non-standard: https://github.com/christopherburke/tess-point
from tess_stars2px import tess_stars2px_function_entry

######################
# variable arguments #
######################

# unique string identifying the reduction, to be used in directories, projid
# strings, etc.
reduc_id = 'Meingast_2021_n100'

# how far to search for stars on silicon.
MAX_SECTOR = 26

######################
# validate arguments #
######################

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
        'Dec[deg][3]':'dec'
        'phot_g_mean_mag[20]','phot_g_mean_mag',
        'phot_bp_mean_mag[25]','phot_bp_mean_mag',
        'phot_rp_mean_mag[30]','phot_rp_mean_mag',
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

#
# these are the stars that do not already have light curves, but need them.
#
df = pd.read_csv(availpath)
sel = ~df.hasmatch.astype(bool)
sdf = df[sel]

import IPython; IPython.embed()

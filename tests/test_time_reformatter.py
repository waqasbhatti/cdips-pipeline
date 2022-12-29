"""
How long does it take to run rlm.reformat_headers, for a given set of test
light curves?

Date: Wed 28 Dec 2022 07:32:16 PM EST
Author: LGB

The answer for a period after the Cycle 2 reductions, and during some of the
Cycle 3 reductions, was 5-6 seconds per light curve.  This was because of
inefficiencies described in
https://github.com/lgbouma/cdips/commit/35d0c1cee57fdb022888ae63e53c0dc28f804aa1
and has since been improved by a factor of ~4x.
"""
from glob import glob
import os
from astrobase import imageutils as iu
from cdips.lcproc import detrend as dtr
from cdips.lcproc import reformat_lcs_for_mast as rlm
from datetime import datetime

symlinkdir = "/nfs/phtess2/ar0/TESS/REREDUC/test_data/"
outdirnew = "/nfs/phtess2/ar0/TESS/REREDUC/test_output/"
cam = 2
ccd = 3
sector = 18
projid = 1820
OC_MG_CAT_ver = 0.6
cdipsvnum = 1

sectordir = os.path.join(outdirnew, f'sector-{sector}')
if not os.path.exists(sectordir): os.mkdir(sectordir)

camccddir = os.path.join(sectordir, f'cam{cam}_ccd{ccd}')
if not os.path.exists(camccddir): os.mkdir(camccddir)

lcpaths = glob(os.path.join(
    symlinkdir,
    f'sector-{sector}', f'cam{cam}_ccd{ccd}', '*_llc.fits'
))

if len(lcpaths) > 0:

    projid = iu.get_header_keyword(lcpaths[0], 'PROJID')

    eigveclist, smooth_eigveclist, n_comp_df = dtr.prepare_pca(
        cam, ccd, sector, projid
    )

    rlm.reformat_headers(lcpaths, camccddir, sector,
                         cdipsvnum, OC_MG_CAT_ver,
                         eigveclist=eigveclist,
                         smooth_eigveclist=smooth_eigveclist,
                         n_comp_df=n_comp_df)

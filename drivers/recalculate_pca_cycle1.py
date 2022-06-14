"""
The cycle 1 CDIPS light curves were created with default PCA detrending
parameters that typically overfit stellar rotation signals.  This script
recalculates the PCA light curves for all of cycle 1, and writes the
corresponding v0002 light curves to a dedicated directory.
"""
import os, shutil
from glob import glob
import multiprocessing as mp

from cdips.lcproc import trex_lc_to_mast_lc as tlml

def main():

    outdir = '/nfs/phtess3/ar0/TESS/PROJ/lbouma/CYCLE1PCAV2/LC'
    symlinkdir = '/nfs/phtess3/ar0/TESS/PROJ/lbouma/CYCLE1PCAV2/symlink'
    overwrite = 1 # NOTE: needed for secondary runs
    cams = [1,2,3,4]
    ccds = [1,2,3,4]
    OC_MG_CAT_ver = 0.6
    cdipsvnum = 2
    nworkers = mp.cpu_count()

    #for sector in range(1,14):
    for sector in [6]:

        lcpaths = glob(os.path.join(outdir, 'sector-{}'.format(sector),
                                    'cam?_ccd?', 'hlsp*.fits'))

        # turn cdips-pipeline light curves to HLSP light curves
        if len(lcpaths) == 0 or overwrite:
            tlml.trex_lc_to_mast_lc(sectors=[sector], cams=cams, ccds=ccds,
                                    make_symlinks=1, reformat_lcs=1,
                                    OC_MG_CAT_ver=OC_MG_CAT_ver,
                                    cdipsvnum=cdipsvnum, outdir=outdir,
                                    symlinkdir=symlinkdir)
        else:
            print('found {} HLSP LCs; wont reformat'.format(len(lcpaths)))

if __name__ == "__main__":
    main()

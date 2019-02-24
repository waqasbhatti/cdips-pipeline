# -*- coding: utf-8 -*-
'''
programmatically generate all the bash files in /scripts/tess_tuning_scripts/
that are executed to explore the parameter space of the photometric reduction.

Usage:

    $ cd $PIPETREXDIR/scripts
    $ python make_tess_tuning_scripts.py
'''
from __future__ import division, print_function

import os
from parse import parse, search
import re

import numpy as np

def set_line(matchstr, value, templatelines):
    """
    Given templatelines, and a string to match (matchstr), set the line value
    using "value" and the return templatelines.

    Args:
        matchstr (str): e.g., "epdsmooth" if you want to change epdsmooth.
        value (int or str): e.g., "7" if you want to set epdsmooth to 7.
        templatelines (list of strings): read in from a TEMPLATE shell file,
        e.g., `TEMPLATE_single_tess_ccd.sh`.
    """

    matchstr+='='

    for ix, l in enumerate(templatelines):
        if l.startswith(matchstr):
            lind, line = ix, l
        else:
            pass

    currentval = parse(matchstr+'{}', line)[0]
    try:
        line = line.replace(currentval, str(value))
    except Exception as e:
        print('caught error: {}'.format(e))
        import IPython; IPython.embed()

    templatelines[lind] = line

    return templatelines


def main():

    fname = 'tess_tuning_scripts/TEMPLATE_single_tess_ccd.sh'

    with open(fname,'r') as fhandle:
        templatelines = fhandle.readlines()

    # default kernelspec: i/3;d=3/2
    # change delta kernel half-size
    # change delta kernel order
    # change identity order.
    kernellist = [
        "\"i/3;d=3/2\"",
        "\"i/3;d=1/2\"",
        "\"i/3;d=2/2\"",
        "\"i/3;d=4/2\"",
        "\"i/3;d=5/2\"",
        "\"i/3;d=3/1\"",
        "\"i/3;d=3/3\"",
        "\"i/3;d=3/4\"",
        "\"i/3;d=3/5\"",
        "\"i/1;d=3/2\"",
        "\"i/2;d=3/2\"",
        "\"i/4;d=3/2\"",
        "\"i/5;d=3/2\""
    ]

    apertureslist = [
        "\"1.0:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0\""
    ]

    n_kernels_to_run = len(kernellist)

    # all TUNE reductions for camN, ccdN (3 different cam/ccd pairs). cam 2,
    # ccd 2 has benefit of wasp-4 being on chip!
    camccdnums = [2]
    camnumlist = camccdnums*n_kernels_to_run
    ccdnumlist = camccdnums*n_kernels_to_run

    n_runs = len(camnumlist)

    # make the paramdict for the first 3*13 kernel list runs
    projid = 1194
    varyparamdict = {}

    cam, ccd = 2, 2
    fixedtfatemplate = "\"/nfs/phtess1/ar1/TESS/FFI/LC/FULL/s0002/STATS_FILES/ISP_2-2-1096_stats_files\""
    for kernel in kernellist:
        varyparamdict[projid] = {
            'camnum':cam,
            'ccdnum':ccd,
            'projectid':projid,        # increment this whenever new things go to PSQL database.
            'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':'false',
            'nworkers':32,
            'aperturelist':apertureslist[0],
            'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
            'epdsigclip':10000,
            'photdisjointradius':2,
            'anetfluxthreshold':50000,
            'anettweak':6,
            'anetradius':30,
            'initccdextent':"\"0:2048,0:2048\"",
            'kernelspec':kernel,
            'catalog_faintrmag':14,      ## catalog_faintrmag=16
            'fiphotfluxthreshold':1000,  ## fiphotfluxthreshold=300
            'photreffluxthreshold':1000, ## photreffluxthreshold=300
            'extractsources':0,
            'binlightcurves':0,
            'translateimages':1,
            'reversesubtract':1,
            'skipepd':0,
            'fixedtfatemplate':fixedtfatemplate
        }
        projid += 1

    # then cam1, ccd2
    cam, ccd = 1, 2
    fixedtfatemplate = "\"/nfs/phtess1/ar1/TESS/FFI/LC/FULL/s0002/STATS_FILES/ISP_1-2-1163_stats_files\""
    for kernel in kernellist:
        varyparamdict[projid] = {
            'camnum':cam,
            'ccdnum':ccd,
            'projectid':projid,        # increment this whenever new things go to PSQL database.
            'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':'false',
            'nworkers':32,
            'aperturelist':apertureslist[0],
            'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
            'epdsigclip':10000,
            'photdisjointradius':2,
            'anetfluxthreshold':50000,
            'anettweak':6,
            'anetradius':30,
            'initccdextent':"\"0:2048,0:2048\"",
            'kernelspec':kernel,
            'catalog_faintrmag':14,      ## catalog_faintrmag=16
            'fiphotfluxthreshold':1000,  ## fiphotfluxthreshold=300
            'photreffluxthreshold':1000, ## photreffluxthreshold=300
            'extractsources':0,
            'binlightcurves':0,
            'translateimages':1,
            'reversesubtract':1,
            'skipepd':0,
            'fixedtfatemplate':fixedtfatemplate
        }
        projid += 1

    projids = np.sort(list(varyparamdict.keys()))

    for projid in projids:

        for key in list(varyparamdict[projid].keys()):

            templatelines = set_line(key, varyparamdict[projid][key], templatelines)

        outname = 'tess_tuning_scripts/projid_{:s}.sh'.format(repr(projid))

        with open(outname,'w') as fhandle:
            fhandle.writelines(templatelines)

        # make the output scripts easily executable
        CHMODCMD = 'chmod +x {path}'
        cmdtorun = CHMODCMD.format(path=outname)
        returncode = os.system(cmdtorun)
        if returncode == 0:
            print('wrote {} and converted to executable'.
                  format(outname))
        else:
            print('ERR! wrote {} but was not able to convert to executable'.
                  format(outname))


if __name__=="__main__":
    main()

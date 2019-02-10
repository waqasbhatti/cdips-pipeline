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

    # default: "i/5;b/5;d=3/2"
    # just do identity term
    # change background order
    # change delta kernel order
    # change delta kernel half-size
    # change "identity" term order
    kernellist = [
        "\"i/3\"",
        "\"i/3;d=3/2\"",
        "\"i/3;b/3;d=3/2\"",
        "\"i/3;b/5;d=3/2\"",
        "\"i/5;b/3;d=5/5\"",
        "\"i/5;b/3;d=5/3\"",
        "\"i/5;b/3;d=5/2\"",
        "\"i/3;b/3\"",
        "\"i/3;b/3;d=5/2\"",
        "\"i/3;b/3;d=3/2\"",
        "\"i/5;b/3;d=3/2\"",
        "\"i/3;b/3;d=3/2\"",
        "\"i/2;b/3;d=3/2\""
    ]
    nobkgkernellist = [
        "\"i/3;d=3/2\"", # best from 2019/01/30
        "\"i/1;d=3/2\"",
        "\"i/2;d=3/2\"",
        "\"i/5;d=3/2\"",
        "\"i/3;d=5/2\"",
        "\"i/2;d=3/3\"",
        "\"i/2;d=2/2\"",
        "\"i/3;d=3/4\"",
        "\"i/3;d=3/3\""
    ]

    # top is default
    # change the inner and outer boundary of background annulus
    # vary the aperture sizes 
    apertureslist = [
        "\"0.71:7.0:6.0,1.41:7.0:6.0,2.82:7.0:6.0\"",
        "\"0.71:5.0:6.0,1.41:5.0:6.0,2.82:5.0:6.0\"",
        "\"0.71:5.0:4.0,1.41:5.0:4.0,2.82:5.0:4.0\"",
        "\"1:7.0:6.0,2:7.0:6.0,3:7.0:6.0\"",
        "\"1.7:7.0:6.0,3.5:7.0:6.0,4.5:7.0:6.0\""
    ]

    n_kernels_aps = len(kernellist) + len(apertureslist) # 18

    # all TUNE reductions for camN, ccdN (3 different cam/ccd pairs). cam 2,
    # ccd 2 has benefit of wasp-4 being on chip!
    camccdnums = [1,2,4]
    camnumlist = camccdnums*n_kernels_aps
    ccdnumlist = camccdnums*n_kernels_aps

    n_runs = len(camnumlist) # 54 runs, from the above

    # make the paramdict for the first 3*13 kernel list runs
    projid = 1089
    varyparamdict = {}
    for kernel in kernellist:
        for camccd in camccdnums:

            varyparamdict[projid] = {
                'camnum':camccd,
                'ccdnum':camccd,
                'projectid':projid,     # increment this whenever new things go to PSQL database.
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
                'translateimages':1
            }

            projid += 1

    # now ditto, for the subsequent 3*5 aperture size lists
    for aperturelist in apertureslist:
        for camccd in camccdnums:

            varyparamdict[projid] = {
                'camnum':camccd,
                'ccdnum':camccd,
                'projectid':projid,     # increment this whenever new things go to PSQL database.
                'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
                'tuneparameters':'false',
                'nworkers':32,
                'aperturelist':aperturelist,
                'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
                'epdsigclip':10000,
                'photdisjointradius':2,
                'anetfluxthreshold':50000,
                'anettweak':6,
                'anetradius':30,
                'initccdextent':"\"0:2048,0:2048\"",
                'kernelspec':"i/5;b/5;d=3/2",
                'catalog_faintrmag':14,      ## catalog_faintrmag=16
                'fiphotfluxthreshold':1000,  ## fiphotfluxthreshold=300
                'photreffluxthreshold':1000, ## photreffluxthreshold=300
                'extractsources':0,
                'binlightcurves':0,
                'translateimages':1
            }

            projid += 1

    # sector 2, cam1ccd2 -- for blanco 1 initial analysis. a full run down to
    # G_Rp=16. might take a while.
    varyparamdict[projid] = {
        'camnum':1,
        'ccdnum':2,
        'projectid':projid,     # increment this whenever new things go to PSQL database.
        'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
        'tuneparameters':'false',
        'nworkers':32,
        'aperturelist':apertureslist[3],
        'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
        'epdsigclip':10000,
        'photdisjointradius':2,
        'anetfluxthreshold':50000,
        'anettweak':6,
        'anetradius':30,
        'initccdextent':"\"0:2048,0:2048\"",
        'kernelspec':"i/3;d=3/2", # best from 20190130
        'catalog_faintrmag':16,
        'fiphotfluxthreshold':500,
        'photreffluxthreshold':500,
        'extractsources':0,
        'binlightcurves':0,
        'translateimages':1
    }

    projid += 1


    # 2019/01/30...
    # NOTE: BELOW WAS AN ERROR THAT COST A DAY OR SO OF COMPUTE TIME ...
    # SPOT IT ...
    # now ditto, for cam2ccd2, no bkg kernels. omit the one that was already
    # done!
    for kernel in nobkgkernellist[1:]:

        varyparamdict[projid] = {
            'camnum':2,
            'ccdnum':2,
            'projectid':projid,     # increment this whenever new things go to PSQL database.
            'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':'false',
            'nworkers':32,
            'aperturelist':aperturelist,
            'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
            'epdsigclip':10000,
            'photdisjointradius':2,
            'anetfluxthreshold':50000,
            'anettweak':6,
            'anetradius':30,
            'initccdextent':"\"0:2048,0:2048\"",
            'kernelspec':"i/5;b/5;d=3/2", # NOTE ERROR!
            'catalog_faintrmag':14,      ## catalog_faintrmag=16
            'fiphotfluxthreshold':1000,  ## fiphotfluxthreshold=300
            'photreffluxthreshold':1000, ## photreffluxthreshold=300
            'extractsources':0,
            'binlightcurves':0,
            'translateimages':1
        }

        projid += 1

    # 2019/02/02...
    # test out nsub. timing test for cam1ccd2 full (for pending kernel tests.)
    for tuneparam, faintmag in zip(['true','false'],[12,14]):
        varyparamdict[projid] = {
            'camnum':1,
            'ccdnum':2,
            'projectid':projid,     # increment this whenever new things go to PSQL database.
            'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':tuneparam,
            'nworkers':32,
            'aperturelist':apertureslist[0],
            'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
            'epdsigclip':10000,
            'photdisjointradius':2,
            'anetfluxthreshold':50000,
            'anettweak':6,
            'anetradius':30,
            'initccdextent':"\"0:2048,0:2048\"",
            'kernelspec':"i/3;d=3/2", # best from 20190130
            'catalog_faintrmag':faintmag,
            'fiphotfluxthreshold':500,
            'photreffluxthreshold':500,
            'extractsources':0,
            'binlightcurves':0,
            'translateimages':1,
            'reversesubtract':0
        }

        projid += 1

    # cam2ccd2, no bkg kernels. omit the one that was already done!
    for kernel in nobkgkernellist[1:]:
        varyparamdict[projid] = {
            'camnum':2,
            'ccdnum':2,
            'projectid':projid,     # increment this whenever new things go to PSQL database.
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
            'reversesubtract':1 # for cam2ccd2, everything else is rsub
        }

        projid += 1

    # cam1ccd2, full kernel list. 13+9=22 kernels. at 2.5hr each, runtime is
    # 2.5 clock days.
    for kernel in list(kernellist + nobkgkernellist):
        varyparamdict[projid] = {
            'camnum':1,
            'ccdnum':2,
            'projectid':projid,     # increment this whenever new things go to PSQL database.
            'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':'false',
            'nworkers':32,
            'aperturelist':"\"1:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0\"",
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
            'reversesubtract':1 # for cam2ccd2, everything else is rsub
        }

        projid += 1
    # note: rsub on above was an error. but c'est la vie.

    #FIXME
    projid = 1187 # MANUAL, bc i did some manual projids

    # cam3ccd1
    for kernel in list(kernellist + nobkgkernellist):
        varyparamdict[projid] = {
            'camnum':4,
            'ccdnum':3, #southern beehive
            'projectid':projid,     # increment this whenever new things go to PSQL database.
            'sector':"'s0001'",        # match SPOC syntax, zfill to 4.
            'tuneparameters':'false',
            'nworkers':32,
            'aperturelist':"\"1:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0\"",
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
            'reversesubtract':0
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

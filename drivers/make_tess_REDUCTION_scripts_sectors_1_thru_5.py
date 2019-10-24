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

    # start at projid 1300 for sector 1-5 reduction.
    projid = 1300
    varyparamdict = {}

    for snum in [2,3,4,5,1]:
        for cam in range(1,5):
            for ccd in range(1,5):
                varyparamdict[projid] = {
                    'camnum':cam,
                    'ccdnum':ccd,
                    'projectid':projid,        # increment this whenever new things go to PSQL database.
                    'sector':"'s000{}'".format(snum), # match SPOC syntax, zfill to 4.
                    'tuneparameters':'false',
                    'nworkers':32,
                    'aperturelist':"\"1.0:7.0:6.0,1.5:7.0:6.0,2.25:7.0:6.0\"",
                    'epdsmooth':11,            # not used
                    'epdsigclip':10000,        # ditto
                    'photdisjointradius':2,
                    'anetfluxthreshold':50000,
                    'anettweak':6,
                    'anetradius':30,
                    'initccdextent':"\"0:2048,0:2048\"",
                    'kernelspec':"\"i/2;d=3/2\"",
                    'cluster_faint_Rp_mag':16,
                    'field_faint_Rp_mag':13,
                    'fiphotfluxthreshold':300,
                    'photreffluxthreshold':300,
                    'extractsources':0,
                    'binlightcurves':0,
                    'translateimages':1,
                    'reversesubtract':1,
                    'skipepd':1
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

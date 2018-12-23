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

    # TODO: iterate over a list of parameters to generate a large number of
    # shell scripts, once the metrics for assessment are clarified.
    projectid = 1002

    paramdict = {
        'camnum':4,
        'ccdnum':4,
        'projectid':projectid,     # increment this whenever new things go to PSQL database.
        'sector':"'s0002'",        # match SPOC syntax, zfill to 4.
        'tuneparameters':'true',
        'nworkers':20,
        'aperturelist':"\"1.45:7.0:6.0,2.2:7.0:6.0,2.95:7.0:6.0\"",
        'epdsmooth':11,            # 11*30min = 5.5 hr median smooth in EPD pre-processing.
        'epdsigclip':10,
        'photdisjointradius':2,
        'anetfluxthreshold':50000,
        'anettweak':6,
        'anetradius':30,
        'initccdextent':"\"0:2048,0:2048\"",
        'kernelspec':"\"b/4;i/4;d=5/2\"",
        'catalog_faintrmag':13,      ## catalog_faintrmag=16
        'fiphotfluxthreshold':3000,  ## fiphotfluxthreshold=300
        'photreffluxthreshold':3000, ## photreffluxthreshold=300
        'extractsources':0,
        'binlightcurves':0
    }

    for key in list(paramdict.keys()):

        templatelines = set_line(key, paramdict[key], templatelines)

    outname = 'tess_tuning_scripts/projid_{:s}.sh'.format(repr(projectid))

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

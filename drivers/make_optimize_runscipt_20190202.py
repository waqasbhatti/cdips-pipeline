import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sequential_projids.sh & #
####################################

# remaining kernel experiments for cam2ccd2:
# * a bunch of "no background" kernels, similar as above.

# full kernel experiments for cam1ccd2:
# * default: "i/5;b/5;d=3/2"
# * just do identity term
# * change background order
# * change delta kernel order
# * change delta kernel half-size
# * change "identity" term order
# * a bunch of "no background" kernels, similar as above.

"""
)

projidlines = []
for projid in range(1154,1183):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/20190202_optimize_sequential_projids_full.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

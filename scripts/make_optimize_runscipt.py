import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sequential_projids.sh & #
####################################

# kernel experiments:
# * default: "i/5;b/5;d=3/2"
# * just do identity term
# * change background order
# * change delta kernel order
# * change delta kernel half-size
# * change "identity" term order

# aperture experiments:
# * default: "0.71:7.0:6.0,1.41:7.0:6.0,2.82:7.0:6.0"
# * change the inner and outer boundary of background annulus
# * vary the aperture sizes

"""
)

projidlines = []
for projid in range(1089,1143):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/20190128_optimize_sequential_projids_full.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

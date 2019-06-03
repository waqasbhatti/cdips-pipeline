import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sector7_reduction.sh & #
####################################

"""
)

projidlines = []
for projid in range(1516,1532):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/sector7_reduction.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

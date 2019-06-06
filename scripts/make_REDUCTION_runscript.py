import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sector9_reduction.sh & #
####################################

"""
)

projidlines = []
for projid in range(1548,1564):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/sector9_reduction.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

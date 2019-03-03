import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./REDUCTION_sectors_1_thru_5_projids_1300_thru_1380.sh & #
####################################

"""
)

projidlines = []
for projid in range(1300,1380):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/REDUCTION_sectors_1_thru_5_projids_1300_thru_1380.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

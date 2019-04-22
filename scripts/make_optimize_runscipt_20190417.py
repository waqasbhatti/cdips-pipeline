import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sequential_projids.sh & #
####################################

# sector 6, first galactic field reduction

"""
)

projidlines = []
for projid in range(1500,1517):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/sector6_reduction.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

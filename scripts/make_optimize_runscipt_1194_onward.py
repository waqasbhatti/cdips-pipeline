import os

outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sequential_projids.sh & #
####################################

# default kernelspec: i/3;d=3/2
# change delta kernel half-size
# change delta kernel order
# change identity order.
"""
)

projidlines = []
for projid in range(1194,1220):
    projidlines.append(
        '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
        format(projid)
    )

outname = 'tess_tuning_scripts/optimize_sequential_projids_1194_thru_1220.sh'
with open(outname, 'w') as f:
    f.writelines(outlines)
    f.writelines(projidlines)

print('wrote {}'.format(outname))

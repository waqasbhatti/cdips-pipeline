import os

sectors = list(range(14,27)) # [10,11,12,13]

startid = 1750
projidtuples = [ (startid+ix*16, startid+(ix+1)*16)  for ix in
                 range(0,len(sectors))  ]

for ix, sector in enumerate(sectors):

    outlines = (
"""#!/usr/bin/env bash

####################################
# USAGE: ./sector{}_reduction.sh & #
####################################

""".format(sector)
    )

    projidlines = []
    for projid in range(projidtuples[ix][0],projidtuples[ix][1]):
        projidlines.append(
            '( source activate trex_37; source projid_{}.sh; wait ) & wait\n'.
            format(projid)
        )

    outname = 'tess_tuning_scripts/sector{}_reduction.sh'.format(sector)
    with open(outname, 'w') as f:
        f.writelines(outlines)
        f.writelines(projidlines)

    print('wrote {}'.format(outname))

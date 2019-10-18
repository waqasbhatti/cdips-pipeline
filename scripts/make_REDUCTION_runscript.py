import os

sectors = [10,11,12,13]

projidtuples = [(1564,1580),
                (1564+16,1580+16),
                (1564+32,1580+32),
                (1564+48,1580+48)
               ]

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

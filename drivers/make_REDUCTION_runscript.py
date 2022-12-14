import os

# cycle 4
sectors = list(range(40, 56))

startid = 4000
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

    outname = f'tess_tuning_scripts/sector{sector}_reduction.sh'
    with open(outname, 'w') as f:
        f.writelines(outlines)
        f.writelines(projidlines)

    print(f'wrote {outname}')

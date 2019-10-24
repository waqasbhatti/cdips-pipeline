import imageutils as iu
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from glob import glob
import os

cam = 2
ccd = 3

caldir = '/nfs/phtess2/ar0/TESS/FFI/CAL/sector-10/'
ffiglob = 'tess*-s0010-2-3-*_ffic.fits'

useglob = os.path.join(caldir, ffiglob)
fitsfiles = glob(useglob)

ds = []
for ix, f in enumerate(fitsfiles):
    print('{}/{}'.format(ix, len(fitsfiles)))
    d = iu.get_header_keyword_list(f, ['DQUALITY', 'A_DMAX', 'B_DMAX'], ext=1)
    ds.append(d)

df = pd.DataFrame(ds)

df['fitsfile'] = fitsfiles

df = df.sort_values(by='fitsfile')

df.to_csv('s0010-2-3_quality_info.csv', index=False)


f, ax = plt.subplots(figsize=(6,3))
ax.scatter(range(len(df)), df['DQUALITY'])
ax.set_ylabel('dquality')
f.savefig('dquality_flags_s0010-2-3.png', dpi=300, bbox_inches='tight')


f, ax = plt.subplots(figsize=(6,3))
ax.scatter(range(len(df)), df['A_DMAX'])
ax.set_ylabel('A_DMAX')
f.savefig('A_DMAX_flags_s0010-2-3.png', dpi=300, bbox_inches='tight')

f, ax = plt.subplots(figsize=(6,3))
ax.scatter(range(len(df)), df['B_DMAX'])
ax.set_ylabel('B_DMAX')
f.savefig('B_DMAX_flags_s0010-2-3.png', dpi=300, bbox_inches='tight')

"""
How long does opening a fits file take?
"""
from datetime import datetime
from astropy.io import fits

lcpath = "/nfs/phtess2/ar0/TESS/REREDUC/test_data/sector-18/cam2_ccd3/348892464178559104_llc.fits"

# How long does it take to open/close a fits file?
t0 = datetime.utcnow().isoformat()
print(f'{t0}: beginning open for {lcpath}')

hdul = fits.open(lcpath)
primaryhdr, hdr, data = (
    hdul[0].header, hdul[1].header, hdul[1].data
)
hdul.close()

t1 = datetime.utcnow().isoformat()
print(f'{t1}: finish open & close for {lcpath}')


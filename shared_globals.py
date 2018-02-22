'''
These variables are shared across `cron_transfer_calibrate.py`, and
`autocal_direct.py`

LOCAL_IMGBASE and LOCAL_GLOBPATTERN define the local image directories.

FIELDCAT_DIR, _FOV, _BRIGHT, _FAINT, and _BANDS define the field catalog
location and properties

FIELD_REGEX and FIELD_CCDS defines the field string and CCDs
'''
import re

LOCAL_IMGBASE = '/nfs/phtess1/ar1/HATPI/HP0/RAW/',
LOCAL_GLOBPATTERN = '?-????????',
FIELDCAT_DIR = '/nfs/phtess1/ar1/HATPI/HP0/BASE/CAT',
FIELDCAT_FOV = 14.0,
FIELDCAT_BRIGHT = 0.0,
FIELDCAT_FAINT = 15.5,
FIELDCAT_BANDS = ['g,r,i','r,i,z'],
FIELD_REGEX = re.compile('^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$'),
FIELD_CCDS = [5,6,7,8]


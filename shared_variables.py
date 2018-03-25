'''
This file contains path variables that are set across the pipeline.
'''
import re, os

###############################################################
# The following variables are called in autocal_direct.py and #
# cron_transfer_calibrate                                     #
###############################################################

# LOCAL_IMGBASE and LOCAL_GLOBPATTERN define the local image directories.
LOCAL_IMGBASE = '/nfs/phtess1/ar1/TESS/SIMFFI/'
LOCAL_GLOBPATTERN = 'tess?????????????-1-1-0016_cal_img.fits'

FITS_TAIL = '.fits' # ".fits.fz" for HAT work.

# FIELDCAT_DIR, _FOV, _BRIGHT, _FAINT, and _BANDS define the field catalog
# location and properties
FIELDCAT_DIR = LOCAL_IMGBASE + 'CAT'
FIELDCAT_FOV = 14.0
FIELDCAT_BRIGHT = 0.0
FIELDCAT_FAINT = 15.5
FIELDCAT_BANDS = ['g,r,i','r,i,z']

#  FIELD_REGEX and FIELD_CCDS defines the field string and CCDs
FIELD_REGEX = re.compile('^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$')
FIELD_CCDS = [5,6,7,8]

#############################################
# The following are called in framecalib.py #
#############################################
# config files
TREXBASE = '/home/lbouma/proj/pipe-trex/' # user must set
CONFBASE = os.path.expanduser(TREXBASE+'/config-files')

# some paths
RAWPATH="/nfs/phtess1/ar1/HATPI/HP0/RAW/"
REDPATH="/nfs/phtess1/ar1/HATPI/HP0/RED/"
CALPATH='/nfs/phtess1/ar1/HATPI/HP0/BASE/CAL/'
TEMPPATH='/nfs/phtess1/ar1/HATPI/HP0/REDTEMP'

# HATPipepy script paths
HATPIPEPATH = '/home/wbhatti/HATpipebin_R3186/include/HATpipepy/Actions'
HATPIPEBIAS = os.path.join(HATPIPEPATH, 'do_masterbias.py')
HATPIPEFLAT = os.path.join(HATPIPEPATH, 'do_masterflat.py')
HATPIPEDARK = os.path.join(HATPIPEPATH, 'do_masterdark.py')
HATPIPEREDUCE = os.path.join(HATPIPEPATH, 'do_object.py')

REDUCTIONCONF = os.path.join(CONFBASE, 'RedCfg.cfg')
MASTERBIASCONF = os.path.join(CONFBASE, 'MasterBiasCfg.cfg')
MASTERFLATCONF = os.path.join(CONFBASE, 'MasterFlatCfg.cfg')
MASTERDARKCONF = os.path.join(CONFBASE, 'MasterDarkCfg.cfg')

# CCDS to look out for
CCDLIST = FIELD_CCDS

# FFMPEG path to make movies
FFMPEG = os.path.expanduser('~/bin/ffmpeg')

# FFMPEG commandline
FFMPEGCMD = ("{ffmpeg} -framerate {framerate} "
             "-pattern_type glob -i '{jpegglob}' "
             "-c:v libx264 -preset veryslow {outfile}")

"""
This file contains path variables that are set across the pipeline.
"""
import re, os

############################################################
# The following variables are called in autocal_direct.py, #
# cron_transfer_calibrate, and autoimagesub.py             #
############################################################

# LOCAL_IMGBASE and LOCAL_GLOBPATTERN define the local image directories.
LOCAL_IMGBASE = "/nfs/phtess1/ar1/HATPI/HP0/"
LOCAL_GLOBPATTERN = "tess?????????????-1-1-0016_cal_img.fits"

# main paths
RAWPATH=LOCAL_IMGBASE+"RAW/"
REDPATH=LOCAL_IMGBASE+"RED/projid12-G577-ccd8-sdssr/" #FIXME: this form of pathid is ugly
LCPATH=LOCAL_IMGBASE+"LC/"
CALPATH=LOCAL_IMGBASE+"BASE/CAL/"
CATPATH=LOCAL_IMGBASE+"BASE/CAT/"
TEMPPATH=LOCAL_IMGBASE+"REDTEMP/"

REFBASEDIR=LOCAL_IMGBASE+"BASE/reference-frames/"
REFINFO=os.path.join(REFBASEDIR,'refinfo.sqlite') # formerly TM-refinfo.sqlite

FRAMEINFOCACHEDIR=LOCAL_IMGBASE+"BASE/frameinfo-cache/"

FITS_TAIL = ".fits" # occasionally ".fits.fz" for HAT work. ??? do fitsh routines work on .fits.fz?

# define the field catalog location and properties
FIELDCAT_DIR = CATPATH
FIELDCAT_FOV = 14.0     # box half-width (in degrees)
FIELDCAT_BRIGHT = 0     # brightest magnitude
FIELDCAT_FAINT = 15.0   # faintest magnitude
FIELDCAT_BANDS = ["g,r,i","r,i,z"]

#  FIELD_REGEX and FIELD_CCDS defines the field string and CCDs
FIELD_REGEX = re.compile("^G(\d{2})(\d{2})([\+\-]\d{2})(\d{2})_(\w{3})$")
FIELD_CCDS = [5,6,7,8]

# CCDS to look out for
CCDLIST = FIELD_CCDS


#################################
# Paths called in framecalib.py #
#################################
# config files
TREXBASE = "/home/lbouma/proj/pipe-trex/" # user must set
CONFBASE = os.path.expanduser(TREXBASE+"/config-files")

# HATPipepy script paths
HATPIPEPATH = "/home/wbhatti/HATpipebin_R3186/include/HATpipepy/Actions"
HATPIPEBIAS = os.path.join(HATPIPEPATH, "do_masterbias.py")
HATPIPEFLAT = os.path.join(HATPIPEPATH, "do_masterflat.py")
HATPIPEDARK = os.path.join(HATPIPEPATH, "do_masterdark.py")
HATPIPEREDUCE = os.path.join(HATPIPEPATH, "do_object.py")

REDUCTIONCONF = os.path.join(CONFBASE, "RedCfg.cfg")
MASTERBIASCONF = os.path.join(CONFBASE, "MasterBiasCfg.cfg")
MASTERFLATCONF = os.path.join(CONFBASE, "MasterFlatCfg.cfg")
MASTERDARKCONF = os.path.join(CONFBASE, "MasterDarkCfg.cfg")

###################################
# Paths called in aperturephot.py #
###################################
# on LCO: "/nfs/lcohpsrv1/ar0/P/HP0/CAT/2MASS/2MASS_JH_AP/data"
TWOMASSPATH = "/nfs/phn12/ar0/H/CAT/2MASS/2MASS_JH_AP/data"
# on LCO: "/nfs/lcohpsrv1/ar0/P/HP0/CAT/UCAC4"
UCAC4PATH = "/nfs/phn12/ar0/H/HP0/CAT/UCAC4"

#####################################
# Postgres database credential info #
#####################################
PGPASSFILE = os.path.expanduser('~/.pgpass')
PGUSER = 'hpx'
PGDATABASE = 'hpx'
PGHOST = 'localhost'

# FFMPEG path to make movies
FFMPEG = os.path.expanduser("~/bin/ffmpeg")

# FFMPEG commandline
FFMPEGCMD = ("{ffmpeg} -framerate {framerate} "
             "-pattern_type glob -i '{jpegglob}' "
             "-c:v libx264 -preset veryslow {outfile}")

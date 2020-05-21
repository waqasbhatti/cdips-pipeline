# `cdips-pipeline`

[![DOI](https://zenodo.org/badge/202777810.svg)](https://zenodo.org/badge/latestdoi/202777810)

This is a time-series photometry pipeline that reduces images to light curves.
It is a port of a pipeline originally developed for the [HATPI
project](http://hatpi.org).

This pipeline has been used for the Cluster Difference Imaging Photometric
Survey (CDIPS) image-subtraction reductions.  The CDIPS project has made and is
making light curves for stars in open clusters, moving groups, etc.  It has [a
stand-alone repo](https://github.com/lgbouma/cdips).

In theory, if one wished to reproduce the CDIPS reductions, this pipeline would
be the place to start. More practically, the code-base should provide readers
and users with an entry-point into a particular set of programs that can be
used for photometric and time-series analyses.

We're releasing it for general reproducibility and method-sharing and
improvement reasons.  (Without the expectation of converting you into a user).

## 0. How is this structured?

The pipeline is a collection of python functions that can be called from a
single "driver script" to go through steps described in [the
paper](https://github.com/lgbouma/cdips/blob/master/paper_I/ms.pdf). For the
CDIPS-I reduction, the driver script is
[`drivers/TESS_reduction.py`](https://github.com/waqasbhatti/cdips-pipeline/blob/master/drivers/TESS_reduction.py).
The idea is that you run this one program (calling the correct options),  and
you get light curves from images.  This can be done from [a shell
script](https://github.com/waqasbhatti/cdips-pipeline/blob/master/drivers/tess_tuning_drivers/TEMPLATE_single_tess_ccd.sh).
State-awareness (i.e., whether previous reduction attempts succeeded or failed)
is minimal, and based on the pre-existence of files and if-else logic.  A few
similar driver scripts are also in /drivers/, for example to reduce the ETE6
images.  Most of the intermediate files in the pipe-trex reduction (e.g.,
fistar, fiphot, ficonv, etc files) are stored on-disk. A few pieces of metadata
(e.g., image quality diagnostics) are collected to a PostgreSQL database and
used for the reference frame selection steps.

The most important collections of sub-scripts are in `aperturephot.py` and
`imagesubphot.py`, which wrap the `fitsh` tools used to do the aperture
photometry and image subtraction. `autoimagesub.py` rolls up many of these
functions to enable pipeline-like functionality.

## 1. Install

### 1.1 Operating system

This installation is only tested on linux boxes (Debian/Ubuntu).  Compiling the
`fitsh` binaries (see below) might be challenging, but if you do it on MacOS,
please submit a PR describing your experience.


### 1.2 Environment basics

First, clone this repo into a working directory, that we will call $TREX.

  ```
  git clone https://github.com/waqasbhatti/cdips-pipeline
  ```

#### 1.2.1 `python3` environment

Have a working version of [conda](https://conda.io/miniconda.html) on your path
(the 64-bit linux version, presumably).  Then

  ```
  cd $PIPE_TREX
  conda env create -f environment_trex_37.yml -n trex_37
  source activate trex_37
  ```

This deals with most of the details listed below, and is the recommended way to
use pipe-trex. (A variety of "new" features since the 2018 development effort
depend on using a 3.X environment).

An extra step is that to make `pipe-trex` accessible to the virtual
environment, you currently need to add a `.pth` file to the approriate
site-packages directory.  See "Making pipe-trex accessible within virtual
environment" below.


#### 1.2.2 `python2.7` environment

If you opt instead for a 2.7 environment, make a virtual environment in some
local directory, for instance `~/local`. We will call the environment trex\_27:

  ```
  pip install virtualenv; virtualenv --python=/usr/bin/python2.7 trex_27
  ```

Make sure your $PATH and $PYTHON\_PATH environmental variables do not have any
python-related things that will break your venv.

Active the empty python 2.7 environment:

  ```
  source trex_27/bin/activate
  (trex_27) cd $PIPE_TREX
  (trex_27) pip install -r requirements.txt
  ```

This latter step will take some time.  Then, ensure pyeebls and bleeding-edge
astrobase are installed (they are commented out by default because they require
a number of dependencies):

  ```
  (trex_27) pip install pyeebls
  (trex_27) cd $ASTROBASE_DIR
  (trex_27) python setup.py develop
  ```

where I have assumed you may want to make contributions to astrobase as you
develop, which you should!

For transit-fitting, you will want `batman`, `corner`, and `emcee` installed:
  ```
  (trex_27) pip install batman-package
  (trex_27) pip install corner
  (trex_27) cd $SOME_DIR
  (trex_27) git clone https://github.com/dfm/emcee
  (trex_27) cd emcee
  (trex_27) python setup.py install
  ```

### 1.3 catalog reading dependencies

In order to perform photometry, we project known-star catalogs onto images in
order to know where the stars are. (This is more reliable than source
extraction.)

The best catalog in town is Gaia DR2. To access it quickly, `gaia2read` is a
useful program. Jason Kim wrote it for his junior thesis, and his source is
[available here](https://github.com/jkim117/gaia2read).
Sam Yee added an important piece of functionality: cutting on different
magnitudes, his fork is [available
here](https://github.com/samuelyeewl/gaia2read).

To install `gaia2read` to the command line, do the following:
  ```
  (trex_27) cd $SOME_DIRECTORY
  (trex_27) git clone https://github.com/samuelyeewl/gaia2read
  (trex_27) cd gaia2read/gaialib2
  (trex_27) make
  (trex_27) mv gaia2read $BINARY_DIRECTORY_ON_YOUR_PATH
  (trex_27) echo "/nfs/phn15/ar0/H/CAT/GaiaDR2" > ~/.gaia2readrc
  ```
where `$BINARY_DIRECTORY_ON_YOUR_PATH$` is for example `~/bin/`, or some other
directory from which your path can read binaries.

The assumption of the last line is that you are doing this on the
`phn`/`phs`/`phtess` NFS filesystem, where the "DataPreparation" steps have
already been performed to download and sort the Gaia DR2 catalog.


### 1.4 anet and astrometry.net dependencies

You must use either `anet` or `astrometry.net`. The latter is strongly
recommended, since it's free.  To install, follow [this
page](http://astrometry.net/doc/build.html#build). If you're doing wide-field
work, be sure to get both the 4100 and 4200 indexes.


### 1.5 fitsh and HATpipe dependencies

This code inherits from the [`fitsh`](https://fitsh.net/wiki/Main_Page)
project, developed mostly by Andras Pal. Much of `fitsh` was inherited by
`HATpipe` circa 2010, when it forked.  Again, because they are free, we opt for
the public `fitsh` versions, rather than the closed `HATpipe` fork versions.

The utilities we want working on our path include: `ficalib`, `fistar`,
`fiphot`, `grmatch`, `grtrans`, `ficonv`, `fitrans`, and `ficombine`,
`astrometry.net`, and `gaia2read` (deprecated: `2massread`).

Most of these are `fitsh` tasks.  The `fitsh` installation instructions are
[here](https://fitsh.net/wiki/Installation), and they are simple:

  ```
  cd ~/local
  wget http://fitsh.net/download/fitsh/fitsh-0.9.2.tar.gz
  tar xvzf fitsh-0.9.2.tar.gz
  cd fitsh-0.9.2
  ./configure
  make
  make install
  ```
Check to make sure this gives you `ficalib`, `fistar`, `fiphot`, `grmatch`,
`grtrans`, `ficonv`, `fitrans`, and `ficombine`.


### 1.6 PostgreSQL installation

For bookkeeping, you will also need a [PostgreSQL
database](https://www.postgresql.org/files/documentation/pdf/10/postgresql-10-US.pdf).

To install for macs, see [here](https://www.postgresql.org/download/macosx/).
For linux boxes, see
[here](https://wiki.postgresql.org/wiki/Detailed_installation_guides#MacOS).

(If you're installing on Mac OS X, it is a good idea to change your kernel state
by modifying your `/etc/sysctl.conf` file to include things discussed in the
READMEs from the above links.)

Once you've done this:
```
$ psql -U postgres      # login as master user
postgres=# create user hpx with password 'pwgoeshere' createdb;
$ createdb -U hpx hpx
```
and add the appropriate password and info to your `~/.pgpass` file.

To access the database: `psql -U hpx hpx` launches the PostgreSQL database
named `hpx` run by user `hpx`. Or `psql -U hpx -h xphtess1 hpx` does the same,
for a database run on xphtess1, rather than localhost.

To create the tables, run the following:

  ```
  psql$  \i photometry.sql
  psql$  \i xtrnsconvsub.sql
  psql$  \i imagesub-refino.sql
  ```

Beware that these also remove any information you already had in them. The
relations include:

  ```
                       List of relations
   Schema |              Name              |   Type   | Owner
  --------+--------------------------------+----------+-------
   public | ap_photometry                  | table    | hpx
   public | arefshiftedframes              | table    | hpx
   public | arefshiftedframes_framekey_seq | sequence | hpx
   public | astromrefs                     | table    | hpx
   public | calibratedframes               | table    | hpx
   public | calibratedframes_framekey_seq  | sequence | hpx
   public | filters                        | table    | hpx
   public | frameinfo                      | table    | hpx
   public | frameinfo_framekey_seq         | sequence | hpx
   public | iphotfiles                     | table    | hpx
   public | iphotobjects                   | table    | hpx
   public | ism_photometry                 | table    | hpx
   public | lcinfo                         | table    | hpx
   public | objectinfo                     | table    | hpx
   public | photindex_iphots               | table    | hpx
   public | photometryinfo                 | table    | hpx
   public | photrefs                       | table    | hpx
   public | subtractedframes               | table    | hpx
   public | subtractedframes_framekey_seq  | sequence | hpx
  ```

### 1.7 Making pipe-trex accessible within virtual environment

For the moment, go into the venv's `usr/lib/python2.7/site-packages` directory
and create a `.pth` file, e.g.  `pipe-trex.pth` with the location of the local
git cloned repository in it: `/path/to/where/you/cloned/pipe-trex`.

Then activate the virtualenv, and see if you can import a module:

```
py> import imagesubphot as ism
```

For a conda environment, do the same thing, but the site-packages directory
will instead be at a path like
`/home/lbouma/miniconda3/envs/trex_37/lib/python3.7/site-packages`.


## 2 Getting Started

Some usage examples are given in the `drivers/` directory.


### 2.1 Concepts: Directory structure

Everything must be in its right place for your photometry to Just Work. 
During an initial installation, you will need to make a directory structure as
follows, where "." refers to a base directory of your chosing.  (The pipeline
assumes this directory structure has been made).

  ```
  .
  ├── CAL            # You must download calibrated fits images from MAST into here
  ├──                # e.g., ./CAL/sector-4/tess2018307055940-s0004-1-2-0124-s_ffic.fits
  ├── ENGINEERING    # You must download engineering files from MAST into here
  ├──                # e.g., ./ENGINEERING/tess2018344132117_sector04-eng.fits
  ├──                # 
  ├──                # The following are automatically populated:
  ├── BASE           # Reference frames will go here
  ├── RED            # Trimmed and calibrated frames will go here
  ├── LC             # Lightcurves will go here
  ├── PROJ           # Files used to document different reductions will go here
  ├── MOVIES         # Movies will get made and will be put here
  └── REDTEMP        # Temporary space for reduced files
  ```

At the second level:

  ```
  .
  ├── CAL
  │   ├── sector-3          # e.g., ./CAL/sector-4/tess2018307055940-s0004-1-2-0124-s_ffic.fits
  │   └── sector-4          # As above.
  ├── ENGINEERING           # As in the level-1 example.
  ├── BASE
  │   ├── frameinfo-cache   # Once created, these two subdirectories are...
  │   └── reference-frames  # ...automatically populated
  ├── LC
  │   ├── FULL              # For LCs from a full sector
  │   └── TUNE              # For short LCs from a subset of a sector
  ├── BASE                  # These directories have subdirectories that get made automatically.
  ├── RED                   # These directories have subdirectories that get made automatically.
  ├── LC                    # These directories have subdirectories that get made automatically.
  ├── PROJ                  # These directories have subdirectories that get made automatically.
  ├── MOVIES                # These directories have subdirectories that get made automatically.
  └── REDTEMP               # These directories have subdirectories that get made automatically.
  ```

Maintaining this structure is essential. Commands can be run from anywhere,
provided that this structure is maintained.

## 3. Authors

Waqas Bhatti, Luke Bouma, Samuel Yee

## 4. License

MIT

## Appendix. Notes on mac installation

**Aside**: compiling the HATpipe source on a mac is not advised, because many
of the libraries are linux-specific. The entire pipeline is untested on macs,
and the following are some notes from a _failed_ attempt at getting the
pipeline to work on a mac.

  To compile `2massread`:
    ```
    cd /Users/luke/local/HATpipe_R3186/source/2mass
    cp /usr/include/malloc/malloc.h .
    # change 2massread.c's malloc include statement to be `#include "malloc.h"`.
    ./hatconf.sh
    make
    ```
  you then need to put the appropriately formatted ~150Gb of 2MASS index files
  somewhere accessible, and point to them in your `~/.2massreadrc` file.

  To compile `anrd2xy`:
    ```
    cd /Users/luke/local/HATpipe_R3186/source/odoncontrib/anet
    brew install gsl  # this is not a default on macs
    ```

  If on a mac, you then must edit all six makefiles,
  ```
    /Users/luke/local/HATpipe_R3186/source/odoncontrib/anet/Makefile
    /Users/luke/local/HATpipe_R3186/source/odoncontrib/anet/libc*/Makefile
    /Users/luke/local/HATpipe_R3186/source/odoncontrib/anet/libc*/*/Makefile
  ```
  to use GNU `gcp`, not `cp`, because mac `cp` has different options.  Even then
  though, linking on my mac fails because of architecture problems that I don't
  understand. This is perhaps a waste of time, and you should just develop on
  linux, if you have a very good internet connection, or do not develop.

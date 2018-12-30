# `pipe-trex`

Timeseries photometry pipeline used for HATPI prototype, TESS simulated
data, and real TESS data.

## Install

### Environment basics

First, clone this repo into a working directory, that we will call $TREX.

  ```
  git clone https://github.com/waqasbhatti/pipe-trex
  ```

Make a virtual environment in some local directory, for instance `~/local`. We
will call the environment trex\_27:

  ```
  pip install virtualenv; virtualenv --python=/usr/bin/python2.7 trex_27
  ```

Make sure your $PATH and $PYTHON\_PATH environmental variables do not have any
python-related things that will break your venv.

Active the empty python 2.7 environment: 

  ```
  source trex_27/bin/activate
  ```

Change directories to $TREX, and install everything pip can find:

  ```
  (trex_27) pip install -r requirements.txt
  ```

This will take some time.  Then, ensure pyeebls and bleeding-edge astrobase are
installed (they are commented out by default because they require a number of
dependencies):

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
  (trex_27) pip corner
  (trex_27) cd $SOME_DIR
  (trex_27) git clone https://github.com/dfm/emcee
  (trex_27) cd emcee
  (trex_27) python setup.py install
  ```

## catalog reading dependencies

At various stages in the photometry process, we project known-star catalogs
onto images in order to know where the stars are. (This is more reliable than
source extraction.)

The best catalog in town is Gaia DR2. To access it quickly, `gaia2read` is a
useful program. Jason Kim wrote it for his junior thesis, and his source is
[available here](https://github.com/jkim117/gaia2read).
Sam Yee added an important piece of functionality: cutting on different
magnitudes, his fork is [available
here](https://github.com/samuelyeewl/gaia2read).

To install `gaia2read`, you will want to do the following:
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


### anet and astrometry.net dependencies

You must use either `anet` or `astrometry.net`. I recommend the later, since
it's free.  To install, follow [this
page](http://astrometry.net/doc/build.html#build). If you're doing wide-field
work, be sure to get both the 4100 and 4200 indexes.


### fitsh and HATpipe dependencies

This code inherits a large amount from the
[`fitsh`](https://fitsh.net/wiki/Main_Page) project, developed mostly by Andras
Pal. Much of `fitsh` was inherited by `HATpipe` circa 2010, when it forked.
`HATpipe` now has top-secret special-sauce added, that you are not allowed to
see or reuse without express written permission.

The utilities we want working on our path include: `ficalib`, `fistar`,
`fiphot`, `grmatch`, `grtrans`, `ficonv`, `fitrans`, and `ficombine`, `anet`,
and `2massread`.

If you had access to the top-secret HATpipe svn repo, you would do something
like this:
  ```
  cd ~/local
  svn checkout --username lbouma https://hat.astro.princeton.edu/svn/HATreduc/HATpipe
  cd HATpipe
  svn up -r3186   # we froze dependencies for pipe-trex
  cd ../ ; cp -r HATpipe HATpipe_R3186
  ```

Compiling the HATpipe source on a mac is challenging, because many of the
libraries are linux-specific.


Most of these are `fitsh` tasks. We will proceed using only free software
unless absolutely necessary, which means using Andras' public version.
The `fitsh` installation instructions are
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

### PostgreSQL installation

For macs, see [here](https://www.postgresql.org/download/macosx/).
For linux boxes, see
[here](https://wiki.postgresql.org/wiki/Detailed_installation_guides#MacOS).

NB if you're installing on Mac OS X, it is a good idea to change your kernel state
by modifying your `/etc/sysctl.conf` file to include things discussed in the
READMEs from the above links.

Once you've done this:
```
$ psql -U postgres      # login as master user
postgres=# create user hpx with password 'pwgoeshere' createdb;
$ createdb -U hpx hpx
```
and add the appropriate password and info to your `~/.pgpass` file.
  

### Making pipe-trex accessible within virtual environment
For the moment, go into the venv's `usr/lib/python2.7/site-packages` directory
and create a `.pth` file, e.g.  `pipe-trex.pth` with the location of the local
git cloned repository in it: `/path/to/where/you/cloned/pipe-trex`.

Then activate the virtualenv, and see if you can import a module:

```
py> import imagesubphot as ism
```


## Getting Started

Some usage examples are given in the `examples/` directory.


### Concepts: Directory and database structure

As Thom Yorke understood, everything must be in its right place for your
photometry to Just Work. The assumed directory structure at its first level
must look like this:

  ```
  .
  ├── BASE      # 
  ├── LC        # lightcurves go here
  ├── PROJ      # project directories for specific users
  ├── RAW       # raw fits images (or fits.fz images) go here
  ├── RED       # calibrated fits images go here
  └── REDTEMP   # temporary space for reduced files
  ```

At the second level:

  ```
  .
  ├── BASE
  │   ├── CAL                 # calibration frames (bias, flat, dark) go here
  │   ├── CAT                 # FOV source catalogs for your field go here
  │   └── reference-frames    # astrometric reference frames; sqlite3 databases
  ├── LC
  │   ├── projid12-G577       # lightcurves specific to particular projectids
  │   ├── stats_files         # for lightcurve statistics
  │   └── projid8-9-10-11-12
  ├── PROJ
  │   ├── lbouma
  │   └── wbhatti
  ├── RAW
  │   ├── 1-20161109          # ??? each night gets its own directory under RAW
  │   ├── 1-20161110
  │   ├── 1-20161111
  │   ├── 1-20161130
  │   └── tmp
  ├── RED
  │   ├── projid12            # each project gets its own subdirectory under RED
  │   └── wbhatti-red-frames
  ```

Maintaining this structure is essential. Commands can be run from anywhere,
provided that this structure is maintained.

For bookkeeping, you will also need a [PostgreSQL
database](https://www.postgresql.org/files/documentation/pdf/10/postgresql-10-US.pdf).

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

### Usage

* `shared_variables.py`: contains path variables that are set across the
  pipeline.


# Useful things

Most of the tasks above produce JPEGs that can be examined to see how
everything is going. Use the following command (if you have ffmpeg installed
and libx264 installed) to make a movie of these frames:

`ffmpeg -framerate 60 -pattern_type glob -i '*.jpg' -c:v libx264 -preset fast out.mp4`

replace the \*.jpg above with the appropriate glob pattern to use for the files
in question

# Authors

Waqas Bhatti

_Contributors:_
Luke Bouma

# License

MIT

# Notes on mac installation

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


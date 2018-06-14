# `pipe-trex`

Timeseries photometry pipeline used for HATPI prototype, and TESS simulated
data.

# Install

TODO FOR INSTALL INSTRUCTIONS:
* figure out `anet` substitute (can it be replaced by just using
  astrometry.net?)
* check this actually works.


__Environment basics__

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

This will take some time.  Then, ensure pyeebls and astrobase are installed
(they are commented out by default because they require a number of
dependencies):

  ```
  (trex_27) pip install pyeebls
  (trex_27) pip install astrobase
  ```

__fitsh (and HATpipe) dependencies__

This code inherits a large amount from the
[`fitsh`](https://fitsh.net/wiki/Main_Page) project, developed mostly by Andras
Pal. Much of `fitsh` was inherited by `HATpipe` circa 2010, when it forked.
`HATpipe` now has top-secret special-sauce added, that you are not allowed to
reuse without express written permission.

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

------------------------------

The utilities we want working on our path include: `anet`, `ficalib`, `fistar`,
`fiphot`, `grmatch`, `grtrans`, `ficonv`, `fitrans`, and `ficombine`.

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

__PostgreSQL installation__

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
  

__Making pipe-trex accessible within virtual environment__
For the moment, go into the venv's `usr/lib/python2.7/site-packages` directory
and create a `.pth` file, e.g.  `pipe-trex.pth` with the location of the local
git cloned repository in it: `/path/to/where/you/cloned/pipe-trex`.

Then activate the virtualenv, and see if you can import a module:

```
py> import imagesubphot as ism
```


# Getting Started

Congratulations on installing things.

Some usage examples are given in the `examples/` directory.


## Concepts: Directory and database structure

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

## Usage

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

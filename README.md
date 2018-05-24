# `pipe-trex`
Timeseries photometry pipeline for the HATPI prototype

**this readme is currently under construction**

# Prerequisites

# Installation

# Getting Started

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
│   └── reference-frames    # astrometric reference frames, and their databases
├── LC
│   ├── projid12-G577       # lightcurves specific to particular projectids
│   ├── projid16-18
│   └── projid8-9-10-11-12
├── PROJ
│   ├── lbouma
│   └── wbhatti
├── RAW
│   ├── 1-20161109          # ??? each night gets its own directory under RAW
│   ├── 1-20161110
│   ├── 1-20161111
│   ├── 1-20161112
│   ├── 1-20161129
│   ├── 1-20161130
│   └── tmp
├── RED
│   ├── projid12            # each project gets its own subdirectory under RED
│   └── wbhatti-red-frames
```

Maintaining this structure is essential. Commands can be run from anywhere,
provided that this structure is maintained.

* `shared_variables.py`: contains path variables that are set across the
  pipeline.

# Authors

# License

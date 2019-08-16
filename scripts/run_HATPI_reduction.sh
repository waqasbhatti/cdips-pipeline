#!/usr/bin/env bash

# Run aperture tuning experiments for HATPI
BASEDIR="/nfs/phtess1/ar1/HATPI/HP0/RED/projid12-G577-ccd8-sdssr/"
FITSGLOB="?-???????_?"
PROJDIR="/home/syee/work/reduction-G577-proj1299/"

fitsfiles=("$BASEDIR$FITSGLOB.fits")
for f in $fitsfiles; do
    ln $f $PROJDIR
done


# Experiment parameters
aperturelist="1.95:6.0:5.0,2.45:6.0:5.0,2.95:6.0:5.0"
projectid=1207
subsample="tuning"
catalog="GAIADR2-RA277.5-DEC-22.5-SIZE14.0"

# Create new directory and link files that don't have to be changed
basesubdir="$BASEDIR$subsample-base/"
projdir="$BASEDIR$subsample-$projectid/"

if [ ! -d "$projdir" ]; then
    mkdir "$projdir"
fi

echo "Linking .fits, .fistar, .wcs, .xysdk, .itrans, .-xtrns.fits files to $projdir"
fitsfiles=("$basesubdir$FITSGLOB.fits")
for f in $fitsfiles; do
    bn=$(basename -- "$f" .fits)
    if [ ! -e "$projdir$bn.fits" ]; then
        ln $f $projdir
    fi
    if [ ! -e "$projdir$bn.fistar" ]; then
        ln $basesubdir$bn.fistar $projdir
    fi
    if [ ! -e "$projdir$bn.wcs" ]; then
        ln $basesubdir$bn.wcs $projdir
    fi
    if [ ! -e "$projdir$bn.xysdk" ]; then
        ln $basesubdir$bn.xysdk $projdir
    fi
    if [ ! -e "$projdir$bn.itrans" ]; then
        ln $basesubdir$bn.itrans $projdir
    fi
    if [ ! -e "$projdir$bn-xtrns.fits" ]; then
        ln $basesubdir$bn-xtrns.fits $projdir
    fi
    if [ ! -e "$projdir$bn.projcatalog" ]; then
        ln $basesubdir$bn.projcatalog $projdir
    fi
done

echo "Linking catalog to $projdir"
catfile="$basesubdir$catalog.catalog"
refcatfile="$basesubdir$catalog.reformed-catalog"
if [ ! -e "$projdir$catalog.catalog" ]; then
    ln $catfile $projdir
fi
if [ ! -e "$projdir$catalog.reformed-catalog" ]; then
    ln $refcatfile $projdir
fi

echo "Running aperture experiment $aperturelist, projid $projectid"
echo "Experiment $projid parameters:" > "$projdir/$projectid-exp.txt"
echo "aperturelist=$aperturelist" >> "$projdir/$projectid-exp.txt"

python HATPI_reduction.py tune -d $projdir --profile --projid $projectid --nworkers 28 --photrelaxcriteria --aperturelist $aperturelist | tee "$subsample-$projectid.log"

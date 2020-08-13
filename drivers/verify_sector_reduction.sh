#
# Symlink the plots one has to manually look at in order to verify the output
# of a reduction to a single directory, to enable faster geeqie click-throughs.
#

##########
# config #
##########
snum=0001

##########
# script #
##########
echo 'VERIFYING s'$snum

# check number of light-curves is reasonable
cat `tail -n1 -q logs/s$snum-cam*log | awk '{print $2}'`

# symlink stats_files directories to the "verification" directory
verificationdir='/nfs/phn12/ar0/H/PROJ/lbouma/cdips-pipeline/drivers/verification/s'$snum

if [ ! -d $verificationdir ]; then
  mkdir $verificationdir
fi

for cam in {1..4}; do
  for ccd in {1..4}; do
    src=/nfs/phtess2/ar0/TESS/FFI/LC/FULL/s$snum/ISP_$cam-$ccd-????/stats_files
    projid=$(basename $(dirname $src) | awk -F "-" '{print $3}')
    dst=$verificationdir/cam${cam}_ccd${ccd}_projid${projid}_stats_files
    ln -s $src $dst
    echo 'symlinked' $src 'to' $dst
    rm $dst/*stats_files
  done
done


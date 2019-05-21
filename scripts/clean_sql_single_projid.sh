projectid=1500

for projectid in {1500..1515}; do

  psql -U hpx -h xphtess1 hpx -c \
    "DELETE from calibratedframes where (fitsheader->'PROJID' ='"${projectid}"');"

  psql -U hpx -h xphtess1 hpx -c \
    "delete from astromrefs where projectid = "${projectid}";"

  psql -U hpx -h xphtess1 hpx -c \
    "delete from photrefs where projectid = "${projectid}";"

done

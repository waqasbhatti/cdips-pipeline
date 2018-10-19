#!/usr/bin/env bash

# PURPOSE: USE WITH CAUTION.
#
# wipe and recreate the following postgresql tables:
#
# xtrnsconvsub.sql: 
#   arefshiftedframes
#   subtractedframes
# ismphot-database.sql:
#   iphotfiles
#   iphotobjects
#   astromrefs
#   photrefs
# photometry.sql:
#   objectinfo
#   lcinfo
#   filters
#   ism_photometry
#   ap_photometry
#   calibrateframes
#   photometryinfo
#
# One would do such a thing if you didn't believe values you had input into
# them.


psql -U hpx -h xphtess1 hpx -c "\i ../xtrnsconvsub.sql;"
psql -U hpx -h xphtess1 hpx -c "\i ../ismphot-database.sql;"
psql -U hpx -h xphtess1 hpx -c "\i ../photometry.sql;"

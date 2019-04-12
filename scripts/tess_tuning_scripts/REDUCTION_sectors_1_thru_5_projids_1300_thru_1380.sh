#!/usr/bin/env bash

####################################
# USAGE: ./REDUCTION_sectors_1_thru_5_projids_1300_thru_1380.sh & #
####################################

# sector 2, all cam/ccds. started 2019/03/03

# ( source activate trex_37; source projid_1300.sh; wait ) & wait
# ( source activate trex_37; source projid_1301.sh; wait ) & wait
# ( source activate trex_37; source projid_1302.sh; wait ) & wait
# ( source activate trex_37; source projid_1303.sh; wait ) & wait
# ( source activate trex_37; source projid_1304.sh; wait ) & wait
# ( source activate trex_37; source projid_1305.sh; wait ) & wait
# ( source activate trex_37; source projid_1306.sh; wait ) & wait

# completed 2019/03/06 (& in between)

# # stop to get cam 4 early. 2019/03/06
# ( source activate trex_37; source projid_1312.sh; wait ) & wait
# ( source activate trex_37; source projid_1313.sh; wait ) & wait
# ( source activate trex_37; source projid_1314.sh; wait ) & wait
# ( source activate trex_37; source projid_1315.sh; wait ) & wait
# 
# # restart 1307, which was interrupted!
# ( source activate trex_37; source projid_1307.sh; wait ) & wait
# ( source activate trex_37; source projid_1308.sh; wait ) & wait
# ( source activate trex_37; source projid_1309.sh; wait ) & wait

# 1307 and 1308 completed.  1309 was interrupted by Sam requesting a stop.
# 20190308 restarting 1309. Note 2019/03/06 cam 4 runs were buggy, perhaps
# because of the LMC subtraction issues.
# ( source activate trex_37; source projid_1309.sh; wait ) & wait

# 20190309: above seems to have mysteriously failed. fix it later. 1378
# contains the cluster you want...
# ( source activate trex_37; source projid_1378.sh; wait ) & wait

# 20190313: 1378 completed. Try these (sector 2, cam3.):
# is photometric performance changed b/c of the photref change?
# proj1310 is the first unbiased check of this.
# ( source activate trex_37; source projid_1309.sh; wait ) & wait
# ( source activate trex_37; source projid_1310.sh; wait ) & wait

# 20190315: above basically worked. Run the following. (MISSING 1312-1315 b/c
# they're camera4 fields). Start 190315.

# ( source activate trex_37; source projid_1311.sh; wait ) & wait

# 20190319 above finished. below for some reason didnt start. omg.

# 20190320 sector 3, cameras 1 through 3.
# ( source activate trex_37; source projid_1316.sh; wait ) & wait
# ( source activate trex_37; source projid_1317.sh; wait ) & wait
# ( source activate trex_37; source projid_1318.sh; wait ) & wait
# ( source activate trex_37; source projid_1319.sh; wait ) & wait
# ( source activate trex_37; source projid_1320.sh; wait ) & wait
# ( source activate trex_37; source projid_1321.sh; wait ) & wait
# ( source activate trex_37; source projid_1322.sh; wait ) & wait
# ( source activate trex_37; source projid_1323.sh; wait ) & wait
# ( source activate trex_37; source projid_1324.sh; wait ) & wait
# ( source activate trex_37; source projid_1325.sh; wait ) & wait
# ( source activate trex_37; source projid_1326.sh; wait ) & wait
# ( source activate trex_37; source projid_1327.sh; wait ) & wait

# 20190320. above died at making the ccd temperature pickles, b/c i hadn't
# downloaded the engineering files to drive. but ran
# 20190320 sector 3, cameras 1 through 3.
( source activate trex_37; source projid_1316.sh; wait ) & wait
( source activate trex_37; source projid_1317.sh; wait ) & wait
( source activate trex_37; source projid_1318.sh; wait ) & wait
( source activate trex_37; source projid_1319.sh; wait ) & wait
( source activate trex_37; source projid_1320.sh; wait ) & wait
( source activate trex_37; source projid_1321.sh; wait ) & wait
( source activate trex_37; source projid_1322.sh; wait ) & wait
( source activate trex_37; source projid_1323.sh; wait ) & wait
( source activate trex_37; source projid_1324.sh; wait ) & wait
( source activate trex_37; source projid_1325.sh; wait ) & wait
( source activate trex_37; source projid_1326.sh; wait ) & wait
( source activate trex_37; source projid_1327.sh; wait ) & wait

# # sector 4, cameras 1 through 3. NOTE: before running automatically, might
# need to fine-tune the orbitgap!!
# ( source activate trex_37; source projid_1332.sh; wait ) & wait
# ( source activate trex_37; source projid_1333.sh; wait ) & wait
# ( source activate trex_37; source projid_1334.sh; wait ) & wait
# ( source activate trex_37; source projid_1335.sh; wait ) & wait
# ( source activate trex_37; source projid_1336.sh; wait ) & wait
# ( source activate trex_37; source projid_1337.sh; wait ) & wait
# ( source activate trex_37; source projid_1338.sh; wait ) & wait
# ( source activate trex_37; source projid_1339.sh; wait ) & wait
# ( source activate trex_37; source projid_1340.sh; wait ) & wait
# ( source activate trex_37; source projid_1341.sh; wait ) & wait
# ( source activate trex_37; source projid_1342.sh; wait ) & wait
# ( source activate trex_37; source projid_1343.sh; wait ) & wait
# 
# # sector 5, cameras 1 through 3
# ( source activate trex_37; source projid_1348.sh; wait ) & wait
# ( source activate trex_37; source projid_1349.sh; wait ) & wait
# ( source activate trex_37; source projid_1350.sh; wait ) & wait
# ( source activate trex_37; source projid_1351.sh; wait ) & wait
# ( source activate trex_37; source projid_1352.sh; wait ) & wait
# ( source activate trex_37; source projid_1353.sh; wait ) & wait
# ( source activate trex_37; source projid_1354.sh; wait ) & wait
# ( source activate trex_37; source projid_1355.sh; wait ) & wait
# ( source activate trex_37; source projid_1356.sh; wait ) & wait
# ( source activate trex_37; source projid_1357.sh; wait ) & wait
# ( source activate trex_37; source projid_1358.sh; wait ) & wait
# ( source activate trex_37; source projid_1359.sh; wait ) & wait
# 
# # sector 1, cameras 1 through 3
# ( source activate trex_37; source projid_1364.sh; wait ) & wait
# ( source activate trex_37; source projid_1365.sh; wait ) & wait
# ( source activate trex_37; source projid_1366.sh; wait ) & wait
# ( source activate trex_37; source projid_1367.sh; wait ) & wait
# # ( source activate trex_37; source projid_1368.sh; wait ) & wait
# ( source activate trex_37; source projid_1369.sh; wait ) & wait
# ( source activate trex_37; source projid_1370.sh; wait ) & wait
# ( source activate trex_37; source projid_1371.sh; wait ) & wait
# ( source activate trex_37; source projid_1372.sh; wait ) & wait
# # ( source activate trex_37; source projid_1373.sh; wait ) & wait
# ( source activate trex_37; source projid_1374.sh; wait ) & wait
# ( source activate trex_37; source projid_1375.sh; wait ) & wait

# sectors 2-5, cam 4. (a) 
# ( source activate trex_37; source projid_1312.sh; wait ) & wait
# ( source activate trex_37; source projid_1313.sh; wait ) & wait
# ( source activate trex_37; source projid_1314.sh; wait ) & wait
# ( source activate trex_37; source projid_1315.sh; wait ) & wait

# ( source activate trex_37; source projid_1328.sh; wait ) & wait
# ( source activate trex_37; source projid_1329.sh; wait ) & wait
# ( source activate trex_37; source projid_1330.sh; wait ) & wait
# ( source activate trex_37; source projid_1331.sh; wait ) & wait

# ( source activate trex_37; source projid_1344.sh; wait ) & wait
# ( source activate trex_37; source projid_1345.sh; wait ) & wait
# ( source activate trex_37; source projid_1346.sh; wait ) & wait
# ( source activate trex_37; source projid_1347.sh; wait ) & wait

# ( source activate trex_37; source projid_1360.sh; wait ) & wait
# ( source activate trex_37; source projid_1361.sh; wait ) & wait
# ( source activate trex_37; source projid_1362.sh; wait ) & wait
# ( source activate trex_37; source projid_1363.sh; wait ) & wait

# ( source activate trex_37; source projid_1376.sh; wait ) & wait
# ( source activate trex_37; source projid_1377.sh; wait ) & wait
### ( source activate trex_37; source projid_1378.sh; wait ) & wait
# ( source activate trex_37; source projid_1379.sh; wait ) & wait

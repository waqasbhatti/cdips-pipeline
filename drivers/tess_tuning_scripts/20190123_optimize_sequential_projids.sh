#!/usr/bin/env bash

##########################################
# MUST BE CALLED TO BACKGROUND!! ./sequential_projids.sh & #
##########################################

# kernel experiments:
# * default: "i/5;b/5;d=3/2"
# * just do identity term
# * change background order
# * change delta kernel order
# * change delta kernel half-size
# * change "identity" term order

# aperture experiments:
# * default: "0.71:7.0:6.0,1.41:7.0:6.0,2.82:7.0:6.0"
# * change the inner and outer boundary of background annulus
# * vary the aperture sizes

( source activate trex_37; source projid_1033.sh; wait ) & wait
( source activate trex_37; source projid_1034.sh; wait ) & wait
( source activate trex_37; source projid_1035.sh; wait ) & wait
( source activate trex_37; source projid_1036.sh; wait ) & wait
( source activate trex_37; source projid_1037.sh; wait ) & wait
( source activate trex_37; source projid_1038.sh; wait ) & wait
( source activate trex_37; source projid_1039.sh; wait ) & wait
( source activate trex_37; source projid_1040.sh; wait ) & wait
( source activate trex_37; source projid_1041.sh; wait ) & wait
( source activate trex_37; source projid_1042.sh; wait ) & wait
( source activate trex_37; source projid_1043.sh; wait ) & wait
( source activate trex_37; source projid_1044.sh; wait ) & wait
( source activate trex_37; source projid_1045.sh; wait ) & wait
( source activate trex_37; source projid_1046.sh; wait ) & wait
( source activate trex_37; source projid_1047.sh; wait ) & wait
( source activate trex_37; source projid_1048.sh; wait ) & wait
( source activate trex_37; source projid_1049.sh; wait ) & wait
( source activate trex_37; source projid_1050.sh; wait ) & wait
( source activate trex_37; source projid_1051.sh; wait ) & wait
( source activate trex_37; source projid_1052.sh; wait ) & wait
( source activate trex_37; source projid_1053.sh; wait ) & wait
( source activate trex_37; source projid_1054.sh; wait ) & wait
( source activate trex_37; source projid_1055.sh; wait ) & wait
( source activate trex_37; source projid_1056.sh; wait ) & wait
( source activate trex_37; source projid_1057.sh; wait ) & wait
( source activate trex_37; source projid_1058.sh; wait ) & wait
( source activate trex_37; source projid_1059.sh; wait ) & wait
( source activate trex_37; source projid_1060.sh; wait ) & wait
( source activate trex_37; source projid_1061.sh; wait ) & wait
( source activate trex_37; source projid_1062.sh; wait ) & wait
( source activate trex_37; source projid_1063.sh; wait ) & wait
( source activate trex_37; source projid_1064.sh; wait ) & wait
( source activate trex_37; source projid_1065.sh; wait ) & wait
( source activate trex_37; source projid_1066.sh; wait ) & wait
( source activate trex_37; source projid_1067.sh; wait ) & wait
( source activate trex_37; source projid_1068.sh; wait ) & wait
( source activate trex_37; source projid_1069.sh; wait ) & wait
( source activate trex_37; source projid_1070.sh; wait ) & wait
( source activate trex_37; source projid_1071.sh; wait ) & wait
( source activate trex_37; source projid_1072.sh; wait ) & wait
( source activate trex_37; source projid_1073.sh; wait ) & wait
( source activate trex_37; source projid_1074.sh; wait ) & wait
( source activate trex_37; source projid_1075.sh; wait ) & wait
( source activate trex_37; source projid_1076.sh; wait ) & wait
( source activate trex_37; source projid_1077.sh; wait ) & wait
( source activate trex_37; source projid_1078.sh; wait ) & wait
( source activate trex_37; source projid_1079.sh; wait ) & wait
( source activate trex_37; source projid_1080.sh; wait ) & wait
( source activate trex_37; source projid_1081.sh; wait ) & wait
( source activate trex_37; source projid_1082.sh; wait ) & wait
( source activate trex_37; source projid_1083.sh; wait ) & wait
( source activate trex_37; source projid_1084.sh; wait ) & wait
( source activate trex_37; source projid_1085.sh; wait ) & wait
( source activate trex_37; source projid_1086.sh; wait ) & wait
echo 'submitted optimization jobs' > 20190123_run_status.txt

#!/usr/bin/env bash

####################################
# USAGE: ./sector13_reduction.sh & #
####################################

( source activate trex_37; source projid_1612.sh; wait ) & wait
( source activate trex_37; source projid_1613.sh; wait ) & wait
( source activate trex_37; source projid_1614.sh; wait ) & wait
( source activate trex_37; source projid_1616.sh; wait ) & wait
( source activate trex_37; source projid_1618.sh; wait ) & wait
( source activate trex_37; source projid_1619.sh; wait ) & wait
( source activate trex_37; source projid_1620.sh; wait ) & wait
( source activate trex_37; source projid_1621.sh; wait ) & wait
( source activate trex_37; source projid_1622.sh; wait ) & wait
( source activate trex_37; source projid_1623.sh; wait ) & wait
( source activate trex_37; source projid_1624.sh; wait ) & wait
( source activate trex_37; source projid_1625.sh; wait ) & wait
( source activate trex_37; source projid_1626.sh; wait ) & wait
( source activate trex_37; source projid_1627.sh; wait ) & wait

# NOTE: requires rest of cluster to be down
( source activate trex_37; source projid_1617.sh; wait ) & wait
( source activate trex_37; source projid_1615.sh; wait ) & wait

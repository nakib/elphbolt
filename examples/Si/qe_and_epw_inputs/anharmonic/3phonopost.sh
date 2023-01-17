#!/bin/bash

WDIR=anharm_sc555_6nn
cd $WDIR

find DISP_3RD* -name scf.out | sort -n | ~/thirdorder/thirdorder_espresso.py si.in reap 5 5 5 -6

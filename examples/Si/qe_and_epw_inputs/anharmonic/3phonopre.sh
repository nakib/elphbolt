#!/bin/bash

WDIR=anharm_sc555_6nn
cd $WDIR

~/thirdorder/thirdorder_espresso.py si.in sow 5 5 5 -6 si_sc.in

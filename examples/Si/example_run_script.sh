#!/bin/bash

workdir="./Si_10r3_300K_CBM/"
inputdir="./input"

mkdir $workdir
cd $workdir

#copy input file
cp ../$inputdir/input.nml .

#make soft link to the rest of the input data
ln -s ../$inputdir/rcells_g .
ln -s ../$inputdir/rcells_k .
ln -s ../$inputdir/rcells_q .
ln -s ../$inputdir/wsdeg_g .
ln -s ../$inputdir/wsdeg_k .
ln -s ../$inputdir/wsdeg_q .
ln -s ../$inputdir/epwdata.fmt .
ln -s ../$inputdir/epmatwp1 .
ln -s ../$inputdir/FORCE_CONSTANTS_3RD .
ln -s ../$inputdir/espresso.ifc2 .

#gcc+opencoarrays
#Call elphbolt, for example, like this to run with 6 coarray images:
cafrun -n 8 elphbolt
##

cd ..

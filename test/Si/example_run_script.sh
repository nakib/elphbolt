#!/bin/bash

workdir="<enter suitable working directory>"
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
#ln -s ../$inputdir/highsympath.txt .
#ln -s ../$inputdir/initialk.txt .

#Call elphbolt.x, for example, like this:
#~/OpenCoarrays/opencoarrays-install/bin/cafrun -n 2 ~/elphbolt/elphbolt.x

#You need to adapt this script for submission to HPC queueing systems.

cd ..

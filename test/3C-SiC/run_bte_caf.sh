#!/bin/bash

# Download data from public repo
curl -X GET "https://nomad-lab.eu/prod/v1/api/v1/uploads/5b-c4JYRSG6qitMFjFUUzQ/raw/?offset=0&length=-1&compress=true" -o 3C-SiC.zip

#Extract data
unzip 3C-SiC.zip
rm 3C-SiC.zip

workdir="./3C-SiC_test_output"
datadir="./3C-SiC/elphbolt_input_data"
inputdir="../test/3C-SiC"

mkdir $workdir
cp $datadir/* $workdir #symbolic links don't seem to work with CMake :(
cp $inputdir/input.nml $workdir
rm -r 3C-SiC
cd $workdir

#gcc+opencoarrays
# if test fails, propagate the error
cafrun -n 2 ../bin/test_bte_regression || exit -1

cd ..

#!/bin/bash

# Download data from public repo
#curl -X GET "https://nomad-lab.eu/prod/v1/api/v1/uploads/5b-c4JYRSG6qitMFjFUUzQ/raw/?offset=0&length=-1&compress=true" -o 3C-SiC.zip

#Extract data
#unzip 3C-SiC.zip
#rm 3C-SiC.zip

workdir="./3C-SiC_test_output"
datadir="./3C-SiC/elphbolt_input_data"
inputdir="./test/3C-SiC"

#mkdir -p $workdir
#cp $datadir/* $workdir
#cp $inputdir/input.nml $workdir
#rm -r 3C-SiC
cd $workdir

#gcc+opencoarrays
cafrun -n 1 ../build/caf_*/test/check_interactions_symmetries | tee 3C_SiC_test.output

cd ..

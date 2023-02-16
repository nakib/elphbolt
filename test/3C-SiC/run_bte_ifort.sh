#!/bin/bash

workdir="./3C-SiC_test_output"
inputdir="../test/3C-SiC/input"

mkdir $workdir
cp $inputdir/* $workdir #symbolic links won't work with CMake
cd $workdir

#intel
export FOR_COARRAY_NUM_IMAGES=4
# if test fails, propagate the error
~/elphbolt/build/bin/test_bte_regression || exit -1

cd ..

#rm -r $workdir

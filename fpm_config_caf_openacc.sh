#!/bin/bash

#By default, try to use the shipped forlapack and forblas libraries.
#These will have been placed in the following directory when install_blas+lapack.sh was executed.
FORBLASLAPACK_DIR="~/.local/lib/"

export FPM_FC=caf
#bounds check does not work with -fopenacc...not sure why...
#export FPM_FFLAGS="-g -O2 -cpp -fopenacc -DOPENACC=1 -foffload=-lm -fopt-info-omp -fPIC -Wunused -Wconversion -Wunderflow -Wdo-subscript -L $FORBLASLAPACK_DIR"
export FPM_FFLAGS="-g -O2 -cpp -no-pie -fopenacc -DOPENACC=1 -foffload=nvptx-none -foffload=-lm -fopt-info-omp -fPIC -Wunused -Wconversion -Wunderflow -Wdo-subscript -L $FORBLASLAPACK_DIR"

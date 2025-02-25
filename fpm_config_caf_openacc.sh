#!/bin/bash

#By default, try to use the shipped forlapack and forblas libraries.
#These will have been placed in the following directory when install_blas+lapack.sh was executed.
FORBLASLAPACK_DIR="~/.local/lib/"

export FPM_FC=caf
#bounds check does not work with -fopenacc...not sure why...
export FPM_FFLAGS="-g -O3 -cpp -fopenacc -fopt-info-omp -DOPENACC=1 -foffload=-lm -foffload=nvptx-none -no-pie -fPIC -Wunused -Wconversion -Wunderflow -Wdo-subscript -L $FORBLASLAPACK_DIR"

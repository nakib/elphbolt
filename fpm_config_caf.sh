#!/bin/bash

#By default, try to use the shipped forlapack and forblas libraries.
#These will have been placed in the following directory when install_blas+lapack.sh was executed.
FORBLASLAPACK_DIR="~/.local/lib/"

export FPM_FC=caf
export FPM_FFLAGS="-g -O2 -fPIC -fbounds-check -Wunused -Wconversion -Wunderflow -Wdo-subscript -L $FORBLASLAPACK_DIR"

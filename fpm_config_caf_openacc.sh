#!/bin/bash

export FPM_FC=caf
#bounds check does not work with -fopenacc...not sure why...
#-fno-fast-math needed to use fortran intrinsic matmul
export FPM_FFLAGS="-g -O2 -fopenacc -foffload=-lm -foffload=-fno-fast-math -fopt-info-omp -fPIC -Wunused -Wconversion -Wunderflow -Wdo-subscript"

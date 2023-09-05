#!/bin/bash

export FPM_FC=caf
#bounds check does not work with -fopenacc...not sure why...
export FPM_FFLAGS="-g -O2 -cpp -fopenacc -DOPENACC=1 -foffload=-lm -fopt-info-omp -fPIC -Wunused -Wconversion -Wunderflow -Wdo-subscript"

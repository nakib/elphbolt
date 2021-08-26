export FFLAGS=-g -O2 -fPIC -fbounds-check -Wunused -Wconversion -Wunderflow -Wdo-subscript 
export LDFLAGS=-L/usr/public/lapack/3.4.2/lib -L/usr/public/openblas/OpenBLAS
export FC=~/OpenCoarrays/opencoarrays-install/bin/caf
export LIBS=-llapack -lblas -lsymspg

#libsymspg-dev and libsymspg1 version 1.14.1 on Ubuntu 20.04.2 LTS

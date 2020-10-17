#export FFLAGS=-g -O2 -fPIC -fbounds-check
export LDFLAGS=-L/usr/public/lapack/3.4.2/lib -L/usr/public/openblas/OpenBLAS -L/usr/public/spglib/1.6.0/lib -lsymspg
export FC=~/OpenCoarrays/opencoarrays-install/bin/caf
#export LIBS=-lgfortran -llapack -lblas -lsymspg
export LIBS=-llapack -lblas -lsymspg

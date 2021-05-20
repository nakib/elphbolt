export FFLAGS=-g -O2 -fbounds-check
export LDFLAGS=-L/home/icn85/icn85265/spglib-1.6.0/build/lib -lsymspg
export FC=caf
export LIBS=-lblas -llapack

#module load gcc/7.2.0 openmpi/3.1.1 opencoarrays/2.6.3 lapack/3.8.0 openblas/0.3.6
#spglib-1.6.0 built with gcc and path to ~/spglib-1.6.0/build/lib/ added to environment

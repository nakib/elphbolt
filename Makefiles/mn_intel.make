export FFLAGS=-g -O2 -coarray=distributed
export LDFLAGS=-L/home/icn85/icn85265/spglib-1.6.0/build/lib -lsymspg
export FC=ifort
export LIBS=-mkl

#module load intel/2020.1 impi/2018.4 mkl/2017.4
#spglib-1.6.0 built with gcc and path to ~/spglib-1.6.0/build/lib/ added to environment

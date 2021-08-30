export FFLAGS=-g -O2 -fPIC -no-wrap-margin -coarray=distributed
export LDFLAGS=-L/usr/local/lib -lsymspg #v1.6.0
export FC=ifort
export LIBS=-mkl

#On msi machine need to do the following before compiling
# source ~/intel/oneapi/compiler/2021.1.1/env/vars.sh
# source ~/intel/oneapi/mkl/2021.1.1/env/vars.sh
# source ~/intel/oneapi/mpi/2021.1.1/env/vars.sh

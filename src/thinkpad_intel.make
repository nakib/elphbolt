export FFLAGS=-g -O2 -fPIC -no-wrap-margin -coarray=distributed
#export LDFLAGS=-L/usr/public/spglib/1.6.0/lib -lsymspg
export FC=ifort
export LIBS=-mkl

#libsymspg-dev and libsymspg1 version 1.14.1 on Ubuntu 20.04.2 LTS

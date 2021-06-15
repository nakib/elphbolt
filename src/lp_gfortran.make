#This works on the La Palma cluster.
export FFLAGS=-g -O2 -fbounds-check
export LDFLAGS=-L/home/icn46/icn46491/spglib-1.6.0/lib -lsymspg
export FC=caf
export LIBS=-lmkl_gf_lp64 -lmkl_core -lmkl_sequential

#module load gcc/7.2.0 openmpi/gnu/3.0.1 mkl/2018.2
#spglib-1.6.0 built with gcc and path to ~/spglib-1.6.0/lib/ added to environment
#OpenCoarrays-2.0.0 built with gcc, openmpi, and cmake and added path to environment

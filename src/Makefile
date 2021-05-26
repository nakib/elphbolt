include thinkpad_gcc.make
#include thinkpad_intel.make
#include msi_intel.make

objects:= params.o misc.o spglib_wrapper.o numerics.o crystal.o symmetry.o wannier.o delta.o electron.o phonon.o bz_sums.o interactions.o bte.o elphbolt.o

all: ../elphbolt.x
../elphbolt.x: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
delta.o delta.mod: params.o misc.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
numerics.o numerics.mod: params.mod misc.mod
crystal.o crystal.mod: params.mod misc.mod
symmetry.o symmetry.mod: params.mod misc.mod numerics.mod crystal.mod spglib_wrapper.mod
wannier.o wannier.mod: params.mod misc.mod crystal.mod numerics.mod
electron.o electron.mod: params.mod misc.mod numerics.mod crystal.mod symmetry.o wannier.mod delta.mod
phonon.o phonon.mod: params.mod misc.mod numerics.mod crystal.mod symmetry.o wannier.mod delta.mod
bz_sums.mod bz_sums.mod: params.o misc.mod crystal.mod electron.mod phonon.mod symmetry.mod
interactions.o interactions.mod: params.mod misc.mod wannier.mod crystal.mod electron.mod phonon.mod numerics.mod delta.mod
bte.o bte.mod: params.mod misc.mod numerics.mod crystal.mod symmetry.mod phonon.mod electron.mod interactions.mod bz_sums.mod
elphbolt.o elphbolt.mod: misc.mod params.mod numerics.mod crystal.mod wannier.mod symmetry.mod electron.mod phonon.mod bz_sums.mod interactions.mod bte.mod
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../*.x *.mod *.o

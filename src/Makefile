include example.make

objects:= params.o misc.o spglib_wrapper.o numerics.o migel_sc.o crystal.o symmetry.o wannier.o delta.o particle.o electron.o phonon.o bz_sums.o green.o interactions.o eliashberg.o bte.o elphbolt.o

all: ../elphbolt.x
../elphbolt.x: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
particle.o particle.mod: params.mod
delta.o delta.mod: params.mod misc.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
crystal.o crystal.mod: params.mod misc.mod
numerics.o numerics.mod: params.mod misc.mod crystal.mod
symmetry.o symmetry.mod: params.mod misc.mod numerics.mod crystal.mod spglib_wrapper.mod
wannier.o wannier.mod: params.mod misc.mod crystal.mod numerics.mod
electron.o electron.mod: params.mod particle.mod misc.mod numerics.mod crystal.mod symmetry.mod wannier.mod delta.mod
phonon.o phonon.mod: params.mod particle.mod misc.mod numerics.mod crystal.mod symmetry.mod wannier.mod delta.mod
bz_sums.o bz_sums.mod: params.mod misc.mod crystal.mod delta.mod electron.mod phonon.mod symmetry.mod
green.o green.mod: params.mod misc.mod electron.mod phonon.mod
interactions.o interactions.mod: params.mod misc.mod wannier.mod crystal.mod electron.mod phonon.mod numerics.mod delta.mod
eliashberg.o eliashberg.mod: params.mod misc.mod wannier.mod crystal.mod electron.mod phonon.mod numerics.mod delta.mod
migel_sc.o migel_sc.mod: params.mod misc.mod numerics.mod wannier.mod electron.mod eliashberg.mod
bte.o bte.mod: params.mod misc.mod numerics.mod crystal.mod symmetry.mod phonon.mod electron.mod interactions.mod bz_sums.mod
elphbolt.o elphbolt.mod: misc.mod params.mod numerics.mod crystal.mod wannier.mod symmetry.mod electron.mod phonon.mod bz_sums.mod migel_sc.o eliashberg.o interactions.mod bte.mod

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../*.x *.mod *.o

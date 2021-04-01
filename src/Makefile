include thinkpad_arch.make

objects:= params.o misc.o spglib_wrapper.o numerics.o crystal.o symmetries.o wannier.o delta.o electrons.o phonons.o interactions.o elphbolt.o

all: ../elphbolt
../elphbolt: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
delta.o delta.mod: params.o misc.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
numerics.o numerics.mod: params.mod misc.mod
crystal.o crystal.mod: params.mod misc.mod
symmetries.o symmetries.mod: params.mod misc.mod numerics.mod crystal.mod spglib_wrapper.mod
wannier.o wannier.mod: params.mod misc.mod crystal.mod numerics.mod
electrons.o electrons.mod: params.mod misc.mod numerics.mod crystal.mod symmetries.o wannier.mod delta.mod
phonons.o phonons.mod: params.mod misc.mod numerics.mod crystal.mod symmetries.o wannier.mod
interactions.o interactions.mod: params.mod misc.mod wannier.mod crystal.mod electrons.mod phonons.mod numerics.mod
elphbolt.o elphbolt.mod: params.mod numerics.mod crystal.mod wannier.mod symmetries.mod electrons.mod phonons.mod interactions.mod
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../elphbolt *.mod *.o

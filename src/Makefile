include thinkpad_arch.make

#objects:= params.o data.o misc.o spglib_wrapper.o config.o wannier.o mesh.o elphbolt.o
#objects:= params.o data.o misc.o spglib_wrapper.o config.o wannier.o elphbolt.o
objects:= params.o misc.o spglib_wrapper.o numerics.o crystal.o symmetries.o wannier.o elphbolt.o

all: ../elphbolt
../elphbolt: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
numerics.o numerics.mod: params.mod misc.mod
crystal.o crystal.mod: params.mod misc.mod
wannier.o wannier.mod: params.mod misc.mod crystal.mod numerics.mod
symmetries.o symmetries.mod: params.mod misc.mod numerics.mod crystal.mod spglib_wrapper.mod
elphbolt.o elphbolt.mod: params.mod numerics.mod crystal.mod wannier.mod symmetries.mod
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../elphbolt *.mod *.o

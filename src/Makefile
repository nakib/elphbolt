include thinkpad_arch.make

objects:= params.o derived_types.o misc.o spglib_wrapper.o config.o wannier.o elphbolt.o

all: ../elphbolt
../elphbolt: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
derived_types.o derived_types.mod: params.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
wannier.o wannier.mod: params.mod derived_types.mod misc.mod
config.o config.mod: params.mod misc.mod derived_types.mod spglib_wrapper.mod
elphbolt.o elphbolt.mod: params.mod derived_types.mod config.mod wannier.mod
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../elphbolt *.mod *.o

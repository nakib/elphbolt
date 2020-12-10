include thinkpad_arch.make

objects:= params.o data.o misc.o spglib_wrapper.o config.o wannier.o mesh.o elphbolt.o

all: ../elphbolt
../elphbolt: $(objects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
misc.o misc.mod: params.mod
spglib_wrapper.o spglib_wrapper.mod: params.mod
data.o data.mod: params.mod
mesh.o mesh.mod: params.mod data.mod misc.mod wannier.mod
wannier.o wannier.mod: params.mod data.mod misc.mod
config.o config.mod: params.mod misc.mod data.mod spglib_wrapper.mod wannier.mod
elphbolt.o elphbolt.mod: params.mod data.mod config.mod mesh.mod wannier.mod
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../elphbolt *.mod *.o

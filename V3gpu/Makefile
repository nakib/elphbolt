# libgpu_v3.so built from gpu_3ph_kernels.cuf, gpu_driver.f90, gpu_wrapper.f90
CUDA_HOME ?= /usr/local/cuda

.PHONY: all build-so clean

all: build-so

build-so:
	nvfortran -shared -cuda -gpu=cc89 -fPIC \
	  -o libgpu_v3.so \
	  gpu_3ph_kernels.f90 gpu_driver.f90 gpu_wrapper.f90 \
	  -L"$(CUDA_HOME)/lib64" -lcudart -Wl,-rpath,"$(CUDA_HOME)/lib64"

clean:
	rm -f libgpu_v3.so *.o *.mod


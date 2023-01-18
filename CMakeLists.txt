#Set minimum version of CMake
cmake_minimum_required(VERSION 3.10)

#Set project name and language
project(elphbolt Fortran)

#Set compiler
set(CMAKE_Fortran_COMPILER "ifort")
#set(CMAKE_Fortran_COMPILER "caf")

#Set flags and mkl library for ifort
if(CMAKE_Fortran_COMPILER MATCHES "ifort")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2 -fPIC -warn -no-wrap-margin -coarray=distributed")

  find_library(MKL_LIB REQUIRED)
endif()

#Set flags and blas+lapack libraries for caf
if(CMAKE_Fortran_COMPILER MATCHES "caf")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2 -fPIC -fbounds-check -Wunused -Wconversion -Wunderflow -Wdo-subscript")

  find_library(BLAS_LIB NAMES blas REQUIRED)
  find_library(LAPACK_LIB NAMES lapack REQUIRED)
endif()

#Set spg library
find_library(SPG_LIB NAMES symspg REQUIRED)

#Set directory for executable binary/ies
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#Set source directory/ies
add_subdirectory(src)


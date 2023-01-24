#Set minimum version of CMake
cmake_minimum_required(VERSION 3.10)

#Set project name and language
project(elphbolt Fortran)

#Set directory for executable binary/ies
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#Set compiler
#set(CMAKE_Fortran_COMPILER "ifort")
set(CMAKE_Fortran_COMPILER "caf")

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

#Set compiler dependent library links
#For ifort
if(CMAKE_Fortran_COMPILER MATCHES "ifort")
  set(LINKS "-mkl ${SPG_LIB}")
endif()

#For caf
if(CMAKE_Fortran_COMPILER MATCHES "caf")
  set(LINKS "-lblas -llapack ${SPG_LIB}")
endif()

#Set elphbolt source
set(SRC
  src/params.f90
  src/misc.f90
  src/spglib_wrapper.f90
  src/numerics.f90
  src/migel_sc.f90
  src/crystal.f90
  src/symmetry.f90
  src/wannier.f90
  src/delta.f90
  src/particle.f90
  src/electron.f90
  src/phonon.f90
  src/bz_sums.f90
  src/green.f90
  src/phonon_defect.f90
  src/interactions.f90
  src/eliashberg.f90
  src/bte.f90)

#Set elphbolt driver program
add_executable(elphbolt app/elphbolt.f90 ${SRC})

#Link libraries
target_link_libraries(elphbolt ${LINKS})

#Set test library
add_library(testify thirdparty/testify/src/testify.f90)

enable_testing()

#Unit tests
foreach(component misc)
  #Set test program
  add_executable(test_${component} test/test_${component}.f90 ${SRC})

  #Link libraries
  target_link_libraries(test_${component} testify ${LINKS})

  #Add to test battery
  add_test(test_${component} bin/test_${component})
endforeach()

#Regression tests
find_program (bash /bin/bash)
if(bash)
  #For ifort
  if(CMAKE_Fortran_COMPILER MATCHES "ifort")
    add_test (bte_regression ${script} ../test/3C-SiC/run_bte_ifort.sh)
  endif()

  #For caf
  if(CMAKE_Fortran_COMPILER MATCHES "caf")
    add_test (bte_regression ${script} ../test/3C-SiC/run_bte_caf.sh)
  endif()
endif(bash)

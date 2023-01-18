include_directories(${CMAKE_BINARY_DIR}/src)

add_executable(elphbolt
  params.f90
  misc.f90
  spglib_wrapper.f90
  numerics.f90
  migel_sc.f90
  crystal.f90
  symmetry.f90
  wannier.f90
  delta.f90
  particle.f90
  electron.f90
  phonon.f90
  bz_sums.f90
  green.f90
  phonon_defect.f90
  interactions.f90
  eliashberg.f90
  bte.f90
  elphbolt.f90
  )

#For ifort
if(CMAKE_Fortran_COMPILER MATCHES "ifort")
  target_link_libraries(elphbolt -mkl "${SPG_LIB}")
endif()

#For caf
if(CMAKE_Fortran_COMPILER MATCHES "caf")
  target_link_libraries(elphbolt lapack blas "${SPG_LIB}")
endif()


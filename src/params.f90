module params
  !! Module containing various parameters and constants.
  
  implicit none

  integer, parameter :: dp = selected_real_kind(14,200)
  !! Custom real precision.
  integer, parameter :: k8 = selected_real_kind(8)
  !! Custom integer precision.

  !Physical constants:
  real(dp), parameter :: qe = 1.602176634e-19_dp
  !! Electron charge magnitude (C)
  real(dp), parameter :: me = 9.1093837015e-31_dp
  !! Electron mass (Kg)
  real(dp), parameter :: amu = 1.66053906660e-27_dp
  !! Atomic mass unit (Kg)
  real(dp), parameter :: hbar = 1.05457172647e-22_dp
  !! Reduced Planck's constant (J/THz = J.ps)
  real(dp), parameter :: hbar_eVps = hbar/qe
  !! Reduced Planck's constant (eV/THz = eV.ps)
  real(dp), parameter :: perm0 =  8.854187817e-12_dp
  !! Permittivity of free space (F/m)
  real(dp), parameter :: kB = 1.380649e-23_dp/qe
  !! Boltzmann constant (eV/K)
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  !! Value of pi
  real(dp), parameter :: twopi = 2.0_dp*pi
  !! Value of 2pi

  !Conversion factors:
  real(dp), parameter :: Hartree2radTHz = 27.2116_dp*qe/hbar !Hartree to rad.THz
  real(dp), parameter :: Hartree2eV = 27.2116_dp !Hartree to eV
  real(dp), parameter :: Ryd2radTHz = 0.5_dp*Hartree2radTHz !Rydberg to rad.THz
  real(dp), parameter :: Ryd2eV = 0.5_dp*Hartree2eV !Rydberg to eV
  real(dp), parameter :: Ryd2meV = Ryd2eV*1.0e3_dp !Rydberg to meV
  real(dp), parameter :: Ryd2amu = 2.0_dp*me/amu !Rydberg mass to amu
  real(dp), parameter :: bohr2nm=0.052917721092_dp !Bohr to nm
  
  !Miscellaneous
  complex(dp), parameter :: oneI = (0.0_dp,1.0_dp)
  complex(dp), parameter :: twopiI = twopi*oneI
end module params

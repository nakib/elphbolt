! Copyright 2020 elphbolt contributors.
! This file is part of elphbolt <https://github.com/nakib/elphbolt>.
!
! elphbolt is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! elphbolt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with elphbolt. If not, see <http://www.gnu.org/licenses/>.

module params
  !! Module containing various parameters and constants.

  use iso_fortran_env, only: real64, real128, int64
  
  implicit none

  integer, parameter :: r64 = real64
  integer, parameter :: r128 = real128
  integer, parameter :: i64 = int64

  !Physical constants:
  real(r64), parameter :: qe = 1.602176634e-19_r64
  !! Electron charge magnitude (C)
  real(r64), parameter :: me = 9.1093837015e-31_r64
  !! Electron mass (Kg)
  real(r64), parameter :: amu = 1.66053906660e-27_r64
  !! Atomic mass unit (Kg)
  real(r64), parameter :: hbar = 1.05457172647e-22_r64
  !! Reduced Planck's constant (J/THz = J.ps)
  real(r64), parameter :: hbar_eVps = hbar/qe
  !! Reduced Planck's constant (eV/THz = eV.ps)
  real(r64), parameter :: perm0 =  8.854187817e-12_r64
  !! Permittivity of free space (F/m)
  real(r64), parameter :: kB = 1.380649e-23_r64/qe
  !! Boltzmann constant (eV/K)
  real(r64), parameter :: pi = 4.0_r64*atan(1.0_r64)
  !! Value of pi
  real(r64), parameter :: twopi = 2.0_r64*pi
  !! Value of 2pi

  !Conversion factors:
  real(r64), parameter :: Hartree2radTHz = 27.2116_r64*qe/hbar !Hartree to rad.THz
  real(r64), parameter :: Hartree2eV = 27.2116_r64 !Hartree to eV
  real(r64), parameter :: Ryd2radTHz = 0.5_r64*Hartree2radTHz !Rydberg to rad.THz
  real(r64), parameter :: Ryd2eV = 0.5_r64*Hartree2eV !Rydberg to eV
  real(r64), parameter :: Ryd2meV = Ryd2eV*1.0e3_r64 !Rydberg to meV
  real(r64), parameter :: Ryd2amu = 2.0_r64*me/amu !Rydberg mass to amu
  real(r64), parameter :: bohr2nm=0.052917721092_r64 !Bohr to nm
  
  !Miscellaneous
  complex(r64), parameter :: oneI = (0.0_r64,1.0_r64)
  complex(r64), parameter :: twopiI = twopi*oneI

  !The code below is adapted from ShengBTE (file data.f90):
  
  !Periodic table from ShengBTE
  character(len = 3), parameter :: periodic_table(114)=[character(len=3) :: &
       "H","He","Li","Be","B","C","N","O", &
       "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V", &
       "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb", &
       "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb", &
       "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb", &
       "Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg", &
       "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am", &
       "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds" ,&
       "Rg","Cn","Uuq","Uuh"]
end module params

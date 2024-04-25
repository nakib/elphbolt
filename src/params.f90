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

  use precision
  use fhash, only: fhash_tbl_t, key => fhash_key
  use isotopes_module, only: isotopes

  implicit none

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
  real(r64), parameter :: perm0 = 8.854187817e-12_r64
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
  real(r64), parameter :: bohr2nm = 0.052917721092_r64 !Bohr to nm
  
  !Miscellaneous
  complex(r64), parameter :: oneI = (0.0_r64,1.0_r64)
  complex(r64), parameter :: twopiI = twopi*oneI

  !Not parameters, but effectively behaves as parameters
  logical, private :: periodic_table_exists = .false.
  type(fhash_tbl_t), private :: periodic_table

contains

  function lookup_periodic_table(element) result(data)
      character(*), intent(in) :: element
      type(isotopes) :: data

      class(*), allocatable :: raw_data

      if (.not. periodic_table_exists) then
         call create_periodic_table(periodic_table)
         periodic_table_exists = .true.
      end if

      call periodic_table%get_raw(key(element), raw_data)

      select type (raw_data)
      type is (isotopes)
         data = raw_data
      class default
         print *, 'Error in reading periodic table for element ' // element // '.'
      end select
   end function lookup_periodic_table

   subroutine create_periodic_table(periodic_table)
      type(fhash_tbl_t), intent(out) :: periodic_table

      call periodic_table%set( & 
         key('Ag'), value=isotopes( & 
            [106.905095_r64, 108.904754_r64], & 
            [51.84_r64, 48.16_r64]))
      call periodic_table%set( & 
         key('Al'), value=isotopes( & 
            [26.981541_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Ar'), value=isotopes( & 
            [35.967546_r64, 37.962732_r64, 39.962383_r64], & 
            [0.34_r64, 0.063_r64, 99.6_r64]))
      call periodic_table%set( & 
         key('As'), value=isotopes( & 
            [74.921596_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Au'), value=isotopes( & 
            [196.96656_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('B'), value=isotopes( & 
            [10.012938_r64, 11.009305_r64], & 
            [19.8_r64, 80.2_r64]))
      call periodic_table%set( & 
         key('Ba'), value=isotopes( & 
            [129.906277_r64, 131.905042_r64, 133.90449_r64, 134.905668_r64, 135.904556_r64, 136.905816_r64, 137.905236_r64], & 
            [0.11_r64, 0.1_r64, 2.42_r64, 6.59_r64, 7.85_r64, 11.23_r64, 71.7_r64]))
      call periodic_table%set( & 
         key('Be'), value=isotopes( & 
            [9.012183_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Bi'), value=isotopes( & 
            [208.980388_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Br'), value=isotopes( & 
            [78.918336_r64, 80.91629_r64], & 
            [50.69_r64, 49.31_r64]))
      call periodic_table%set( & 
         key('C'), value=isotopes( & 
            [12.0_r64, 13.003355_r64], & 
            [98.9_r64, 1.1_r64]))
      call periodic_table%set( & 
         key('Ca'), value=isotopes( & 
            [39.962591_r64, 41.958622_r64, 42.95877_r64, 43.955485_r64, 45.953689_r64, 47.952532_r64], & 
            [96.95_r64, 0.65_r64, 0.14_r64, 2.086_r64, 0.004_r64, 0.19_r64]))
      call periodic_table%set( & 
         key('Cd'), value=isotopes( & 
            [105.906461_r64, 107.904186_r64, 109.903007_r64, 110.904182_r64, 111.902761_r64, 112.904401_r64, 113.903361_r64, & 
             115.904758_r64], & 
            [1.25_r64, 0.89_r64, 12.49_r64, 12.8_r64, 24.13_r64, 12.22_r64, 28.73_r64, 7.49_r64]))
      call periodic_table%set( & 
         key('Ce'), value=isotopes( & 
            [135.90714_r64, 137.905996_r64, 139.905442_r64, 141.909249_r64], & 
            [0.19_r64, 0.25_r64, 88.48_r64, 11.08_r64]))
      call periodic_table%set( & 
         key('Cl'), value=isotopes( & 
            [34.968853_r64, 36.965903_r64], & 
            [75.77_r64, 24.23_r64]))
      call periodic_table%set( & 
         key('Co'), value=isotopes( & 
            [58.933198_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Cr'), value=isotopes( & 
            [49.946046_r64, 51.94051_r64, 52.940651_r64, 53.938882_r64], & 
            [4.35_r64, 83.79_r64, 9.5_r64, 2.36_r64]))
      call periodic_table%set( & 
         key('Cs'), value=isotopes( & 
            [132.905433_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Cu'), value=isotopes( & 
            [62.929599_r64, 64.927792_r64], & 
            [69.17_r64, 30.83_r64]))
      call periodic_table%set( & 
         key('Dy'), value=isotopes( & 
            [155.924287_r64, 157.924412_r64, 159.925203_r64, 160.926939_r64, 161.926805_r64, 162.928737_r64, 163.929183_r64], & 
            [0.06_r64, 0.1_r64, 2.34_r64, 18.9_r64, 25.5_r64, 24.9_r64, 28.2_r64]))
      call periodic_table%set( & 
         key('Er'), value=isotopes( & 
            [161.928787_r64, 163.929211_r64, 165.930305_r64, 166.932061_r64, 167.932383_r64, 169.935476_r64], & 
            [0.14_r64, 1.61_r64, 33.6_r64, 22.95_r64, 26.8_r64, 14.9_r64]))
      call periodic_table%set( & 
         key('Eu'), value=isotopes( & 
            [150.91986_r64, 152.921243_r64], & 
            [47.8_r64, 52.2_r64]))
      call periodic_table%set( & 
         key('F'), value=isotopes( & 
            [18.998403_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Fe'), value=isotopes( & 
            [53.939612_r64, 55.934939_r64, 56.935396_r64, 57.933278_r64], & 
            [5.8_r64, 91.72_r64, 2.2_r64, 0.28_r64]))
      call periodic_table%set( & 
         key('Ga'), value=isotopes( & 
            [68.925581_r64, 70.924701_r64], & 
            [60.1_r64, 39.9_r64]))
      call periodic_table%set( & 
         key('Gd'), value=isotopes( & 
            [151.919803_r64, 153.920876_r64, 154.822629_r64, 155.92213_r64, 156.923967_r64, 157.924111_r64, 159.927061_r64], & 
            [0.2_r64, 2.18_r64, 14.8_r64, 20.47_r64, 15.65_r64, 24.84_r64, 21.86_r64]))
      call periodic_table%set( & 
         key('Ge'), value=isotopes( & 
            [69.92425_r64, 71.92208_r64, 72.923464_r64, 73.921179_r64, 75.921403_r64], & 
            [20.5_r64, 27.4_r64, 7.8_r64, 36.5_r64, 7.8_r64]))
      call periodic_table%set( & 
         key('H'), value=isotopes( & 
            [1.007825_r64, 2.014102_r64], & 
            [99.99_r64, 0.015_r64]))
      call periodic_table%set( & 
         key('He'), value=isotopes( & 
            [3.016029_r64, 4.002603_r64], & 
            [0.0001_r64, 100.0_r64]))
      call periodic_table%set( & 
         key('Hf'), value=isotopes( & 
            [173.940065_r64, 175.94142_r64, 176.943233_r64, 177.94371_r64, 178.945827_r64, 179.946561_r64], & 
            [0.16_r64, 5.2_r64, 18.6_r64, 27.1_r64, 13.74_r64, 35.2_r64]))
      call periodic_table%set( & 
         key('Hg'), value=isotopes( & 
            [195.965812_r64, 197.96676_r64, 198.968269_r64, 199.968316_r64, 200.970293_r64, 201.970632_r64, 203.973481_r64], & 
            [0.15_r64, 10.1_r64, 17.0_r64, 23.1_r64, 13.2_r64, 29.65_r64, 6.8_r64]))
      call periodic_table%set( & 
         key('Ho'), value=isotopes( & 
            [164.930332_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('I'), value=isotopes( & 
            [126.904477_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('In'), value=isotopes( & 
            [112.904056_r64, 114.903875_r64], & 
            [4.3_r64, 95.7_r64]))
      call periodic_table%set( & 
         key('Ir'), value=isotopes( & 
            [190.960603_r64, 192.962942_r64], & 
            [37.3_r64, 62.7_r64]))
      call periodic_table%set( & 
         key('K'), value=isotopes( & 
            [38.963708_r64, 39.963999_r64, 40.961825_r64], & 
            [93.2_r64, 0.012_r64, 6.73_r64]))
      call periodic_table%set( & 
         key('Kr'), value=isotopes( & 
            [77.920397_r64, 79.916375_r64, 81.913483_r64, 82.914134_r64, 83.911506_r64, 85.910614_r64], & 
            [0.35_r64, 2.25_r64, 11.6_r64, 11.5_r64, 57.0_r64, 17.3_r64]))
      call periodic_table%set( & 
         key('La'), value=isotopes( & 
            [137.907114_r64, 138.906355_r64], & 
            [0.09_r64, 99.91_r64]))
      call periodic_table%set( & 
         key('Li'), value=isotopes( & 
            [6.015123_r64, 7.016005_r64], & 
            [7.42_r64, 92.58_r64]))
      call periodic_table%set( & 
         key('Lu'), value=isotopes( & 
            [174.940785_r64, 175.942694_r64], & 
            [97.4_r64, 2.6_r64]))
      call periodic_table%set( & 
         key('Mg'), value=isotopes( & 
            [23.985045_r64, 24.985839_r64, 25.982595_r64], & 
            [78.9_r64, 10.0_r64, 11.1_r64]))
      call periodic_table%set( & 
         key('Mn'), value=isotopes( & 
            [54.938046_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Mo'), value=isotopes( & 
            [91.906809_r64, 93.905086_r64, 94.905838_r64, 95.904676_r64, 96.906018_r64, 97.905405_r64, 99.907473_r64], & 
            [14.84_r64, 9.25_r64, 15.92_r64, 16.68_r64, 9.55_r64, 24.13_r64, 9.63_r64]))
      call periodic_table%set( & 
         key('N'), value=isotopes( & 
            [14.003074_r64, 15.000109_r64], & 
            [99.63_r64, 0.37_r64]))
      call periodic_table%set( & 
         key('Na'), value=isotopes( & 
            [22.98977_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Nb'), value=isotopes( & 
            [92.906378_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Nd'), value=isotopes( & 
            [141.907731_r64, 142.909823_r64, 143.910096_r64, 144.912582_r64, 145.913126_r64, 147.916901_r64, 149.9209_r64], & 
            [27.13_r64, 12.18_r64, 23.8_r64, 8.3_r64, 17.19_r64, 5.76_r64, 5.64_r64]))
      call periodic_table%set( & 
         key('Ne'), value=isotopes( & 
            [19.992439_r64, 20.993845_r64, 21.991384_r64], & 
            [90.6_r64, 0.26_r64, 9.2_r64]))
      call periodic_table%set( & 
         key('Ni'), value=isotopes( & 
            [57.935347_r64, 59.930789_r64, 60.931059_r64, 61.928346_r64, 63.927968_r64], & 
            [68.27_r64, 26.1_r64, 1.13_r64, 3.59_r64, 0.91_r64]))
      call periodic_table%set( & 
         key('O'), value=isotopes( & 
            [15.994915_r64, 16.999131_r64, 17.999159_r64], & 
            [99.76_r64, 0.038_r64, 0.2_r64]))
      call periodic_table%set( & 
         key('Os'), value=isotopes( & 
            [183.952514_r64, 185.953852_r64, 186.955762_r64, 187.95585_r64, 188.958156_r64, 189.958455_r64, 191.961487_r64], & 
            [0.02_r64, 1.58_r64, 1.6_r64, 13.3_r64, 16.1_r64, 26.4_r64, 41.0_r64]))
      call periodic_table%set( & 
         key('P'), value=isotopes( & 
            [30.973763_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Pb'), value=isotopes( & 
            [203.973037_r64, 205.974455_r64, 206.975885_r64, 207.976641_r64], & 
            [1.4_r64, 24.1_r64, 22.1_r64, 52.4_r64]))
      call periodic_table%set( & 
         key('Pd'), value=isotopes( & 
            [101.905609_r64, 103.904026_r64, 104.905075_r64, 105.903475_r64, 107.903894_r64, 109.905169_r64], & 
            [1.02_r64, 11.14_r64, 22.33_r64, 27.33_r64, 26.46_r64, 11.72_r64]))
      call periodic_table%set( & 
         key('Pr'), value=isotopes( & 
            [140.907657_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Pt'), value=isotopes( & 
            [189.959937_r64, 191.961049_r64, 193.962679_r64, 194.964785_r64, 195.964947_r64, 197.967879_r64], & 
            [0.01_r64, 0.79_r64, 32.9_r64, 33.8_r64, 25.3_r64, 7.2_r64]))
      call periodic_table%set( & 
         key('Rb'), value=isotopes( & 
            [84.9118_r64, 86.909184_r64], & 
            [72.17_r64, 27.84_r64]))
      call periodic_table%set( & 
         key('Re'), value=isotopes( & 
            [184.952977_r64, 186.955765_r64], & 
            [37.4_r64, 62.6_r64]))
      call periodic_table%set( & 
         key('Rh'), value=isotopes( & 
            [102.905503_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Ru'), value=isotopes( & 
            [95.907596_r64, 97.905287_r64, 98.905937_r64, 99.904218_r64, 100.905581_r64, 101.904348_r64, 103.905422_r64], & 
            [5.52_r64, 1.88_r64, 12.7_r64, 12.6_r64, 17.0_r64, 31.6_r64, 18.7_r64]))
      call periodic_table%set( & 
         key('S'), value=isotopes( & 
            [31.972072_r64, 32.971459_r64, 33.967868_r64, 35.967079_r64], & 
            [95.02_r64, 0.75_r64, 4.21_r64, 0.02_r64]))
      call periodic_table%set( & 
         key('Sb'), value=isotopes( & 
            [120.903824_r64, 122.904222_r64], & 
            [57.3_r64, 42.7_r64]))
      call periodic_table%set( & 
         key('Sc'), value=isotopes( & 
            [44.955914_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Se'), value=isotopes( & 
            [73.922477_r64, 75.919207_r64, 76.919908_r64, 77.917304_r64, 79.916521_r64, 81.916709_r64], & 
            [0.9_r64, 9.0_r64, 7.6_r64, 23.5_r64, 49.6_r64, 9.4_r64]))
      call periodic_table%set( & 
         key('Si'), value=isotopes( & 
            [27.976928_r64, 28.976496_r64, 29.973772_r64], & 
            [92.23_r64, 4.67_r64, 3.1_r64]))
      call periodic_table%set( & 
         key('Sm'), value=isotopes( & 
            [143.912009_r64, 146.914907_r64, 147.914832_r64, 148.917193_r64, 149.917285_r64, 151.919741_r64, 153.922218_r64], & 
            [3.1_r64, 15.0_r64, 11.3_r64, 13.8_r64, 7.4_r64, 26.7_r64, 22.7_r64]))
      call periodic_table%set( & 
         key('Sn'), value=isotopes( & 
            [111.904826_r64, 113.902784_r64, 114.903348_r64, 115.901744_r64, 116.902954_r64, 117.901607_r64, 118.90331_r64, & 
             119.902199_r64, 121.90344_r64, 123.905271_r64], & 
            [0.97_r64, 0.65_r64, 0.36_r64, 14.7_r64, 7.7_r64, 24.3_r64, 8.6_r64, 32.4_r64, 4.6_r64, 5.6_r64]))
      call periodic_table%set( & 
         key('Sr'), value=isotopes( & 
            [83.913428_r64, 85.909273_r64, 86.908902_r64, 87.905625_r64], & 
            [0.56_r64, 9.86_r64, 7.0_r64, 82.58_r64]))
      call periodic_table%set( & 
         key('Ta'), value=isotopes( & 
            [179.947489_r64, 180.948014_r64], & 
            [0.012_r64, 99.99_r64]))
      call periodic_table%set( & 
         key('Tb'), value=isotopes( & 
            [158.92535_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Te'), value=isotopes( & 
            [119.904021_r64, 121.903055_r64, 122.904278_r64, 123.902825_r64, 124.904435_r64, 125.90331_r64, 127.904464_r64, &
             129.906229_r64], & 
            [0.096_r64, 2.6_r64, 0.91_r64, 4.82_r64, 7.14_r64, 18.95_r64, 31.69_r64, 33.8_r64]))
      call periodic_table%set( & 
         key('Th'), value=isotopes( & 
            [232.038054_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Ti'), value=isotopes( & 
            [45.952633_r64, 46.951765_r64, 47.947947_r64, 48.947871_r64, 49.944786_r64], & 
            [8.0_r64, 7.3_r64, 73.8_r64, 5.5_r64, 5.4_r64]))
      call periodic_table%set( & 
         key('Tl'), value=isotopes( & 
            [202.972336_r64, 204.97441_r64], & 
            [29.52_r64, 70.48_r64]))
      call periodic_table%set( & 
         key('Tm'), value=isotopes( & 
            [168.934225_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('U'), value=isotopes( & 
            [234.040947_r64, 235.043925_r64, 238.050786_r64], & 
            [0.006_r64, 0.72_r64, 99.27_r64]))
      call periodic_table%set( & 
         key('V'), value=isotopes( & 
            [49.947161_r64, 50.943963_r64], & 
            [0.25_r64, 99.75_r64]))
      call periodic_table%set( & 
         key('W'), value=isotopes( & 
            [179.946727_r64, 181.948225_r64, 182.950245_r64, 183.950953_r64, 185.954377_r64], & 
            [0.13_r64, 26.3_r64, 14.3_r64, 30.67_r64, 28.6_r64]))
      call periodic_table%set( & 
         key('Xe'), value=isotopes( & 
            [123.905894_r64, 125.904281_r64, 127.903531_r64, 128.90478_r64, 129.90351_r64, 130.905076_r64, 131.904148_r64, &
             133.905395_r64, 135.907219_r64], & 
            [0.1_r64, 0.09_r64, 1.91_r64, 26.4_r64, 4.1_r64, 21.2_r64, 26.9_r64, 10.4_r64, 8.9_r64]))
      call periodic_table%set( & 
         key('Y'), value=isotopes( & 
            [88.905856_r64], & 
            [100.0_r64]))
      call periodic_table%set( & 
         key('Yb'), value=isotopes( & 
            [167.933908_r64, 169.934774_r64, 170.936338_r64, 171.936393_r64, 172.938222_r64, 173.938873_r64, 175.942576_r64], & 
            [0.13_r64, 3.05_r64, 14.3_r64, 21.9_r64, 16.12_r64, 31.8_r64, 12.7_r64]))
      call periodic_table%set( & 
         key('Zn'), value=isotopes( & 
            [63.929145_r64, 65.926035_r64, 66.927129_r64, 67.924846_r64, 69.925325_r64], & 
            [48.6_r64, 27.9_r64, 4.1_r64, 18.8_r64, 0.6_r64]))
      call periodic_table%set( & 
         key('Zr'), value=isotopes( & 
            [89.904708_r64, 90.905644_r64, 91.905039_r64, 93.906319_r64, 95.908272_r64], & 
            [51.45_r64, 11.27_r64, 17.17_r64, 17.33_r64, 2.78_r64]))

   end subroutine create_periodic_table
end module params

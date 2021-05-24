! Copyright (C) 2020- Nakib Haider Protik <nakib.haider.protik@gmail.com>
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

module crystal_module
  !! Module containing type and procedures related to the crystal structure.

  use params, only: dp, k8, twopi
  use misc, only: exit_with_message, print_message, cross_product, demux_vector, &
       subtitle

  implicit none

  private
  public crystal, calculate_wavevectors_full

  type crystal
     !! Data and procedures related to the crystal structure.

     integer(k8) :: numelements
     !! Number of types of basis atoms.
     integer(k8) :: numatoms
     !! Number of basis atoms.
     character(len=100) :: name
     !! Name of material.
     character(len=3), allocatable :: elements(:)
     !! Elements in the in the basis.
     integer(k8), allocatable :: atomtypes(:)
     !! Integer tagging unique elements in the basis.
     real(dp), allocatable :: masses(:)
     !! Masses of the basis atoms.
     logical :: polar
     !! Is the system polar?
     real(dp) :: epsilon(3,3)
     !! Dielectric tensor
     real(dp), allocatable :: born(:,:,:)
     !! Born effective charge
     real(dp), allocatable :: basis(:,:)
     !! Basis vectors (crystal coordinates).
     real(dp), allocatable :: basis_cart(:,:)
     !! Basis vectors (Cartesian coordinates).
     real(dp) :: lattvecs(3,3)
     !! Lattice vectors (nm).
     real(dp) :: volume
     !! Volume of primitive cell (nm^3).
     real(dp) :: reclattvecs(3,3)
     !! Reciprocal lattice vectors.
     real(dp) :: volume_bz
     !! Brillouin zone volume (nm^-3).
     real(dp) :: T
     !! Crystal temperature (T).
     logical :: autoisotopes
     !! Use isotopic mix for masses?
     real(dp), allocatable :: gfactors(:)
     !! g-factors.
     
   contains

     procedure :: initialize=>read_input_and_setup_crystal
  end type crystal

contains

  subroutine read_input_and_setup_crystal(c)
    !! Read input file and initialize crystal data.

    class(crystal), intent(out) :: c

    !Local variables
    integer(k8) :: i, j, k, numelements, numatoms
    integer(k8), allocatable :: atomtypes(:)
    real(dp), allocatable :: masses(:), gfactors(:), born(:,:,:), basis(:,:), basis_cart(:,:)
    real(dp) :: epsilon(3,3), lattvecs(3,3), volume, reclattvecs(3,3), volume_bz, T
    character(len=3), allocatable :: elements(:)
    character(len=100) :: name
    logical :: polar, autoisotopes, phiso
    
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ name, elements, atomtypes, basis, lattvecs, &
         polar, born, epsilon, masses, T, autoisotopes, phiso

    call subtitle("Setting up crystal...")

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Set values from input:
    ! Read allocations
    read(1, nml = allocations)
    if(numelements < 1 .or. numatoms < 1 .or. numatoms > numelements) then
       call exit_with_message('Bad input(s) in allocations.')
    end if
    c%numelements = numelements
    c%numatoms = numatoms

    !Allocate variables
    allocate(elements(numelements), atomtypes(numatoms), born(3,3,numatoms), &
         basis(3,numatoms), masses(numelements), basis_cart(3,numatoms))
    allocate(c%elements(c%numelements), c%atomtypes(c%numatoms), c%born(3,3,c%numatoms), &
         c%masses(c%numatoms), c%gfactors(c%numelements), c%basis(3,c%numatoms), &
         c%basis_cart(3,c%numatoms))
    
    !Read crystal_info
    autoisotopes = .true.
    polar = .false.
    epsilon = 0.0_dp
    born = 0.0_dp
    T = -1.0_dp
    read(1, nml = crystal_info)
    if(any(atomtypes < 1) .or. any(masses < 0) .or. T < 0.0_dp) then
       call exit_with_message('Bad input(s) in crystal_info.')
    end if

    !Close input file
    close(1)
    
    c%name = name
    c%elements = elements
    c%atomtypes = atomtypes
    c%born = born
    c%epsilon = epsilon
    c%basis = basis
    c%polar = polar
    c%lattvecs = lattvecs
    c%T = T
    c%autoisotopes = autoisotopes
    c%masses = masses
    c%gfactors = 0.0_dp

    !If required, calculate isotopic average masses and g-factors
    if(autoisotopes) then
       call calculate_mavg_and_g(c%elements, c%masses, c%gfactors)
    end if

    !Calculate atomic basis in Cartesian coordinates
    c%basis_cart(:,:) = matmul(c%lattvecs,c%basis)

    !Calculate reciprocal lattice vectors and real and reciprocal cell volumes
    do i = 1,3
       j = mod(i,3) + 1
       k = mod(j,3) + 1
       c%reclattvecs(:,i) = &
            cross_product(c%lattvecs(:,j), c%lattvecs(:,k))
    end do
    c%volume = abs(dot_product(c%lattvecs(:,1),c%reclattvecs(:,1)))
    c%volume_bz = twopi/c%volume
    c%reclattvecs(:,:) = c%volume_bz*c%reclattvecs(:,:)

    !Print out crystal and reciprocal lattice information.
    if(this_image() == 1) then
       write(*, "(A, A)") 'Material: ', c%name
       if(c%autoisotopes) write(*,"(A)") 'Isotopic average of masses will be used.'
       do i = 1, c%numatoms
          write(*,"(A, A, 1E16.8, A)") trim(c%elements(i)), " mass = ", c%masses(i), " u"
       end do
       write(*,"(A)") 'Lattice vectors [nm]:'
       write(*,"(3(1E16.8,x))") c%lattvecs(:,1)
       write(*,"(3(1E16.8,x))") c%lattvecs(:,2)
       write(*,"(3(1E16.8,x))") c%lattvecs(:,3)
       write(*,"(A,(1E16.8,x),A)") 'Primitive cell volume =', c%volume, 'nm^3'

       write(*,"(A)") 'Reciprocal lattice vectors [1/nm]:'
       write(*,"(3(1E16.8,x))") c%reclattvecs(:,1)
       write(*,"(3(1E16.8,x))") c%reclattvecs(:,2)
       write(*,"(3(1E16.8,x))") c%reclattvecs(:,3)
       write(*,"(A,(1E16.8,x),A)") 'Brillouin zone volume =', c%volume_bz, '1/nm^3'

       if(c%polar) then
          write(*,"(A)") 'System is polar.'
          write(*,"(A)") 'Dielectric tensor:'
          do i = 1, 3
             write(*,"(3(1E16.8,x))") c%epsilon(:,i)
          end do
          write(*,"(A)") 'Born effective charges:'
          do i = 1, c%numatoms
             write(*,"(A, I3)") "  Atom ", i
             do j = 1, 3
                write(*,"(3(1E16.8,x))") c%born(:,j,i)
             end do
          end do
       end if
       write(*,"(A, F7.2, A)") 'Crystal temperature = ', c%T, ' K'
    end if

    !Below is from ShengBTE:

    if(autoisotopes) then

    end if
  end subroutine read_input_and_setup_crystal

  subroutine calculate_wavevectors_full(mesh, wavevecs, blocks, indexlist)
    !! Calculate wave vectors (crystal coords.) of the full Brillouin zone (FBZ)
    !!
    !! mesh is the array of number of points along the reciprocal lattice vectors
    !! wavevecs is the list of all the wave vectors

    integer(k8), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(k8), optional, intent(in) :: indexlist(:)
    real(dp), allocatable, intent(out) :: wavevecs(:,:)
    integer(k8) :: nwavevecs, ijk(3), i, imux

    if(blocks .and. .not. present(indexlist)) &
         call exit_with_message("If blocks is true then indexlist must be present")

    if(blocks) then
       nwavevecs = size(indexlist)
    else
       nwavevecs = product(mesh)
    end if

    allocate(wavevecs(nwavevecs, 3))
    do i = 1, nwavevecs !run over total number of vectors
       if(blocks) then
          imux = indexlist(i)
       else
          imux = i
       end if
       call demux_vector(imux, ijk, mesh, 0_k8) !get 0-based (i,j,k) indices
       wavevecs(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine calculate_wavevectors_full

  subroutine calculate_mavg_and_g(elements, m, g)
    !! Compute the average mass of each element and its g-factor (Pearson
    !! deviation coefficient of the masses).
    !!
    !! This subroutine is adapted from ShengBTE.
    
    character(len=3), intent(in) :: elements(:)
    real(dp), intent(out) :: m(:), g(:)

    !Local variables
    integer(k8) :: i, niso, nelems, e
    character(len = 3) :: isotope_element(287)
    real(dp) :: isotope_mass(287)
    real(dp) :: isotope_abundance(287)

    nelems = size(elements)
    niso = 287
    
    ! Fill in isotope data.
    isotope_element(1) = "Ag"
    isotope_mass(1) = 106.905095_dp
    isotope_abundance(1) = 51.84_dp
    isotope_element(2) = "Ag"
    isotope_mass(2) = 108.904754_dp
    isotope_abundance(2) = 48.16_dp
    isotope_element(3) = "Al"
    isotope_mass(3) = 26.981541_dp
    isotope_abundance(3) = 100.0_dp
    isotope_element(4) = "Ar"
    isotope_mass(4) = 35.967546_dp
    isotope_abundance(4) = 0.34_dp
    isotope_element(5) = "Ar"
    isotope_mass(5) = 37.962732_dp
    isotope_abundance(5) = 0.063_dp
    isotope_element(6) = "Ar"
    isotope_mass(6) = 39.962383_dp
    isotope_abundance(6) = 99.6_dp
    isotope_element(7) = "As"
    isotope_mass(7) = 74.921596_dp
    isotope_abundance(7) = 100.0_dp
    isotope_element(8) = "Au"
    isotope_mass(8) = 196.96656_dp
    isotope_abundance(8) = 100.0_dp
    isotope_element(9) = "B"
    isotope_mass(9) = 10.012938_dp
    isotope_abundance(9) = 19.8_dp
    isotope_element(10) = "B"
    isotope_mass(10) = 11.009305_dp
    isotope_abundance(10) = 80.2_dp
    isotope_element(11) = "Ba"
    isotope_mass(11) = 129.906277_dp
    isotope_abundance(11) = 0.11_dp
    isotope_element(12) = "Ba"
    isotope_mass(12) = 131.905042_dp
    isotope_abundance(12) = 0.1_dp
    isotope_element(13) = "Ba"
    isotope_mass(13) = 133.90449_dp
    isotope_abundance(13) = 2.42_dp
    isotope_element(14) = "Ba"
    isotope_mass(14) = 134.905668_dp
    isotope_abundance(14) = 6.59_dp
    isotope_element(15) = "Ba"
    isotope_mass(15) = 135.904556_dp
    isotope_abundance(15) = 7.85_dp
    isotope_element(16) = "Ba"
    isotope_mass(16) = 136.905816_dp
    isotope_abundance(16) = 11.23_dp
    isotope_element(17) = "Ba"
    isotope_mass(17) = 137.905236_dp
    isotope_abundance(17) = 71.7_dp
    isotope_element(18) = "Be"
    isotope_mass(18) = 9.012183_dp
    isotope_abundance(18) = 100.0_dp
    isotope_element(19) = "Bi"
    isotope_mass(19) = 208.980388_dp
    isotope_abundance(19) = 100.0_dp
    isotope_element(20) = "Br"
    isotope_mass(20) = 78.918336_dp
    isotope_abundance(20) = 50.69_dp
    isotope_element(21) = "Br"
    isotope_mass(21) = 80.91629_dp
    isotope_abundance(21) = 49.31_dp
    isotope_element(22) = "C"
    isotope_mass(22) = 12.0_dp
    isotope_abundance(22) = 98.9_dp
    isotope_element(23) = "C"
    isotope_mass(23) = 13.003355_dp
    isotope_abundance(23) = 1.1_dp
    isotope_element(24) = "Ca"
    isotope_mass(24) = 39.962591_dp
    isotope_abundance(24) = 96.95_dp
    isotope_element(25) = "Ca"
    isotope_mass(25) = 41.958622_dp
    isotope_abundance(25) = 0.65_dp
    isotope_element(26) = "Ca"
    isotope_mass(26) = 42.95877_dp
    isotope_abundance(26) = 0.14_dp
    isotope_element(27) = "Ca"
    isotope_mass(27) = 43.955485_dp
    isotope_abundance(27) = 2.086_dp
    isotope_element(28) = "Ca"
    isotope_mass(28) = 45.953689_dp
    isotope_abundance(28) = 0.004_dp
    isotope_element(29) = "Ca"
    isotope_mass(29) = 47.952532_dp
    isotope_abundance(29) = 0.19_dp
    isotope_element(30) = "Cd"
    isotope_mass(30) = 105.906461_dp
    isotope_abundance(30) = 1.25_dp
    isotope_element(31) = "Cd"
    isotope_mass(31) = 107.904186_dp
    isotope_abundance(31) = 0.89_dp
    isotope_element(32) = "Cd"
    isotope_mass(32) = 109.903007_dp
    isotope_abundance(32) = 12.49_dp
    isotope_element(33) = "Cd"
    isotope_mass(33) = 110.904182_dp
    isotope_abundance(33) = 12.8_dp
    isotope_element(34) = "Cd"
    isotope_mass(34) = 111.902761_dp
    isotope_abundance(34) = 24.13_dp
    isotope_element(35) = "Cd"
    isotope_mass(35) = 112.904401_dp
    isotope_abundance(35) = 12.22_dp
    isotope_element(36) = "Cd"
    isotope_mass(36) = 113.903361_dp
    isotope_abundance(36) = 28.73_dp
    isotope_element(37) = "Cd"
    isotope_mass(37) = 115.904758_dp
    isotope_abundance(37) = 7.49_dp
    isotope_element(38) = "Ce"
    isotope_mass(38) = 135.90714_dp
    isotope_abundance(38) = 0.19_dp
    isotope_element(39) = "Ce"
    isotope_mass(39) = 137.905996_dp
    isotope_abundance(39) = 0.25_dp
    isotope_element(40) = "Ce"
    isotope_mass(40) = 139.905442_dp
    isotope_abundance(40) = 88.48_dp
    isotope_element(41) = "Ce"
    isotope_mass(41) = 141.909249_dp
    isotope_abundance(41) = 11.08_dp
    isotope_element(42) = "Cl"
    isotope_mass(42) = 34.968853_dp
    isotope_abundance(42) = 75.77_dp
    isotope_element(43) = "Cl"
    isotope_mass(43) = 36.965903_dp
    isotope_abundance(43) = 24.23_dp
    isotope_element(44) = "Co"
    isotope_mass(44) = 58.933198_dp
    isotope_abundance(44) = 100.0_dp
    isotope_element(45) = "Cr"
    isotope_mass(45) = 49.946046_dp
    isotope_abundance(45) = 4.35_dp
    isotope_element(46) = "Cr"
    isotope_mass(46) = 51.94051_dp
    isotope_abundance(46) = 83.79_dp
    isotope_element(47) = "Cr"
    isotope_mass(47) = 52.940651_dp
    isotope_abundance(47) = 9.5_dp
    isotope_element(48) = "Cr"
    isotope_mass(48) = 53.938882_dp
    isotope_abundance(48) = 2.36_dp
    isotope_element(49) = "Cs"
    isotope_mass(49) = 132.905433_dp
    isotope_abundance(49) = 100.0_dp
    isotope_element(50) = "Cu"
    isotope_mass(50) = 62.929599_dp
    isotope_abundance(50) = 69.17_dp
    isotope_element(51) = "Cu"
    isotope_mass(51) = 64.927792_dp
    isotope_abundance(51) = 30.83_dp
    isotope_element(52) = "Dy"
    isotope_mass(52) = 155.924287_dp
    isotope_abundance(52) = 0.06_dp
    isotope_element(53) = "Dy"
    isotope_mass(53) = 157.924412_dp
    isotope_abundance(53) = 0.1_dp
    isotope_element(54) = "Dy"
    isotope_mass(54) = 159.925203_dp
    isotope_abundance(54) = 2.34_dp
    isotope_element(55) = "Dy"
    isotope_mass(55) = 160.926939_dp
    isotope_abundance(55) = 18.9_dp
    isotope_element(56) = "Dy"
    isotope_mass(56) = 161.926805_dp
    isotope_abundance(56) = 25.5_dp
    isotope_element(57) = "Dy"
    isotope_mass(57) = 162.928737_dp
    isotope_abundance(57) = 24.9_dp
    isotope_element(58) = "Dy"
    isotope_mass(58) = 163.929183_dp
    isotope_abundance(58) = 28.2_dp
    isotope_element(59) = "Er"
    isotope_mass(59) = 161.928787_dp
    isotope_abundance(59) = 0.14_dp
    isotope_element(60) = "Er"
    isotope_mass(60) = 163.929211_dp
    isotope_abundance(60) = 1.61_dp
    isotope_element(61) = "Er"
    isotope_mass(61) = 165.930305_dp
    isotope_abundance(61) = 33.6_dp
    isotope_element(62) = "Er"
    isotope_mass(62) = 166.932061_dp
    isotope_abundance(62) = 22.95_dp
    isotope_element(63) = "Er"
    isotope_mass(63) = 167.932383_dp
    isotope_abundance(63) = 26.8_dp
    isotope_element(64) = "Er"
    isotope_mass(64) = 169.935476_dp
    isotope_abundance(64) = 14.9_dp
    isotope_element(65) = "Eu"
    isotope_mass(65) = 150.91986_dp
    isotope_abundance(65) = 47.8_dp
    isotope_element(66) = "Eu"
    isotope_mass(66) = 152.921243_dp
    isotope_abundance(66) = 52.2_dp
    isotope_element(67) = "F"
    isotope_mass(67) = 18.998403_dp
    isotope_abundance(67) = 100.0_dp
    isotope_element(68) = "Fe"
    isotope_mass(68) = 53.939612_dp
    isotope_abundance(68) = 5.8_dp
    isotope_element(69) = "Fe"
    isotope_mass(69) = 55.934939_dp
    isotope_abundance(69) = 91.72_dp
    isotope_element(70) = "Fe"
    isotope_mass(70) = 56.935396_dp
    isotope_abundance(70) = 2.2_dp
    isotope_element(71) = "Fe"
    isotope_mass(71) = 57.933278_dp
    isotope_abundance(71) = 0.28_dp
    isotope_element(72) = "Ga"
    isotope_mass(72) = 68.925581_dp
    isotope_abundance(72) = 60.1_dp
    isotope_element(73) = "Ga"
    isotope_mass(73) = 70.924701_dp
    isotope_abundance(73) = 39.9_dp
    isotope_element(74) = "Gd"
    isotope_mass(74) = 151.919803_dp
    isotope_abundance(74) = 0.2_dp
    isotope_element(75) = "Gd"
    isotope_mass(75) = 153.920876_dp
    isotope_abundance(75) = 2.18_dp
    isotope_element(76) = "Gd"
    isotope_mass(76) = 154.822629_dp
    isotope_abundance(76) = 14.8_dp
    isotope_element(77) = "Gd"
    isotope_mass(77) = 155.92213_dp
    isotope_abundance(77) = 20.47_dp
    isotope_element(78) = "Gd"
    isotope_mass(78) = 156.923967_dp
    isotope_abundance(78) = 15.65_dp
    isotope_element(79) = "Gd"
    isotope_mass(79) = 157.924111_dp
    isotope_abundance(79) = 24.84_dp
    isotope_element(80) = "Gd"
    isotope_mass(80) = 159.927061_dp
    isotope_abundance(80) = 21.86_dp
    isotope_element(81) = "Ge"
    isotope_mass(81) = 69.92425_dp
    isotope_abundance(81) = 20.5_dp
    isotope_element(82) = "Ge"
    isotope_mass(82) = 71.92208_dp
    isotope_abundance(82) = 27.4_dp
    isotope_element(83) = "Ge"
    isotope_mass(83) = 72.923464_dp
    isotope_abundance(83) = 7.8_dp
    isotope_element(84) = "Ge"
    isotope_mass(84) = 73.921179_dp
    isotope_abundance(84) = 36.5_dp
    isotope_element(85) = "Ge"
    isotope_mass(85) = 75.921403_dp
    isotope_abundance(85) = 7.8_dp
    isotope_element(86) = "H"
    isotope_mass(86) = 1.007825_dp
    isotope_abundance(86) = 99.99_dp
    isotope_element(87) = "H"
    isotope_mass(87) = 2.014102_dp
    isotope_abundance(87) = 0.015_dp
    isotope_element(88) = "He"
    isotope_mass(88) = 3.016029_dp
    isotope_abundance(88) = 0.0001_dp
    isotope_element(89) = "He"
    isotope_mass(89) = 4.002603_dp
    isotope_abundance(89) = 100.0_dp
    isotope_element(90) = "Hf"
    isotope_mass(90) = 173.940065_dp
    isotope_abundance(90) = 0.16_dp
    isotope_element(91) = "Hf"
    isotope_mass(91) = 175.94142_dp
    isotope_abundance(91) = 5.2_dp
    isotope_element(92) = "Hf"
    isotope_mass(92) = 176.943233_dp
    isotope_abundance(92) = 18.6_dp
    isotope_element(93) = "Hf"
    isotope_mass(93) = 177.94371_dp
    isotope_abundance(93) = 27.1_dp
    isotope_element(94) = "Hf"
    isotope_mass(94) = 178.945827_dp
    isotope_abundance(94) = 13.74_dp
    isotope_element(95) = "Hf"
    isotope_mass(95) = 179.946561_dp
    isotope_abundance(95) = 35.2_dp
    isotope_element(96) = "Hg"
    isotope_mass(96) = 195.965812_dp
    isotope_abundance(96) = 0.15_dp
    isotope_element(97) = "Hg"
    isotope_mass(97) = 197.96676_dp
    isotope_abundance(97) = 10.1_dp
    isotope_element(98) = "Hg"
    isotope_mass(98) = 198.968269_dp
    isotope_abundance(98) = 17.0_dp
    isotope_element(99) = "Hg"
    isotope_mass(99) = 199.968316_dp
    isotope_abundance(99) = 23.1_dp
    isotope_element(100) = "Hg"
    isotope_mass(100) = 200.970293_dp
    isotope_abundance(100) = 13.2_dp
    isotope_element(101) = "Hg"
    isotope_mass(101) = 201.970632_dp
    isotope_abundance(101) = 29.65_dp
    isotope_element(102) = "Hg"
    isotope_mass(102) = 203.973481_dp
    isotope_abundance(102) = 6.8_dp
    isotope_element(103) = "Ho"
    isotope_mass(103) = 164.930332_dp
    isotope_abundance(103) = 100.0_dp
    isotope_element(104) = "I"
    isotope_mass(104) = 126.904477_dp
    isotope_abundance(104) = 100.0_dp
    isotope_element(105) = "In"
    isotope_mass(105) = 112.904056_dp
    isotope_abundance(105) = 4.3_dp
    isotope_element(106) = "In"
    isotope_mass(106) = 114.903875_dp
    isotope_abundance(106) = 95.7_dp
    isotope_element(107) = "Ir"
    isotope_mass(107) = 190.960603_dp
    isotope_abundance(107) = 37.3_dp
    isotope_element(108) = "Ir"
    isotope_mass(108) = 192.962942_dp
    isotope_abundance(108) = 62.7_dp
    isotope_element(109) = "K"
    isotope_mass(109) = 38.963708_dp
    isotope_abundance(109) = 93.2_dp
    isotope_element(110) = "K"
    isotope_mass(110) = 39.963999_dp
    isotope_abundance(110) = 0.012_dp
    isotope_element(111) = "K"
    isotope_mass(111) = 40.961825_dp
    isotope_abundance(111) = 6.73_dp
    isotope_element(112) = "Kr"
    isotope_mass(112) = 77.920397_dp
    isotope_abundance(112) = 0.35_dp
    isotope_element(113) = "Kr"
    isotope_mass(113) = 79.916375_dp
    isotope_abundance(113) = 2.25_dp
    isotope_element(114) = "Kr"
    isotope_mass(114) = 81.913483_dp
    isotope_abundance(114) = 11.6_dp
    isotope_element(115) = "Kr"
    isotope_mass(115) = 82.914134_dp
    isotope_abundance(115) = 11.5_dp
    isotope_element(116) = "Kr"
    isotope_mass(116) = 83.911506_dp
    isotope_abundance(116) = 57.0_dp
    isotope_element(117) = "Kr"
    isotope_mass(117) = 85.910614_dp
    isotope_abundance(117) = 17.3_dp
    isotope_element(118) = "La"
    isotope_mass(118) = 137.907114_dp
    isotope_abundance(118) = 0.09_dp
    isotope_element(119) = "La"
    isotope_mass(119) = 138.906355_dp
    isotope_abundance(119) = 99.91_dp
    isotope_element(120) = "Li"
    isotope_mass(120) = 6.015123_dp
    isotope_abundance(120) = 7.42_dp
    isotope_element(121) = "Li"
    isotope_mass(121) = 7.016005_dp
    isotope_abundance(121) = 92.58_dp
    isotope_element(122) = "Lu"
    isotope_mass(122) = 174.940785_dp
    isotope_abundance(122) = 97.4_dp
    isotope_element(123) = "Lu"
    isotope_mass(123) = 175.942694_dp
    isotope_abundance(123) = 2.6_dp
    isotope_element(124) = "Mg"
    isotope_mass(124) = 23.985045_dp
    isotope_abundance(124) = 78.9_dp
    isotope_element(125) = "Mg"
    isotope_mass(125) = 24.985839_dp
    isotope_abundance(125) = 10.0_dp
    isotope_element(126) = "Mg"
    isotope_mass(126) = 25.982595_dp
    isotope_abundance(126) = 11.1_dp
    isotope_element(127) = "Mn"
    isotope_mass(127) = 54.938046_dp
    isotope_abundance(127) = 100.0_dp
    isotope_element(128) = "Mo"
    isotope_mass(128) = 91.906809_dp
    isotope_abundance(128) = 14.84_dp
    isotope_element(129) = "Mo"
    isotope_mass(129) = 93.905086_dp
    isotope_abundance(129) = 9.25_dp
    isotope_element(130) = "Mo"
    isotope_mass(130) = 94.905838_dp
    isotope_abundance(130) = 15.92_dp
    isotope_element(131) = "Mo"
    isotope_mass(131) = 95.904676_dp
    isotope_abundance(131) = 16.68_dp
    isotope_element(132) = "Mo"
    isotope_mass(132) = 96.906018_dp
    isotope_abundance(132) = 9.55_dp
    isotope_element(133) = "Mo"
    isotope_mass(133) = 97.905405_dp
    isotope_abundance(133) = 24.13_dp
    isotope_element(134) = "Mo"
    isotope_mass(134) = 99.907473_dp
    isotope_abundance(134) = 9.63_dp
    isotope_element(135) = "N"
    isotope_mass(135) = 14.003074_dp
    isotope_abundance(135) = 99.63_dp
    isotope_element(136) = "N"
    isotope_mass(136) = 15.000109_dp
    isotope_abundance(136) = 0.37_dp
    isotope_element(137) = "Na"
    isotope_mass(137) = 22.98977_dp
    isotope_abundance(137) = 100.0_dp
    isotope_element(138) = "Nb"
    isotope_mass(138) = 92.906378_dp
    isotope_abundance(138) = 100.0_dp
    isotope_element(139) = "Nd"
    isotope_mass(139) = 141.907731_dp
    isotope_abundance(139) = 27.13_dp
    isotope_element(140) = "Nd"
    isotope_mass(140) = 142.909823_dp
    isotope_abundance(140) = 12.18_dp
    isotope_element(141) = "Nd"
    isotope_mass(141) = 143.910096_dp
    isotope_abundance(141) = 23.8_dp
    isotope_element(142) = "Nd"
    isotope_mass(142) = 144.912582_dp
    isotope_abundance(142) = 8.3_dp
    isotope_element(143) = "Nd"
    isotope_mass(143) = 145.913126_dp
    isotope_abundance(143) = 17.19_dp
    isotope_element(144) = "Nd"
    isotope_mass(144) = 147.916901_dp
    isotope_abundance(144) = 5.76_dp
    isotope_element(145) = "Nd"
    isotope_mass(145) = 149.9209_dp
    isotope_abundance(145) = 5.64_dp
    isotope_element(146) = "Ne"
    isotope_mass(146) = 19.992439_dp
    isotope_abundance(146) = 90.6_dp
    isotope_element(147) = "Ne"
    isotope_mass(147) = 20.993845_dp
    isotope_abundance(147) = 0.26_dp
    isotope_element(148) = "Ne"
    isotope_mass(148) = 21.991384_dp
    isotope_abundance(148) = 9.2_dp
    isotope_element(149) = "Ni"
    isotope_mass(149) = 57.935347_dp
    isotope_abundance(149) = 68.27_dp
    isotope_element(150) = "Ni"
    isotope_mass(150) = 59.930789_dp
    isotope_abundance(150) = 26.1_dp
    isotope_element(151) = "Ni"
    isotope_mass(151) = 60.931059_dp
    isotope_abundance(151) = 1.13_dp
    isotope_element(152) = "Ni"
    isotope_mass(152) = 61.928346_dp
    isotope_abundance(152) = 3.59_dp
    isotope_element(153) = "Ni"
    isotope_mass(153) = 63.927968_dp
    isotope_abundance(153) = 0.91_dp
    isotope_element(154) = "O"
    isotope_mass(154) = 15.994915_dp
    isotope_abundance(154) = 99.76_dp
    isotope_element(155) = "O"
    isotope_mass(155) = 16.999131_dp
    isotope_abundance(155) = 0.038_dp
    isotope_element(156) = "O"
    isotope_mass(156) = 17.999159_dp
    isotope_abundance(156) = 0.2_dp
    isotope_element(157) = "Os"
    isotope_mass(157) = 183.952514_dp
    isotope_abundance(157) = 0.02_dp
    isotope_element(158) = "Os"
    isotope_mass(158) = 185.953852_dp
    isotope_abundance(158) = 1.58_dp
    isotope_element(159) = "Os"
    isotope_mass(159) = 186.955762_dp
    isotope_abundance(159) = 1.6_dp
    isotope_element(160) = "Os"
    isotope_mass(160) = 187.95585_dp
    isotope_abundance(160) = 13.3_dp
    isotope_element(161) = "Os"
    isotope_mass(161) = 188.958156_dp
    isotope_abundance(161) = 16.1_dp
    isotope_element(162) = "Os"
    isotope_mass(162) = 189.958455_dp
    isotope_abundance(162) = 26.4_dp
    isotope_element(163) = "Os"
    isotope_mass(163) = 191.961487_dp
    isotope_abundance(163) = 41.0_dp
    isotope_element(164) = "P"
    isotope_mass(164) = 30.973763_dp
    isotope_abundance(164) = 100.0_dp
    isotope_element(165) = "Pb"
    isotope_mass(165) = 203.973037_dp
    isotope_abundance(165) = 1.4_dp
    isotope_element(166) = "Pb"
    isotope_mass(166) = 205.974455_dp
    isotope_abundance(166) = 24.1_dp
    isotope_element(167) = "Pb"
    isotope_mass(167) = 206.975885_dp
    isotope_abundance(167) = 22.1_dp
    isotope_element(168) = "Pb"
    isotope_mass(168) = 207.976641_dp
    isotope_abundance(168) = 52.4_dp
    isotope_element(169) = "Pd"
    isotope_mass(169) = 101.905609_dp
    isotope_abundance(169) = 1.02_dp
    isotope_element(170) = "Pd"
    isotope_mass(170) = 103.904026_dp
    isotope_abundance(170) = 11.14_dp
    isotope_element(171) = "Pd"
    isotope_mass(171) = 104.905075_dp
    isotope_abundance(171) = 22.33_dp
    isotope_element(172) = "Pd"
    isotope_mass(172) = 105.903475_dp
    isotope_abundance(172) = 27.33_dp
    isotope_element(173) = "Pd"
    isotope_mass(173) = 107.903894_dp
    isotope_abundance(173) = 26.46_dp
    isotope_element(174) = "Pd"
    isotope_mass(174) = 109.905169_dp
    isotope_abundance(174) = 11.72_dp
    isotope_element(175) = "Pr"
    isotope_mass(175) = 140.907657_dp
    isotope_abundance(175) = 100.0_dp
    isotope_element(176) = "Pt"
    isotope_mass(176) = 189.959937_dp
    isotope_abundance(176) = 0.01_dp
    isotope_element(177) = "Pt"
    isotope_mass(177) = 191.961049_dp
    isotope_abundance(177) = 0.79_dp
    isotope_element(178) = "Pt"
    isotope_mass(178) = 193.962679_dp
    isotope_abundance(178) = 32.9_dp
    isotope_element(179) = "Pt"
    isotope_mass(179) = 194.964785_dp
    isotope_abundance(179) = 33.8_dp
    isotope_element(180) = "Pt"
    isotope_mass(180) = 195.964947_dp
    isotope_abundance(180) = 25.3_dp
    isotope_element(181) = "Pt"
    isotope_mass(181) = 197.967879_dp
    isotope_abundance(181) = 7.2_dp
    isotope_element(182) = "Rb"
    isotope_mass(182) = 84.9118_dp
    isotope_abundance(182) = 72.17_dp
    isotope_element(183) = "Rb"
    isotope_mass(183) = 86.909184_dp
    isotope_abundance(183) = 27.84_dp
    isotope_element(184) = "Re"
    isotope_mass(184) = 184.952977_dp
    isotope_abundance(184) = 37.4_dp
    isotope_element(185) = "Re"
    isotope_mass(185) = 186.955765_dp
    isotope_abundance(185) = 62.6_dp
    isotope_element(186) = "Rh"
    isotope_mass(186) = 102.905503_dp
    isotope_abundance(186) = 100.0_dp
    isotope_element(187) = "Ru"
    isotope_mass(187) = 95.907596_dp
    isotope_abundance(187) = 5.52_dp
    isotope_element(188) = "Ru"
    isotope_mass(188) = 97.905287_dp
    isotope_abundance(188) = 1.88_dp
    isotope_element(189) = "Ru"
    isotope_mass(189) = 98.905937_dp
    isotope_abundance(189) = 12.7_dp
    isotope_element(190) = "Ru"
    isotope_mass(190) = 99.904218_dp
    isotope_abundance(190) = 12.6_dp
    isotope_element(191) = "Ru"
    isotope_mass(191) = 100.905581_dp
    isotope_abundance(191) = 17.0_dp
    isotope_element(192) = "Ru"
    isotope_mass(192) = 101.904348_dp
    isotope_abundance(192) = 31.6_dp
    isotope_element(193) = "Ru"
    isotope_mass(193) = 103.905422_dp
    isotope_abundance(193) = 18.7_dp
    isotope_element(194) = "S"
    isotope_mass(194) = 31.972072_dp
    isotope_abundance(194) = 95.02_dp
    isotope_element(195) = "S"
    isotope_mass(195) = 32.971459_dp
    isotope_abundance(195) = 0.75_dp
    isotope_element(196) = "S"
    isotope_mass(196) = 33.967868_dp
    isotope_abundance(196) = 4.21_dp
    isotope_element(197) = "S"
    isotope_mass(197) = 35.967079_dp
    isotope_abundance(197) = 0.02_dp
    isotope_element(198) = "Sb"
    isotope_mass(198) = 120.903824_dp
    isotope_abundance(198) = 57.3_dp
    isotope_element(199) = "Sb"
    isotope_mass(199) = 122.904222_dp
    isotope_abundance(199) = 42.7_dp
    isotope_element(200) = "Sc"
    isotope_mass(200) = 44.955914_dp
    isotope_abundance(200) = 100.0_dp
    isotope_element(201) = "Se"
    isotope_mass(201) = 73.922477_dp
    isotope_abundance(201) = 0.9_dp
    isotope_element(202) = "Se"
    isotope_mass(202) = 75.919207_dp
    isotope_abundance(202) = 9.0_dp
    isotope_element(203) = "Se"
    isotope_mass(203) = 76.919908_dp
    isotope_abundance(203) = 7.6_dp
    isotope_element(204) = "Se"
    isotope_mass(204) = 77.917304_dp
    isotope_abundance(204) = 23.5_dp
    isotope_element(205) = "Se"
    isotope_mass(205) = 79.916521_dp
    isotope_abundance(205) = 49.6_dp
    isotope_element(206) = "Se"
    isotope_mass(206) = 81.916709_dp
    isotope_abundance(206) = 9.4_dp
    isotope_element(207) = "Si"
    isotope_mass(207) = 27.976928_dp
    isotope_abundance(207) = 92.23_dp
    isotope_element(208) = "Si"
    isotope_mass(208) = 28.976496_dp
    isotope_abundance(208) = 4.67_dp
    isotope_element(209) = "Si"
    isotope_mass(209) = 29.973772_dp
    isotope_abundance(209) = 3.1_dp
    isotope_element(210) = "Sm"
    isotope_mass(210) = 143.912009_dp
    isotope_abundance(210) = 3.1_dp
    isotope_element(211) = "Sm"
    isotope_mass(211) = 146.914907_dp
    isotope_abundance(211) = 15.0_dp
    isotope_element(212) = "Sm"
    isotope_mass(212) = 147.914832_dp
    isotope_abundance(212) = 11.3_dp
    isotope_element(213) = "Sm"
    isotope_mass(213) = 148.917193_dp
    isotope_abundance(213) = 13.8_dp
    isotope_element(214) = "Sm"
    isotope_mass(214) = 149.917285_dp
    isotope_abundance(214) = 7.4_dp
    isotope_element(215) = "Sm"
    isotope_mass(215) = 151.919741_dp
    isotope_abundance(215) = 26.7_dp
    isotope_element(216) = "Sm"
    isotope_mass(216) = 153.922218_dp
    isotope_abundance(216) = 22.7_dp
    isotope_element(217) = "Sn"
    isotope_mass(217) = 111.904826_dp
    isotope_abundance(217) = 0.97_dp
    isotope_element(218) = "Sn"
    isotope_mass(218) = 113.902784_dp
    isotope_abundance(218) = 0.65_dp
    isotope_element(219) = "Sn"
    isotope_mass(219) = 114.903348_dp
    isotope_abundance(219) = 0.36_dp
    isotope_element(220) = "Sn"
    isotope_mass(220) = 115.901744_dp
    isotope_abundance(220) = 14.7_dp
    isotope_element(221) = "Sn"
    isotope_mass(221) = 116.902954_dp
    isotope_abundance(221) = 7.7_dp
    isotope_element(222) = "Sn"
    isotope_mass(222) = 117.901607_dp
    isotope_abundance(222) = 24.3_dp
    isotope_element(223) = "Sn"
    isotope_mass(223) = 118.90331_dp
    isotope_abundance(223) = 8.6_dp
    isotope_element(224) = "Sn"
    isotope_mass(224) = 119.902199_dp
    isotope_abundance(224) = 32.4_dp
    isotope_element(225) = "Sn"
    isotope_mass(225) = 121.90344_dp
    isotope_abundance(225) = 4.6_dp
    isotope_element(226) = "Sn"
    isotope_mass(226) = 123.905271_dp
    isotope_abundance(226) = 5.6_dp
    isotope_element(227) = "Sr"
    isotope_mass(227) = 83.913428_dp
    isotope_abundance(227) = 0.56_dp
    isotope_element(228) = "Sr"
    isotope_mass(228) = 85.909273_dp
    isotope_abundance(228) = 9.86_dp
    isotope_element(229) = "Sr"
    isotope_mass(229) = 86.908902_dp
    isotope_abundance(229) = 7.0_dp
    isotope_element(230) = "Sr"
    isotope_mass(230) = 87.905625_dp
    isotope_abundance(230) = 82.58_dp
    isotope_element(231) = "Ta"
    isotope_mass(231) = 179.947489_dp
    isotope_abundance(231) = 0.012_dp
    isotope_element(232) = "Ta"
    isotope_mass(232) = 180.948014_dp
    isotope_abundance(232) = 99.99_dp
    isotope_element(233) = "Tb"
    isotope_mass(233) = 158.92535_dp
    isotope_abundance(233) = 100.0_dp
    isotope_element(234) = "Te"
    isotope_mass(234) = 119.904021_dp
    isotope_abundance(234) = 0.096_dp
    isotope_element(235) = "Te"
    isotope_mass(235) = 121.903055_dp
    isotope_abundance(235) = 2.6_dp
    isotope_element(236) = "Te"
    isotope_mass(236) = 122.904278_dp
    isotope_abundance(236) = 0.91_dp
    isotope_element(237) = "Te"
    isotope_mass(237) = 123.902825_dp
    isotope_abundance(237) = 4.82_dp
    isotope_element(238) = "Te"
    isotope_mass(238) = 124.904435_dp
    isotope_abundance(238) = 7.14_dp
    isotope_element(239) = "Te"
    isotope_mass(239) = 125.90331_dp
    isotope_abundance(239) = 18.95_dp
    isotope_element(240) = "Te"
    isotope_mass(240) = 127.904464_dp
    isotope_abundance(240) = 31.69_dp
    isotope_element(241) = "Te"
    isotope_mass(241) = 129.906229_dp
    isotope_abundance(241) = 33.8_dp
    isotope_element(242) = "Th"
    isotope_mass(242) = 232.038054_dp
    isotope_abundance(242) = 100.0_dp
    isotope_element(243) = "Ti"
    isotope_mass(243) = 45.952633_dp
    isotope_abundance(243) = 8.0_dp
    isotope_element(244) = "Ti"
    isotope_mass(244) = 46.951765_dp
    isotope_abundance(244) = 7.3_dp
    isotope_element(245) = "Ti"
    isotope_mass(245) = 47.947947_dp
    isotope_abundance(245) = 73.8_dp
    isotope_element(246) = "Ti"
    isotope_mass(246) = 48.947871_dp
    isotope_abundance(246) = 5.5_dp
    isotope_element(247) = "Ti"
    isotope_mass(247) = 49.944786_dp
    isotope_abundance(247) = 5.4_dp
    isotope_element(248) = "Tl"
    isotope_mass(248) = 202.972336_dp
    isotope_abundance(248) = 29.52_dp
    isotope_element(249) = "Tl"
    isotope_mass(249) = 204.97441_dp
    isotope_abundance(249) = 70.48_dp
    isotope_element(250) = "Tm"
    isotope_mass(250) = 168.934225_dp
    isotope_abundance(250) = 100.0_dp
    isotope_element(251) = "U"
    isotope_mass(251) = 234.040947_dp
    isotope_abundance(251) = 0.006_dp
    isotope_element(252) = "U"
    isotope_mass(252) = 235.043925_dp
    isotope_abundance(252) = 0.72_dp
    isotope_element(253) = "U"
    isotope_mass(253) = 238.050786_dp
    isotope_abundance(253) = 99.27_dp
    isotope_element(254) = "V"
    isotope_mass(254) = 49.947161_dp
    isotope_abundance(254) = 0.25_dp
    isotope_element(255) = "V"
    isotope_mass(255) = 50.943963_dp
    isotope_abundance(255) = 99.75_dp
    isotope_element(256) = "W"
    isotope_mass(256) = 179.946727_dp
    isotope_abundance(256) = 0.13_dp
    isotope_element(257) = "W"
    isotope_mass(257) = 181.948225_dp
    isotope_abundance(257) = 26.3_dp
    isotope_element(258) = "W"
    isotope_mass(258) = 182.950245_dp
    isotope_abundance(258) = 14.3_dp
    isotope_element(259) = "W"
    isotope_mass(259) = 183.950953_dp
    isotope_abundance(259) = 30.67_dp
    isotope_element(260) = "W"
    isotope_mass(260) = 185.954377_dp
    isotope_abundance(260) = 28.6_dp
    isotope_element(261) = "Xe"
    isotope_mass(261) = 123.905894_dp
    isotope_abundance(261) = 0.1_dp
    isotope_element(262) = "Xe"
    isotope_mass(262) = 125.904281_dp
    isotope_abundance(262) = 0.09_dp
    isotope_element(263) = "Xe"
    isotope_mass(263) = 127.903531_dp
    isotope_abundance(263) = 1.91_dp
    isotope_element(264) = "Xe"
    isotope_mass(264) = 128.90478_dp
    isotope_abundance(264) = 26.4_dp
    isotope_element(265) = "Xe"
    isotope_mass(265) = 129.90351_dp
    isotope_abundance(265) = 4.1_dp
    isotope_element(266) = "Xe"
    isotope_mass(266) = 130.905076_dp
    isotope_abundance(266) = 21.2_dp
    isotope_element(267) = "Xe"
    isotope_mass(267) = 131.904148_dp
    isotope_abundance(267) = 26.9_dp
    isotope_element(268) = "Xe"
    isotope_mass(268) = 133.905395_dp
    isotope_abundance(268) = 10.4_dp
    isotope_element(269) = "Xe"
    isotope_mass(269) = 135.907219_dp
    isotope_abundance(269) = 8.9_dp
    isotope_element(270) = "Y"
    isotope_mass(270) = 88.905856_dp
    isotope_abundance(270) = 100.0_dp
    isotope_element(271) = "Yb"
    isotope_mass(271) = 167.933908_dp
    isotope_abundance(271) = 0.13_dp
    isotope_element(272) = "Yb"
    isotope_mass(272) = 169.934774_dp
    isotope_abundance(272) = 3.05_dp
    isotope_element(273) = "Yb"
    isotope_mass(273) = 170.936338_dp
    isotope_abundance(273) = 14.3_dp
    isotope_element(274) = "Yb"
    isotope_mass(274) = 171.936393_dp
    isotope_abundance(274) = 21.9_dp
    isotope_element(275) = "Yb"
    isotope_mass(275) = 172.938222_dp
    isotope_abundance(275) = 16.12_dp
    isotope_element(276) = "Yb"
    isotope_mass(276) = 173.938873_dp
    isotope_abundance(276) = 31.8_dp
    isotope_element(277) = "Yb"
    isotope_mass(277) = 175.942576_dp
    isotope_abundance(277) = 12.7_dp
    isotope_element(278) = "Zn"
    isotope_mass(278) = 63.929145_dp
    isotope_abundance(278) = 48.6_dp
    isotope_element(279) = "Zn"
    isotope_mass(279) = 65.926035_dp
    isotope_abundance(279) = 27.9_dp
    isotope_element(280) = "Zn"
    isotope_mass(280) = 66.927129_dp
    isotope_abundance(280) = 4.1_dp
    isotope_element(281) = "Zn"
    isotope_mass(281) = 67.924846_dp
    isotope_abundance(281) = 18.8_dp
    isotope_element(282) = "Zn"
    isotope_mass(282) = 69.925325_dp
    isotope_abundance(282) = 0.6_dp
    isotope_element(283) = "Zr"
    isotope_mass(283) = 89.904708_dp
    isotope_abundance(283) = 51.45_dp
    isotope_element(284) = "Zr"
    isotope_mass(284) = 90.905644_dp
    isotope_abundance(284) = 11.27_dp
    isotope_element(285) = "Zr"
    isotope_mass(285) = 91.905039_dp
    isotope_abundance(285) = 17.17_dp
    isotope_element(286) = "Zr"
    isotope_mass(286) = 93.906319_dp
    isotope_abundance(286) = 17.33_dp
    isotope_element(287) = "Zr"
    isotope_mass(287) = 95.908272_dp
    isotope_abundance(287) = 2.78_dp
    
    do e = 1, nelems
       m(e) = 0.0_dp
       do i = 1, niso
          if(isotope_element(i) .eq. elements(e)) then
             m(e) = m(e) + isotope_mass(i)*isotope_abundance(i)
          end if
       end do
       m(e) = m(e)/100.0_dp
       g(e) = 0.0_dp
       do i = 1, niso
          if(isotope_element(i) .eq. elements(e)) then
             g(e) = g(e) + isotope_abundance(i)*&
                  (1.0_dp - isotope_mass(i)/m(e))**2
          end if
       end do
       g(e) = g(e)/100.0_dp
    end do
  end subroutine calculate_mavg_and_g
end module crystal_module

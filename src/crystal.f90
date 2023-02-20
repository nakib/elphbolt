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

module crystal_module
  !! Module containing type and procedures related to the crystal structure.

  use params, only: r64, i64, twopi
  use misc, only: exit_with_message, print_message, cross_product, demux_vector, &
       subtitle, trace

  implicit none

  private
  public crystal, calculate_wavevectors_full

  type crystal
     !! Data and procedures related to the crystal structure.

     integer(i64) :: numelements
     !! Number of types of basis atoms.
     integer(i64) :: numatoms
     !! Number of basis atoms.
     character(len=100) :: name
     !! Name of material.
     character(len=3), allocatable :: elements(:)
     !! Elements in the basis.
     integer(i64), allocatable :: atomtypes(:)
     !! Integer tagging unique elements in the basis.
     real(r64), allocatable :: masses(:)
     !! Masses of the basis atoms.
     logical :: polar
     !! Is the system polar?
     real(r64) :: epsilon(3, 3)
     !! Dielectric tensor.
     real(r64), allocatable :: born(:,:,:)
     !! Born effective charge.
     real(r64) :: epsilon0
     !! Static dielectric constant.
     logical :: read_epsiloninf
     !! Read high-frequency dielectric constant?
     real(r64) :: epsiloninf
     !! High frequency dielectric constant.
     real(r64) :: qTF
     !! Thomas-Fermi screening wave vector.
     real(r64), allocatable :: basis(:,:)
     !! Basis vectors (crystal coordinates).
     real(r64), allocatable :: basis_cart(:,:)
     !! Basis vectors (Cartesian coordinates).
     real(r64) :: lattvecs(3, 3)
     !! Lattice vectors (nm).
     real(r64) :: volume
     !! Volume of primitive cell (nm^3).
     real(r64) :: reclattvecs(3,3)
     !! Reciprocal lattice vectors.
     real(r64) :: volume_bz
     !! Brillouin zone volume (nm^-3).
     real(r64) :: T
     !! Crystal temperature (K).
     logical :: autoisotopes
     !! Use isotopic mix for masses?
     real(r64), allocatable :: gfactors(:)
     !! g-factors.
     real(r64), allocatable :: subs_masses(:)
     !! Masses of the substitutional atoms [D]
     real(r64), allocatable :: subs_conc(:)
     !! Concentration of the substitutional atoms in cm^-3 [D]
     real(r64), allocatable :: subs_gfactors(:)
     !! g-factors for the substitutional defects. [D]
     integer(i64), allocatable :: defect_hosts(:)
     !! Basis atom sites that can be a host for an impurity, one for each unique element.
     integer(i64), allocatable :: numdopants_types(:)
     !! Number of dopant types at each host atom site.
     real(r64), allocatable :: dopant_masses(:, :)
     !! Masses of the dopants at each host atom site.
     real(r64), allocatable :: dopant_conc(:, :)
     !! Concentrations [cm^-3] of the dopants at each host atom site.
     logical :: twod
     !! Is the system 2d?
     real(r64) :: dim
     !! Dimension of the system
     real(r64) :: thickness
     !! Thickness of the system
     real(r64) :: bound_length
     !! Characteristic boundary scattering length in mm
     real(r64) :: thinfilm_height
     !! Height of thin-film in mm
     character(1) :: thinfilm_normal
     !! Normal direction of the thin-film: 'x', 'y', or 'z'.
     
   contains

     procedure :: initialize=>read_input_and_setup_crystal
  end type crystal

contains

  subroutine read_input_and_setup_crystal(self)
    !! Read input file and initialize crystal data.

    class(crystal), intent(out) :: self

    !Local variables
    integer(i64) :: i, j, k, numelements, numatoms
    integer(i64), allocatable :: atomtypes(:), num_atomtypes(:), defect_hosts(:), numdopants_types(:)
    real(r64), allocatable :: masses(:), born(:,:,:), basis(:,:), &
         basis_cart(:,:), subs_perc(:), subs_masses(:), subs_conc(:), &
         dopant_masses(:, :), dopant_conc(:, :)
    real(r64) :: epsilon(3,3), lattvecs(3,3), T, &
         epsilon0, epsiloninf, subs_mavg, bound_length, thinfilm_height
    character(len=3), allocatable :: elements(:)
    character(len=100) :: name
    character(1) :: thinfilm_normal
    logical :: polar, autoisotopes, read_epsiloninf, twod
    
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ name, elements, atomtypes, basis, lattvecs, &
         polar, born, epsilon, read_epsiloninf, epsilon0, epsiloninf, &
         masses, T, autoisotopes, twod, subs_masses, subs_conc, bound_length, &
         defect_hosts, numdopants_types, dopant_masses, dopant_conc, thinfilm_height, &
         thinfilm_normal

    call subtitle("Setting up crystal...")

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Set values from input:
    
    !Read allocations
    numelements = 0
    numatoms = 0
    read(1, nml = allocations)
    if(numelements < 1 .or. numatoms < 1 .or. numatoms < numelements) then
       call exit_with_message('Bad input(s) in allocations.')
    end if    
    self%numelements = numelements
    self%numatoms = numatoms
    
    !Allocate variables
    allocate(elements(numelements), atomtypes(numatoms), born(3,3,numatoms), &
         basis(3,numatoms), masses(numelements), basis_cart(3,numatoms), &
         subs_masses(numelements), subs_conc(numelements), subs_perc(numelements), &
         num_atomtypes(numelements), defect_hosts(numelements), &
         numdopants_types(numelements), dopant_masses(10, numelements), dopant_conc(10, numelements))
    allocate(self%elements(self%numelements), self%atomtypes(self%numatoms), self%born(3,3,self%numatoms), &
         self%masses(self%numatoms), self%gfactors(self%numelements), self%basis(3,self%numatoms), &
         self%basis_cart(3,self%numatoms), self%subs_masses(self%numelements), self%subs_conc(self%numelements), &
         self%subs_gfactors(self%numelements), self%defect_hosts(self%numelements))
    
    !Read crystal_info
    name = trim(adjustl('Crystal'))
    elements = 'X'
    atomtypes = 0
    masses = -1.0_r64
    autoisotopes = .true.
    lattvecs = 0.0_r64
    basis = 0.0_r64
    polar = .false.
    read_epsiloninf = .false.
    epsilon = 0.0_r64
    epsilon0 = 0.0_r64
    epsiloninf = 0.0_r64
    born = 0.0_r64
    T = -1.0_r64
    twod = .false.
    subs_masses = 0.0_r64
    subs_conc = 0.0_r64
    defect_hosts = -1_i64
    bound_length = 1.e12_r64 !mm, practically inifinity
    thinfilm_height = 1.e12_r64 !mm, practically inifinity
    thinfilm_normal = 'z'
    numdopants_types = [1, 1]
    dopant_masses = 0.0_r64
    dopant_conc = 0.0_r64
    read(1, nml = crystal_info)
    if(any(atomtypes < 1) .or. T < 0.0_r64) then
       call exit_with_message('Bad input(s) in crystal_info.')
    end if
    if(.not. autoisotopes .and. any(masses < 0)) then
       call exit_with_message('Bad input(s) in crystal_info.')
    end if
    if(bound_length <= 0.0_r64) then
       call exit_with_message('Characteristic length for boundary scattering must be positive.')
    end if
    if(thinfilm_height <= 0.0_r64) then
       call exit_with_message('Height of thin-film must be positive.')
    end if
    if(.not. (thinfilm_normal == 'x' .or.  thinfilm_normal == 'y' .or. thinfilm_normal == 'z')) then
       call exit_with_message("Thin-film normal direction must be 'x', 'y', or 'z'.")
    end if

    !Close input file
    close(1)
    
    self%name = name
    self%elements = elements
    self%atomtypes = atomtypes
    self%born = born
    self%epsilon = epsilon
    self%basis = basis
    self%polar = polar
    self%read_epsiloninf = read_epsiloninf
    self%epsilon0 = epsilon0
    self%lattvecs = lattvecs
    self%T = T
    self%autoisotopes = autoisotopes
    self%masses = masses
    self%gfactors = 0.0_r64
    self%twod = twod
    self%subs_masses = subs_masses
    self%subs_conc = subs_conc
    self%bound_length = bound_length
    self%thinfilm_height = thinfilm_height
    self%thinfilm_normal = thinfilm_normal
    self%defect_hosts = defect_hosts
    self%numdopants_types = numdopants_types
    
    if(product(numdopants_types) <= 0) then
       call exit_with_message('Number of dopant types must be a non-zero integer. Exiting.')   
    end if
    
    if(self%twod) then
       if(lattvecs(1,3) /= 0 .or. lattvecs(2,3) /= 0 .or. lattvecs(3,3) == 0) then
          call exit_with_message('For 2d systems, cross plane lattice vector must be &
               &of the for (0 0 h).')
       end if
       self%thickness = lattvecs(3,3)
       self%dim = 2.0_r64
    else
       self%dim = 3.0_r64
    end if

    !Dopant masses and concentrations
    allocate(self%dopant_masses(maxval(numdopants_types), self%numelements), &
         self%dopant_conc(maxval(numdopants_types), self%numelements))
    self%dopant_masses = 0.0_r64
    self%dopant_conc = 0.0_r64
    self%dopant_masses = dopant_masses(1:size(self%dopant_masses, 1), :)
    self%dopant_conc = dopant_conc(1:size(self%dopant_conc, 1), :)
    
    !Set high-frequency dielectric constant
    if(self%read_epsiloninf) then
       self%epsiloninf = epsiloninf
    else
       self%epsiloninf = trace(self%epsilon)/3.0_r64
    end if
    
    !If required, calculate isotopic average masses and g-factors
    if(autoisotopes) then
       call calculate_mavg_and_g(self%elements, self%masses, self%gfactors)
    end if

    !Calculate atomic basis in Cartesian coordinates
    self%basis_cart(:,:) = matmul(self%lattvecs,self%basis)
    
    !Calculate reciprocal lattice vectors and real and reciprocal cell volumes
    do i = 1, 3
       j = mod(i, 3) + 1
       k = mod(j, 3) + 1
       self%reclattvecs(:,i) = &
            cross_product(self%lattvecs(:, j), self%lattvecs(:, k))
    end do
    self%volume = abs(dot_product(self%lattvecs(:, 1),self%reclattvecs(:, 1)))
    self%volume_bz = twopi/self%volume
    self%reclattvecs(:,:) = self%volume_bz*self%reclattvecs(:,:)

    !Calculate the number of atoms of each type
    num_atomtypes(:) = 0_i64
    do i = 1, self%numelements
       do j = 1, self%numatoms
          if(self%atomtypes(j) == i) num_atomtypes(i) = num_atomtypes(i) + 1 
       end do
    end do
    
    !Convert number concentration of substitutions to percentage
    !of replaced host atoms.
    if(twod) then
       subs_perc = self%subs_conc*(1.0e-14_r64*self%volume/self%thickness)/num_atomtypes*100.0_r64
    else
       subs_perc = self%subs_conc*(1.0e-21_r64*self%volume)/num_atomtypes*100.0_r64
    end if
        
    !Calculate the mass variance parameters for the substitutions
    do i = 1, self%numelements
       !Impurity and host mixed mass
       subs_mavg = (subs_perc(i)*self%subs_masses(i) + &
            (100.0_r64 - subs_perc(i))*self%masses(i))/100.0_r64

       !g-factor
       self%subs_gfactors(i) = subs_perc(i)*(1.0_r64 - self%subs_masses(i)/subs_mavg)**2 + &
            (100.0_r64 - subs_perc(i))*(1.0_r64 - self%masses(i)/subs_mavg)**2
    end do
    self%subs_gfactors = self%subs_gfactors/100.0_r64
    
    !Print out crystal and reciprocal lattice information.
    if(this_image() == 1) then
       write(*, "(A, A)") 'Material: ', self%name

       if(any(self%defect_hosts > 0)) then
          write(*, "(A)") "Basis atom sites ready to host a substitution:"
          do i = 1, self%numelements
             write(*, '(A, A, A, I5)') " ", self%elements(i), " at site ", self%defect_hosts(i)
          end do
       end if
       
       if(self%autoisotopes) write(*,"(A)") 'Isotopic average of masses will be used.'
       do i = 1, self%numelements
          write(*,"(A, A, 1E16.8, A)") trim(self%elements(i)), " mass = ", self%masses(i), " u"
       end do

       if(any(self%subs_conc /= 0.0_r64)) then
          do i = 1, self%numelements
             write(*,"(A, A, 1E16.8, A)") &
                  trim(self%elements(i)), " substitution mass = ", self%subs_masses(i), " u"
          end do
          do i = 1, self%numelements
             write(*,"(A, A, 1E16.8, A)") &
                  trim(self%elements(i)), " substitution amount = ", subs_perc(i), " %"
          end do
       end if
       
       write(*,"(A)") 'Lattice vectors [nm]:'
       write(*,"(3(1E16.8,x))") self%lattvecs(:,1)
       write(*,"(3(1E16.8,x))") self%lattvecs(:,2)
       write(*,"(3(1E16.8,x))") self%lattvecs(:,3)
       write(*,"(A,(1E16.8,x),A)") 'Primitive cell volume =', self%volume, 'nm^3'

       write(*,"(A)") 'Reciprocal lattice vectors [1/nm]:'
       write(*,"(3(1E16.8,x))") self%reclattvecs(:,1)
       write(*,"(3(1E16.8,x))") self%reclattvecs(:,2)
       write(*,"(3(1E16.8,x))") self%reclattvecs(:,3)
       write(*,"(A,(1E16.8,x),A)") 'Brillouin zone volume =', self%volume_bz, '1/nm^3'
       if(self%twod) write(*,"(A)") 'System is 2d.'
       
       if(self%polar) then
          write(*,"(A)") 'System is polar.'
          write(*,"(A)") 'Dielectric tensor:'
          do i = 1, 3
             write(*,"(3(1E16.8,x))") self%epsilon(:,i)
          end do
          write(*,"(A)") 'Born effective charges:'
          do i = 1, self%numatoms
             write(*,"(A)") trim(self%elements(self%atomtypes(i)))
             do j = 1, 3
                write(*,"(3(1E16.8,x))") self%born(:,j,i)
             end do
          end do
          write(*,"(A,1E16.8)") 'Static dielectric (used for screening e-ch. imp. interactions) = ', self%epsilon0
          write(*,"(A,1E16.8)") 'High-frequency dielectric = ', self%epsiloninf
       end if
       write(*,"(A, F7.2, A)") 'Crystal temperature = ', self%T, ' K'
    end if
  end subroutine read_input_and_setup_crystal

  subroutine calculate_wavevectors_full(mesh, wavevecs, blocks, indexlist)
    !! Calculate wave vectors (crystal coords.) of the full Brillouin zone (FBZ)
    !!
    !! mesh is the array of number of points along the reciprocal lattice vectors
    !! wavevecs is the list of all the wave vectors

    integer(i64), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(i64), optional, intent(in) :: indexlist(:)
    real(r64), allocatable, intent(out) :: wavevecs(:,:)
    integer(i64) :: nwavevecs, ijk(3), i, imux

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
       call demux_vector(imux, ijk, mesh, 0_i64) !get 0-based (i,j,k) indices
       wavevecs(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine calculate_wavevectors_full

  subroutine calculate_mavg_and_g(elements, m, g)
    !! Compute the average mass of each element and its g-factor (Pearson
    !! deviation coefficient of the masses).
    !!
    !! This subroutine is adapted from ShengBTE.
    
    character(len=3), intent(in) :: elements(:)
    real(r64), intent(out) :: m(:), g(:)

    !Local variables
    integer(i64) :: i, niso, nelems, e
    character(len = 3) :: isotope_element(287)
    real(r64) :: isotope_mass(287)
    real(r64) :: isotope_abundance(287)

    nelems = size(elements)
    niso = 287
    
    ! Fill in isotope data.
    isotope_element(1) = "Ag"
    isotope_mass(1) = 106.905095_r64
    isotope_abundance(1) = 51.84_r64
    isotope_element(2) = "Ag"
    isotope_mass(2) = 108.904754_r64
    isotope_abundance(2) = 48.16_r64
    isotope_element(3) = "Al"
    isotope_mass(3) = 26.981541_r64
    isotope_abundance(3) = 100.0_r64
    isotope_element(4) = "Ar"
    isotope_mass(4) = 35.967546_r64
    isotope_abundance(4) = 0.34_r64
    isotope_element(5) = "Ar"
    isotope_mass(5) = 37.962732_r64
    isotope_abundance(5) = 0.063_r64
    isotope_element(6) = "Ar"
    isotope_mass(6) = 39.962383_r64
    isotope_abundance(6) = 99.6_r64
    isotope_element(7) = "As"
    isotope_mass(7) = 74.921596_r64
    isotope_abundance(7) = 100.0_r64
    isotope_element(8) = "Au"
    isotope_mass(8) = 196.96656_r64
    isotope_abundance(8) = 100.0_r64
    isotope_element(9) = "B"
    isotope_mass(9) = 10.012938_r64
    isotope_abundance(9) = 19.8_r64
    isotope_element(10) = "B"
    isotope_mass(10) = 11.009305_r64
    isotope_abundance(10) = 80.2_r64
    isotope_element(11) = "Ba"
    isotope_mass(11) = 129.906277_r64
    isotope_abundance(11) = 0.11_r64
    isotope_element(12) = "Ba"
    isotope_mass(12) = 131.905042_r64
    isotope_abundance(12) = 0.1_r64
    isotope_element(13) = "Ba"
    isotope_mass(13) = 133.90449_r64
    isotope_abundance(13) = 2.42_r64
    isotope_element(14) = "Ba"
    isotope_mass(14) = 134.905668_r64
    isotope_abundance(14) = 6.59_r64
    isotope_element(15) = "Ba"
    isotope_mass(15) = 135.904556_r64
    isotope_abundance(15) = 7.85_r64
    isotope_element(16) = "Ba"
    isotope_mass(16) = 136.905816_r64
    isotope_abundance(16) = 11.23_r64
    isotope_element(17) = "Ba"
    isotope_mass(17) = 137.905236_r64
    isotope_abundance(17) = 71.7_r64
    isotope_element(18) = "Be"
    isotope_mass(18) = 9.012183_r64
    isotope_abundance(18) = 100.0_r64
    isotope_element(19) = "Bi"
    isotope_mass(19) = 208.980388_r64
    isotope_abundance(19) = 100.0_r64
    isotope_element(20) = "Br"
    isotope_mass(20) = 78.918336_r64
    isotope_abundance(20) = 50.69_r64
    isotope_element(21) = "Br"
    isotope_mass(21) = 80.91629_r64
    isotope_abundance(21) = 49.31_r64
    isotope_element(22) = "C"
    isotope_mass(22) = 12.0_r64
    isotope_abundance(22) = 98.9_r64
    isotope_element(23) = "C"
    isotope_mass(23) = 13.003355_r64
    isotope_abundance(23) = 1.1_r64
    isotope_element(24) = "Ca"
    isotope_mass(24) = 39.962591_r64
    isotope_abundance(24) = 96.95_r64
    isotope_element(25) = "Ca"
    isotope_mass(25) = 41.958622_r64
    isotope_abundance(25) = 0.65_r64
    isotope_element(26) = "Ca"
    isotope_mass(26) = 42.95877_r64
    isotope_abundance(26) = 0.14_r64
    isotope_element(27) = "Ca"
    isotope_mass(27) = 43.955485_r64
    isotope_abundance(27) = 2.086_r64
    isotope_element(28) = "Ca"
    isotope_mass(28) = 45.953689_r64
    isotope_abundance(28) = 0.004_r64
    isotope_element(29) = "Ca"
    isotope_mass(29) = 47.952532_r64
    isotope_abundance(29) = 0.19_r64
    isotope_element(30) = "Cd"
    isotope_mass(30) = 105.906461_r64
    isotope_abundance(30) = 1.25_r64
    isotope_element(31) = "Cd"
    isotope_mass(31) = 107.904186_r64
    isotope_abundance(31) = 0.89_r64
    isotope_element(32) = "Cd"
    isotope_mass(32) = 109.903007_r64
    isotope_abundance(32) = 12.49_r64
    isotope_element(33) = "Cd"
    isotope_mass(33) = 110.904182_r64
    isotope_abundance(33) = 12.8_r64
    isotope_element(34) = "Cd"
    isotope_mass(34) = 111.902761_r64
    isotope_abundance(34) = 24.13_r64
    isotope_element(35) = "Cd"
    isotope_mass(35) = 112.904401_r64
    isotope_abundance(35) = 12.22_r64
    isotope_element(36) = "Cd"
    isotope_mass(36) = 113.903361_r64
    isotope_abundance(36) = 28.73_r64
    isotope_element(37) = "Cd"
    isotope_mass(37) = 115.904758_r64
    isotope_abundance(37) = 7.49_r64
    isotope_element(38) = "Ce"
    isotope_mass(38) = 135.90714_r64
    isotope_abundance(38) = 0.19_r64
    isotope_element(39) = "Ce"
    isotope_mass(39) = 137.905996_r64
    isotope_abundance(39) = 0.25_r64
    isotope_element(40) = "Ce"
    isotope_mass(40) = 139.905442_r64
    isotope_abundance(40) = 88.48_r64
    isotope_element(41) = "Ce"
    isotope_mass(41) = 141.909249_r64
    isotope_abundance(41) = 11.08_r64
    isotope_element(42) = "Cl"
    isotope_mass(42) = 34.968853_r64
    isotope_abundance(42) = 75.77_r64
    isotope_element(43) = "Cl"
    isotope_mass(43) = 36.965903_r64
    isotope_abundance(43) = 24.23_r64
    isotope_element(44) = "Co"
    isotope_mass(44) = 58.933198_r64
    isotope_abundance(44) = 100.0_r64
    isotope_element(45) = "Cr"
    isotope_mass(45) = 49.946046_r64
    isotope_abundance(45) = 4.35_r64
    isotope_element(46) = "Cr"
    isotope_mass(46) = 51.94051_r64
    isotope_abundance(46) = 83.79_r64
    isotope_element(47) = "Cr"
    isotope_mass(47) = 52.940651_r64
    isotope_abundance(47) = 9.5_r64
    isotope_element(48) = "Cr"
    isotope_mass(48) = 53.938882_r64
    isotope_abundance(48) = 2.36_r64
    isotope_element(49) = "Cs"
    isotope_mass(49) = 132.905433_r64
    isotope_abundance(49) = 100.0_r64
    isotope_element(50) = "Cu"
    isotope_mass(50) = 62.929599_r64
    isotope_abundance(50) = 69.17_r64
    isotope_element(51) = "Cu"
    isotope_mass(51) = 64.927792_r64
    isotope_abundance(51) = 30.83_r64
    isotope_element(52) = "Dy"
    isotope_mass(52) = 155.924287_r64
    isotope_abundance(52) = 0.06_r64
    isotope_element(53) = "Dy"
    isotope_mass(53) = 157.924412_r64
    isotope_abundance(53) = 0.1_r64
    isotope_element(54) = "Dy"
    isotope_mass(54) = 159.925203_r64
    isotope_abundance(54) = 2.34_r64
    isotope_element(55) = "Dy"
    isotope_mass(55) = 160.926939_r64
    isotope_abundance(55) = 18.9_r64
    isotope_element(56) = "Dy"
    isotope_mass(56) = 161.926805_r64
    isotope_abundance(56) = 25.5_r64
    isotope_element(57) = "Dy"
    isotope_mass(57) = 162.928737_r64
    isotope_abundance(57) = 24.9_r64
    isotope_element(58) = "Dy"
    isotope_mass(58) = 163.929183_r64
    isotope_abundance(58) = 28.2_r64
    isotope_element(59) = "Er"
    isotope_mass(59) = 161.928787_r64
    isotope_abundance(59) = 0.14_r64
    isotope_element(60) = "Er"
    isotope_mass(60) = 163.929211_r64
    isotope_abundance(60) = 1.61_r64
    isotope_element(61) = "Er"
    isotope_mass(61) = 165.930305_r64
    isotope_abundance(61) = 33.6_r64
    isotope_element(62) = "Er"
    isotope_mass(62) = 166.932061_r64
    isotope_abundance(62) = 22.95_r64
    isotope_element(63) = "Er"
    isotope_mass(63) = 167.932383_r64
    isotope_abundance(63) = 26.8_r64
    isotope_element(64) = "Er"
    isotope_mass(64) = 169.935476_r64
    isotope_abundance(64) = 14.9_r64
    isotope_element(65) = "Eu"
    isotope_mass(65) = 150.91986_r64
    isotope_abundance(65) = 47.8_r64
    isotope_element(66) = "Eu"
    isotope_mass(66) = 152.921243_r64
    isotope_abundance(66) = 52.2_r64
    isotope_element(67) = "F"
    isotope_mass(67) = 18.998403_r64
    isotope_abundance(67) = 100.0_r64
    isotope_element(68) = "Fe"
    isotope_mass(68) = 53.939612_r64
    isotope_abundance(68) = 5.8_r64
    isotope_element(69) = "Fe"
    isotope_mass(69) = 55.934939_r64
    isotope_abundance(69) = 91.72_r64
    isotope_element(70) = "Fe"
    isotope_mass(70) = 56.935396_r64
    isotope_abundance(70) = 2.2_r64
    isotope_element(71) = "Fe"
    isotope_mass(71) = 57.933278_r64
    isotope_abundance(71) = 0.28_r64
    isotope_element(72) = "Ga"
    isotope_mass(72) = 68.925581_r64
    isotope_abundance(72) = 60.1_r64
    isotope_element(73) = "Ga"
    isotope_mass(73) = 70.924701_r64
    isotope_abundance(73) = 39.9_r64
    isotope_element(74) = "Gd"
    isotope_mass(74) = 151.919803_r64
    isotope_abundance(74) = 0.2_r64
    isotope_element(75) = "Gd"
    isotope_mass(75) = 153.920876_r64
    isotope_abundance(75) = 2.18_r64
    isotope_element(76) = "Gd"
    isotope_mass(76) = 154.822629_r64
    isotope_abundance(76) = 14.8_r64
    isotope_element(77) = "Gd"
    isotope_mass(77) = 155.92213_r64
    isotope_abundance(77) = 20.47_r64
    isotope_element(78) = "Gd"
    isotope_mass(78) = 156.923967_r64
    isotope_abundance(78) = 15.65_r64
    isotope_element(79) = "Gd"
    isotope_mass(79) = 157.924111_r64
    isotope_abundance(79) = 24.84_r64
    isotope_element(80) = "Gd"
    isotope_mass(80) = 159.927061_r64
    isotope_abundance(80) = 21.86_r64
    isotope_element(81) = "Ge"
    isotope_mass(81) = 69.92425_r64
    isotope_abundance(81) = 20.5_r64
    isotope_element(82) = "Ge"
    isotope_mass(82) = 71.92208_r64
    isotope_abundance(82) = 27.4_r64
    isotope_element(83) = "Ge"
    isotope_mass(83) = 72.923464_r64
    isotope_abundance(83) = 7.8_r64
    isotope_element(84) = "Ge"
    isotope_mass(84) = 73.921179_r64
    isotope_abundance(84) = 36.5_r64
    isotope_element(85) = "Ge"
    isotope_mass(85) = 75.921403_r64
    isotope_abundance(85) = 7.8_r64
    isotope_element(86) = "H"
    isotope_mass(86) = 1.007825_r64
    isotope_abundance(86) = 99.99_r64
    isotope_element(87) = "H"
    isotope_mass(87) = 2.014102_r64
    isotope_abundance(87) = 0.015_r64
    isotope_element(88) = "He"
    isotope_mass(88) = 3.016029_r64
    isotope_abundance(88) = 0.0001_r64
    isotope_element(89) = "He"
    isotope_mass(89) = 4.002603_r64
    isotope_abundance(89) = 100.0_r64
    isotope_element(90) = "Hf"
    isotope_mass(90) = 173.940065_r64
    isotope_abundance(90) = 0.16_r64
    isotope_element(91) = "Hf"
    isotope_mass(91) = 175.94142_r64
    isotope_abundance(91) = 5.2_r64
    isotope_element(92) = "Hf"
    isotope_mass(92) = 176.943233_r64
    isotope_abundance(92) = 18.6_r64
    isotope_element(93) = "Hf"
    isotope_mass(93) = 177.94371_r64
    isotope_abundance(93) = 27.1_r64
    isotope_element(94) = "Hf"
    isotope_mass(94) = 178.945827_r64
    isotope_abundance(94) = 13.74_r64
    isotope_element(95) = "Hf"
    isotope_mass(95) = 179.946561_r64
    isotope_abundance(95) = 35.2_r64
    isotope_element(96) = "Hg"
    isotope_mass(96) = 195.965812_r64
    isotope_abundance(96) = 0.15_r64
    isotope_element(97) = "Hg"
    isotope_mass(97) = 197.96676_r64
    isotope_abundance(97) = 10.1_r64
    isotope_element(98) = "Hg"
    isotope_mass(98) = 198.968269_r64
    isotope_abundance(98) = 17.0_r64
    isotope_element(99) = "Hg"
    isotope_mass(99) = 199.968316_r64
    isotope_abundance(99) = 23.1_r64
    isotope_element(100) = "Hg"
    isotope_mass(100) = 200.970293_r64
    isotope_abundance(100) = 13.2_r64
    isotope_element(101) = "Hg"
    isotope_mass(101) = 201.970632_r64
    isotope_abundance(101) = 29.65_r64
    isotope_element(102) = "Hg"
    isotope_mass(102) = 203.973481_r64
    isotope_abundance(102) = 6.8_r64
    isotope_element(103) = "Ho"
    isotope_mass(103) = 164.930332_r64
    isotope_abundance(103) = 100.0_r64
    isotope_element(104) = "I"
    isotope_mass(104) = 126.904477_r64
    isotope_abundance(104) = 100.0_r64
    isotope_element(105) = "In"
    isotope_mass(105) = 112.904056_r64
    isotope_abundance(105) = 4.3_r64
    isotope_element(106) = "In"
    isotope_mass(106) = 114.903875_r64
    isotope_abundance(106) = 95.7_r64
    isotope_element(107) = "Ir"
    isotope_mass(107) = 190.960603_r64
    isotope_abundance(107) = 37.3_r64
    isotope_element(108) = "Ir"
    isotope_mass(108) = 192.962942_r64
    isotope_abundance(108) = 62.7_r64
    isotope_element(109) = "K"
    isotope_mass(109) = 38.963708_r64
    isotope_abundance(109) = 93.2_r64
    isotope_element(110) = "K"
    isotope_mass(110) = 39.963999_r64
    isotope_abundance(110) = 0.012_r64
    isotope_element(111) = "K"
    isotope_mass(111) = 40.961825_r64
    isotope_abundance(111) = 6.73_r64
    isotope_element(112) = "Kr"
    isotope_mass(112) = 77.920397_r64
    isotope_abundance(112) = 0.35_r64
    isotope_element(113) = "Kr"
    isotope_mass(113) = 79.916375_r64
    isotope_abundance(113) = 2.25_r64
    isotope_element(114) = "Kr"
    isotope_mass(114) = 81.913483_r64
    isotope_abundance(114) = 11.6_r64
    isotope_element(115) = "Kr"
    isotope_mass(115) = 82.914134_r64
    isotope_abundance(115) = 11.5_r64
    isotope_element(116) = "Kr"
    isotope_mass(116) = 83.911506_r64
    isotope_abundance(116) = 57.0_r64
    isotope_element(117) = "Kr"
    isotope_mass(117) = 85.910614_r64
    isotope_abundance(117) = 17.3_r64
    isotope_element(118) = "La"
    isotope_mass(118) = 137.907114_r64
    isotope_abundance(118) = 0.09_r64
    isotope_element(119) = "La"
    isotope_mass(119) = 138.906355_r64
    isotope_abundance(119) = 99.91_r64
    isotope_element(120) = "Li"
    isotope_mass(120) = 6.015123_r64
    isotope_abundance(120) = 7.42_r64
    isotope_element(121) = "Li"
    isotope_mass(121) = 7.016005_r64
    isotope_abundance(121) = 92.58_r64
    isotope_element(122) = "Lu"
    isotope_mass(122) = 174.940785_r64
    isotope_abundance(122) = 97.4_r64
    isotope_element(123) = "Lu"
    isotope_mass(123) = 175.942694_r64
    isotope_abundance(123) = 2.6_r64
    isotope_element(124) = "Mg"
    isotope_mass(124) = 23.985045_r64
    isotope_abundance(124) = 78.9_r64
    isotope_element(125) = "Mg"
    isotope_mass(125) = 24.985839_r64
    isotope_abundance(125) = 10.0_r64
    isotope_element(126) = "Mg"
    isotope_mass(126) = 25.982595_r64
    isotope_abundance(126) = 11.1_r64
    isotope_element(127) = "Mn"
    isotope_mass(127) = 54.938046_r64
    isotope_abundance(127) = 100.0_r64
    isotope_element(128) = "Mo"
    isotope_mass(128) = 91.906809_r64
    isotope_abundance(128) = 14.84_r64
    isotope_element(129) = "Mo"
    isotope_mass(129) = 93.905086_r64
    isotope_abundance(129) = 9.25_r64
    isotope_element(130) = "Mo"
    isotope_mass(130) = 94.905838_r64
    isotope_abundance(130) = 15.92_r64
    isotope_element(131) = "Mo"
    isotope_mass(131) = 95.904676_r64
    isotope_abundance(131) = 16.68_r64
    isotope_element(132) = "Mo"
    isotope_mass(132) = 96.906018_r64
    isotope_abundance(132) = 9.55_r64
    isotope_element(133) = "Mo"
    isotope_mass(133) = 97.905405_r64
    isotope_abundance(133) = 24.13_r64
    isotope_element(134) = "Mo"
    isotope_mass(134) = 99.907473_r64
    isotope_abundance(134) = 9.63_r64
    isotope_element(135) = "N"
    isotope_mass(135) = 14.003074_r64
    isotope_abundance(135) = 99.63_r64
    isotope_element(136) = "N"
    isotope_mass(136) = 15.000109_r64
    isotope_abundance(136) = 0.37_r64
    isotope_element(137) = "Na"
    isotope_mass(137) = 22.98977_r64
    isotope_abundance(137) = 100.0_r64
    isotope_element(138) = "Nb"
    isotope_mass(138) = 92.906378_r64
    isotope_abundance(138) = 100.0_r64
    isotope_element(139) = "Nd"
    isotope_mass(139) = 141.907731_r64
    isotope_abundance(139) = 27.13_r64
    isotope_element(140) = "Nd"
    isotope_mass(140) = 142.909823_r64
    isotope_abundance(140) = 12.18_r64
    isotope_element(141) = "Nd"
    isotope_mass(141) = 143.910096_r64
    isotope_abundance(141) = 23.8_r64
    isotope_element(142) = "Nd"
    isotope_mass(142) = 144.912582_r64
    isotope_abundance(142) = 8.3_r64
    isotope_element(143) = "Nd"
    isotope_mass(143) = 145.913126_r64
    isotope_abundance(143) = 17.19_r64
    isotope_element(144) = "Nd"
    isotope_mass(144) = 147.916901_r64
    isotope_abundance(144) = 5.76_r64
    isotope_element(145) = "Nd"
    isotope_mass(145) = 149.9209_r64
    isotope_abundance(145) = 5.64_r64
    isotope_element(146) = "Ne"
    isotope_mass(146) = 19.992439_r64
    isotope_abundance(146) = 90.6_r64
    isotope_element(147) = "Ne"
    isotope_mass(147) = 20.993845_r64
    isotope_abundance(147) = 0.26_r64
    isotope_element(148) = "Ne"
    isotope_mass(148) = 21.991384_r64
    isotope_abundance(148) = 9.2_r64
    isotope_element(149) = "Ni"
    isotope_mass(149) = 57.935347_r64
    isotope_abundance(149) = 68.27_r64
    isotope_element(150) = "Ni"
    isotope_mass(150) = 59.930789_r64
    isotope_abundance(150) = 26.1_r64
    isotope_element(151) = "Ni"
    isotope_mass(151) = 60.931059_r64
    isotope_abundance(151) = 1.13_r64
    isotope_element(152) = "Ni"
    isotope_mass(152) = 61.928346_r64
    isotope_abundance(152) = 3.59_r64
    isotope_element(153) = "Ni"
    isotope_mass(153) = 63.927968_r64
    isotope_abundance(153) = 0.91_r64
    isotope_element(154) = "O"
    isotope_mass(154) = 15.994915_r64
    isotope_abundance(154) = 99.76_r64
    isotope_element(155) = "O"
    isotope_mass(155) = 16.999131_r64
    isotope_abundance(155) = 0.038_r64
    isotope_element(156) = "O"
    isotope_mass(156) = 17.999159_r64
    isotope_abundance(156) = 0.2_r64
    isotope_element(157) = "Os"
    isotope_mass(157) = 183.952514_r64
    isotope_abundance(157) = 0.02_r64
    isotope_element(158) = "Os"
    isotope_mass(158) = 185.953852_r64
    isotope_abundance(158) = 1.58_r64
    isotope_element(159) = "Os"
    isotope_mass(159) = 186.955762_r64
    isotope_abundance(159) = 1.6_r64
    isotope_element(160) = "Os"
    isotope_mass(160) = 187.95585_r64
    isotope_abundance(160) = 13.3_r64
    isotope_element(161) = "Os"
    isotope_mass(161) = 188.958156_r64
    isotope_abundance(161) = 16.1_r64
    isotope_element(162) = "Os"
    isotope_mass(162) = 189.958455_r64
    isotope_abundance(162) = 26.4_r64
    isotope_element(163) = "Os"
    isotope_mass(163) = 191.961487_r64
    isotope_abundance(163) = 41.0_r64
    isotope_element(164) = "P"
    isotope_mass(164) = 30.973763_r64
    isotope_abundance(164) = 100.0_r64
    isotope_element(165) = "Pb"
    isotope_mass(165) = 203.973037_r64
    isotope_abundance(165) = 1.4_r64
    isotope_element(166) = "Pb"
    isotope_mass(166) = 205.974455_r64
    isotope_abundance(166) = 24.1_r64
    isotope_element(167) = "Pb"
    isotope_mass(167) = 206.975885_r64
    isotope_abundance(167) = 22.1_r64
    isotope_element(168) = "Pb"
    isotope_mass(168) = 207.976641_r64
    isotope_abundance(168) = 52.4_r64
    isotope_element(169) = "Pd"
    isotope_mass(169) = 101.905609_r64
    isotope_abundance(169) = 1.02_r64
    isotope_element(170) = "Pd"
    isotope_mass(170) = 103.904026_r64
    isotope_abundance(170) = 11.14_r64
    isotope_element(171) = "Pd"
    isotope_mass(171) = 104.905075_r64
    isotope_abundance(171) = 22.33_r64
    isotope_element(172) = "Pd"
    isotope_mass(172) = 105.903475_r64
    isotope_abundance(172) = 27.33_r64
    isotope_element(173) = "Pd"
    isotope_mass(173) = 107.903894_r64
    isotope_abundance(173) = 26.46_r64
    isotope_element(174) = "Pd"
    isotope_mass(174) = 109.905169_r64
    isotope_abundance(174) = 11.72_r64
    isotope_element(175) = "Pr"
    isotope_mass(175) = 140.907657_r64
    isotope_abundance(175) = 100.0_r64
    isotope_element(176) = "Pt"
    isotope_mass(176) = 189.959937_r64
    isotope_abundance(176) = 0.01_r64
    isotope_element(177) = "Pt"
    isotope_mass(177) = 191.961049_r64
    isotope_abundance(177) = 0.79_r64
    isotope_element(178) = "Pt"
    isotope_mass(178) = 193.962679_r64
    isotope_abundance(178) = 32.9_r64
    isotope_element(179) = "Pt"
    isotope_mass(179) = 194.964785_r64
    isotope_abundance(179) = 33.8_r64
    isotope_element(180) = "Pt"
    isotope_mass(180) = 195.964947_r64
    isotope_abundance(180) = 25.3_r64
    isotope_element(181) = "Pt"
    isotope_mass(181) = 197.967879_r64
    isotope_abundance(181) = 7.2_r64
    isotope_element(182) = "Rb"
    isotope_mass(182) = 84.9118_r64
    isotope_abundance(182) = 72.17_r64
    isotope_element(183) = "Rb"
    isotope_mass(183) = 86.909184_r64
    isotope_abundance(183) = 27.84_r64
    isotope_element(184) = "Re"
    isotope_mass(184) = 184.952977_r64
    isotope_abundance(184) = 37.4_r64
    isotope_element(185) = "Re"
    isotope_mass(185) = 186.955765_r64
    isotope_abundance(185) = 62.6_r64
    isotope_element(186) = "Rh"
    isotope_mass(186) = 102.905503_r64
    isotope_abundance(186) = 100.0_r64
    isotope_element(187) = "Ru"
    isotope_mass(187) = 95.907596_r64
    isotope_abundance(187) = 5.52_r64
    isotope_element(188) = "Ru"
    isotope_mass(188) = 97.905287_r64
    isotope_abundance(188) = 1.88_r64
    isotope_element(189) = "Ru"
    isotope_mass(189) = 98.905937_r64
    isotope_abundance(189) = 12.7_r64
    isotope_element(190) = "Ru"
    isotope_mass(190) = 99.904218_r64
    isotope_abundance(190) = 12.6_r64
    isotope_element(191) = "Ru"
    isotope_mass(191) = 100.905581_r64
    isotope_abundance(191) = 17.0_r64
    isotope_element(192) = "Ru"
    isotope_mass(192) = 101.904348_r64
    isotope_abundance(192) = 31.6_r64
    isotope_element(193) = "Ru"
    isotope_mass(193) = 103.905422_r64
    isotope_abundance(193) = 18.7_r64
    isotope_element(194) = "S"
    isotope_mass(194) = 31.972072_r64
    isotope_abundance(194) = 95.02_r64
    isotope_element(195) = "S"
    isotope_mass(195) = 32.971459_r64
    isotope_abundance(195) = 0.75_r64
    isotope_element(196) = "S"
    isotope_mass(196) = 33.967868_r64
    isotope_abundance(196) = 4.21_r64
    isotope_element(197) = "S"
    isotope_mass(197) = 35.967079_r64
    isotope_abundance(197) = 0.02_r64
    isotope_element(198) = "Sb"
    isotope_mass(198) = 120.903824_r64
    isotope_abundance(198) = 57.3_r64
    isotope_element(199) = "Sb"
    isotope_mass(199) = 122.904222_r64
    isotope_abundance(199) = 42.7_r64
    isotope_element(200) = "Sc"
    isotope_mass(200) = 44.955914_r64
    isotope_abundance(200) = 100.0_r64
    isotope_element(201) = "Se"
    isotope_mass(201) = 73.922477_r64
    isotope_abundance(201) = 0.9_r64
    isotope_element(202) = "Se"
    isotope_mass(202) = 75.919207_r64
    isotope_abundance(202) = 9.0_r64
    isotope_element(203) = "Se"
    isotope_mass(203) = 76.919908_r64
    isotope_abundance(203) = 7.6_r64
    isotope_element(204) = "Se"
    isotope_mass(204) = 77.917304_r64
    isotope_abundance(204) = 23.5_r64
    isotope_element(205) = "Se"
    isotope_mass(205) = 79.916521_r64
    isotope_abundance(205) = 49.6_r64
    isotope_element(206) = "Se"
    isotope_mass(206) = 81.916709_r64
    isotope_abundance(206) = 9.4_r64
    isotope_element(207) = "Si"
    isotope_mass(207) = 27.976928_r64
    isotope_abundance(207) = 92.23_r64
    isotope_element(208) = "Si"
    isotope_mass(208) = 28.976496_r64
    isotope_abundance(208) = 4.67_r64
    isotope_element(209) = "Si"
    isotope_mass(209) = 29.973772_r64
    isotope_abundance(209) = 3.1_r64
    isotope_element(210) = "Sm"
    isotope_mass(210) = 143.912009_r64
    isotope_abundance(210) = 3.1_r64
    isotope_element(211) = "Sm"
    isotope_mass(211) = 146.914907_r64
    isotope_abundance(211) = 15.0_r64
    isotope_element(212) = "Sm"
    isotope_mass(212) = 147.914832_r64
    isotope_abundance(212) = 11.3_r64
    isotope_element(213) = "Sm"
    isotope_mass(213) = 148.917193_r64
    isotope_abundance(213) = 13.8_r64
    isotope_element(214) = "Sm"
    isotope_mass(214) = 149.917285_r64
    isotope_abundance(214) = 7.4_r64
    isotope_element(215) = "Sm"
    isotope_mass(215) = 151.919741_r64
    isotope_abundance(215) = 26.7_r64
    isotope_element(216) = "Sm"
    isotope_mass(216) = 153.922218_r64
    isotope_abundance(216) = 22.7_r64
    isotope_element(217) = "Sn"
    isotope_mass(217) = 111.904826_r64
    isotope_abundance(217) = 0.97_r64
    isotope_element(218) = "Sn"
    isotope_mass(218) = 113.902784_r64
    isotope_abundance(218) = 0.65_r64
    isotope_element(219) = "Sn"
    isotope_mass(219) = 114.903348_r64
    isotope_abundance(219) = 0.36_r64
    isotope_element(220) = "Sn"
    isotope_mass(220) = 115.901744_r64
    isotope_abundance(220) = 14.7_r64
    isotope_element(221) = "Sn"
    isotope_mass(221) = 116.902954_r64
    isotope_abundance(221) = 7.7_r64
    isotope_element(222) = "Sn"
    isotope_mass(222) = 117.901607_r64
    isotope_abundance(222) = 24.3_r64
    isotope_element(223) = "Sn"
    isotope_mass(223) = 118.90331_r64
    isotope_abundance(223) = 8.6_r64
    isotope_element(224) = "Sn"
    isotope_mass(224) = 119.902199_r64
    isotope_abundance(224) = 32.4_r64
    isotope_element(225) = "Sn"
    isotope_mass(225) = 121.90344_r64
    isotope_abundance(225) = 4.6_r64
    isotope_element(226) = "Sn"
    isotope_mass(226) = 123.905271_r64
    isotope_abundance(226) = 5.6_r64
    isotope_element(227) = "Sr"
    isotope_mass(227) = 83.913428_r64
    isotope_abundance(227) = 0.56_r64
    isotope_element(228) = "Sr"
    isotope_mass(228) = 85.909273_r64
    isotope_abundance(228) = 9.86_r64
    isotope_element(229) = "Sr"
    isotope_mass(229) = 86.908902_r64
    isotope_abundance(229) = 7.0_r64
    isotope_element(230) = "Sr"
    isotope_mass(230) = 87.905625_r64
    isotope_abundance(230) = 82.58_r64
    isotope_element(231) = "Ta"
    isotope_mass(231) = 179.947489_r64
    isotope_abundance(231) = 0.012_r64
    isotope_element(232) = "Ta"
    isotope_mass(232) = 180.948014_r64
    isotope_abundance(232) = 99.99_r64
    isotope_element(233) = "Tb"
    isotope_mass(233) = 158.92535_r64
    isotope_abundance(233) = 100.0_r64
    isotope_element(234) = "Te"
    isotope_mass(234) = 119.904021_r64
    isotope_abundance(234) = 0.096_r64
    isotope_element(235) = "Te"
    isotope_mass(235) = 121.903055_r64
    isotope_abundance(235) = 2.6_r64
    isotope_element(236) = "Te"
    isotope_mass(236) = 122.904278_r64
    isotope_abundance(236) = 0.91_r64
    isotope_element(237) = "Te"
    isotope_mass(237) = 123.902825_r64
    isotope_abundance(237) = 4.82_r64
    isotope_element(238) = "Te"
    isotope_mass(238) = 124.904435_r64
    isotope_abundance(238) = 7.14_r64
    isotope_element(239) = "Te"
    isotope_mass(239) = 125.90331_r64
    isotope_abundance(239) = 18.95_r64
    isotope_element(240) = "Te"
    isotope_mass(240) = 127.904464_r64
    isotope_abundance(240) = 31.69_r64
    isotope_element(241) = "Te"
    isotope_mass(241) = 129.906229_r64
    isotope_abundance(241) = 33.8_r64
    isotope_element(242) = "Th"
    isotope_mass(242) = 232.038054_r64
    isotope_abundance(242) = 100.0_r64
    isotope_element(243) = "Ti"
    isotope_mass(243) = 45.952633_r64
    isotope_abundance(243) = 8.0_r64
    isotope_element(244) = "Ti"
    isotope_mass(244) = 46.951765_r64
    isotope_abundance(244) = 7.3_r64
    isotope_element(245) = "Ti"
    isotope_mass(245) = 47.947947_r64
    isotope_abundance(245) = 73.8_r64
    isotope_element(246) = "Ti"
    isotope_mass(246) = 48.947871_r64
    isotope_abundance(246) = 5.5_r64
    isotope_element(247) = "Ti"
    isotope_mass(247) = 49.944786_r64
    isotope_abundance(247) = 5.4_r64
    isotope_element(248) = "Tl"
    isotope_mass(248) = 202.972336_r64
    isotope_abundance(248) = 29.52_r64
    isotope_element(249) = "Tl"
    isotope_mass(249) = 204.97441_r64
    isotope_abundance(249) = 70.48_r64
    isotope_element(250) = "Tm"
    isotope_mass(250) = 168.934225_r64
    isotope_abundance(250) = 100.0_r64
    isotope_element(251) = "U"
    isotope_mass(251) = 234.040947_r64
    isotope_abundance(251) = 0.006_r64
    isotope_element(252) = "U"
    isotope_mass(252) = 235.043925_r64
    isotope_abundance(252) = 0.72_r64
    isotope_element(253) = "U"
    isotope_mass(253) = 238.050786_r64
    isotope_abundance(253) = 99.27_r64
    isotope_element(254) = "V"
    isotope_mass(254) = 49.947161_r64
    isotope_abundance(254) = 0.25_r64
    isotope_element(255) = "V"
    isotope_mass(255) = 50.943963_r64
    isotope_abundance(255) = 99.75_r64
    isotope_element(256) = "W"
    isotope_mass(256) = 179.946727_r64
    isotope_abundance(256) = 0.13_r64
    isotope_element(257) = "W"
    isotope_mass(257) = 181.948225_r64
    isotope_abundance(257) = 26.3_r64
    isotope_element(258) = "W"
    isotope_mass(258) = 182.950245_r64
    isotope_abundance(258) = 14.3_r64
    isotope_element(259) = "W"
    isotope_mass(259) = 183.950953_r64
    isotope_abundance(259) = 30.67_r64
    isotope_element(260) = "W"
    isotope_mass(260) = 185.954377_r64
    isotope_abundance(260) = 28.6_r64
    isotope_element(261) = "Xe"
    isotope_mass(261) = 123.905894_r64
    isotope_abundance(261) = 0.1_r64
    isotope_element(262) = "Xe"
    isotope_mass(262) = 125.904281_r64
    isotope_abundance(262) = 0.09_r64
    isotope_element(263) = "Xe"
    isotope_mass(263) = 127.903531_r64
    isotope_abundance(263) = 1.91_r64
    isotope_element(264) = "Xe"
    isotope_mass(264) = 128.90478_r64
    isotope_abundance(264) = 26.4_r64
    isotope_element(265) = "Xe"
    isotope_mass(265) = 129.90351_r64
    isotope_abundance(265) = 4.1_r64
    isotope_element(266) = "Xe"
    isotope_mass(266) = 130.905076_r64
    isotope_abundance(266) = 21.2_r64
    isotope_element(267) = "Xe"
    isotope_mass(267) = 131.904148_r64
    isotope_abundance(267) = 26.9_r64
    isotope_element(268) = "Xe"
    isotope_mass(268) = 133.905395_r64
    isotope_abundance(268) = 10.4_r64
    isotope_element(269) = "Xe"
    isotope_mass(269) = 135.907219_r64
    isotope_abundance(269) = 8.9_r64
    isotope_element(270) = "Y"
    isotope_mass(270) = 88.905856_r64
    isotope_abundance(270) = 100.0_r64
    isotope_element(271) = "Yb"
    isotope_mass(271) = 167.933908_r64
    isotope_abundance(271) = 0.13_r64
    isotope_element(272) = "Yb"
    isotope_mass(272) = 169.934774_r64
    isotope_abundance(272) = 3.05_r64
    isotope_element(273) = "Yb"
    isotope_mass(273) = 170.936338_r64
    isotope_abundance(273) = 14.3_r64
    isotope_element(274) = "Yb"
    isotope_mass(274) = 171.936393_r64
    isotope_abundance(274) = 21.9_r64
    isotope_element(275) = "Yb"
    isotope_mass(275) = 172.938222_r64
    isotope_abundance(275) = 16.12_r64
    isotope_element(276) = "Yb"
    isotope_mass(276) = 173.938873_r64
    isotope_abundance(276) = 31.8_r64
    isotope_element(277) = "Yb"
    isotope_mass(277) = 175.942576_r64
    isotope_abundance(277) = 12.7_r64
    isotope_element(278) = "Zn"
    isotope_mass(278) = 63.929145_r64
    isotope_abundance(278) = 48.6_r64
    isotope_element(279) = "Zn"
    isotope_mass(279) = 65.926035_r64
    isotope_abundance(279) = 27.9_r64
    isotope_element(280) = "Zn"
    isotope_mass(280) = 66.927129_r64
    isotope_abundance(280) = 4.1_r64
    isotope_element(281) = "Zn"
    isotope_mass(281) = 67.924846_r64
    isotope_abundance(281) = 18.8_r64
    isotope_element(282) = "Zn"
    isotope_mass(282) = 69.925325_r64
    isotope_abundance(282) = 0.6_r64
    isotope_element(283) = "Zr"
    isotope_mass(283) = 89.904708_r64
    isotope_abundance(283) = 51.45_r64
    isotope_element(284) = "Zr"
    isotope_mass(284) = 90.905644_r64
    isotope_abundance(284) = 11.27_r64
    isotope_element(285) = "Zr"
    isotope_mass(285) = 91.905039_r64
    isotope_abundance(285) = 17.17_r64
    isotope_element(286) = "Zr"
    isotope_mass(286) = 93.906319_r64
    isotope_abundance(286) = 17.33_r64
    isotope_element(287) = "Zr"
    isotope_mass(287) = 95.908272_r64
    isotope_abundance(287) = 2.78_r64
    
    do e = 1, nelems
       m(e) = 0.0_r64
       do i = 1, niso
          if(isotope_element(i) .eq. elements(e)) then
             m(e) = m(e) + isotope_mass(i)*isotope_abundance(i)
          end if
       end do
       m(e) = m(e)/100.0_r64
       g(e) = 0.0_r64
       do i = 1, niso
          if(isotope_element(i) .eq. elements(e)) then
             g(e) = g(e) + isotope_abundance(i)*&
                  (1.0_r64 - isotope_mass(i)/m(e))**2
          end if
       end do
       g(e) = g(e)/100.0_r64
    end do
  end subroutine calculate_mavg_and_g
end module crystal_module

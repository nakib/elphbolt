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

  !Captain's log. August 11, 2023.
  !Ideally these should be immutables defined once and for all
  !in params.f90. Also, these should be stored in a dictionary
  !but I'm still shopping for a suitable one.
  character(len = 3), allocatable :: isotope_element(:)
  real(r64), allocatable :: isotope_mass(:)
  real(r64), allocatable :: isotope_abundance(:)
  
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
     logical :: VCA
     !! Use isotopic mix for masses?
     logical :: DIB
     !! Use dominant isotopes for masses?
     real(r64), allocatable :: gfactors_VCA(:)
     !! g-factors in the virtual crystal approximation (VCA).
     real(r64), allocatable :: gfactors_DIB(:)
     !! g-factors in the dominant isotope background approach (DIB).
     real(r64), allocatable :: subs_masses(:)
     !! Masses of the substitutional atoms [D]
     real(r64), allocatable :: subs_conc(:)
     !! Concentration of the substitutional atoms in cm^-3 [D]
     real(r64), allocatable :: subs_gfactors(:)
     !! g-factors for the substitutional defects. [D]
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
    integer(i64), allocatable :: atomtypes(:), num_atomtypes(:), numdopants_types(:)
    real(r64), allocatable :: masses(:), born(:,:,:), basis(:,:), &
         basis_cart(:,:), subs_perc(:), subs_masses(:), subs_conc(:), &
         dopant_masses(:, :), dopant_conc(:, :)
    real(r64) :: epsilon(3,3), lattvecs(3,3), T, &
         epsilon0, epsiloninf, subs_mavg, bound_length, thinfilm_height
    character(len=3), allocatable :: elements(:)
    character(len=100) :: name
    character(1) :: thinfilm_normal
    logical :: polar, VCA, DIB, read_epsiloninf, twod
    
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ name, elements, atomtypes, basis, lattvecs, &
         polar, born, epsilon, read_epsiloninf, epsilon0, epsiloninf, &
         masses, T, VCA, DIB, twod, subs_masses, subs_conc, bound_length, &
         numdopants_types, dopant_masses, dopant_conc, thinfilm_height, &
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
         num_atomtypes(numelements), numdopants_types(numelements), &
         dopant_masses(10, numelements), dopant_conc(10, numelements))
    allocate(self%elements(self%numelements), self%atomtypes(self%numatoms), &
         self%born(3,3,self%numatoms), self%masses(self%numelements), &
         self%gfactors_VCA(self%numelements), self%gfactors_DIB(self%numelements), &
         self%basis(3,self%numatoms), self%basis_cart(3,self%numatoms), &
         self%subs_masses(self%numelements), self%subs_conc(self%numelements), &
         self%subs_gfactors(self%numelements))
    
    !Read crystal_info
    name = trim(adjustl('Crystal'))
    elements = 'X'
    atomtypes = 0
    masses = -1.0_r64
    VCA = .false.
    DIB = .false.
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
    bound_length = 1.e12_r64 !mm, practically inifinity
    thinfilm_height = 1.e12_r64 !mm, practically inifinity
    thinfilm_normal = 'z'
    numdopants_types = [1, 1]
    dopant_masses = 0.0_r64
    dopant_conc = 0.0_r64
    read(1, nml = crystal_info)
    if(VCA .and. DIB) &
         call exit_with_message("'VCA' and 'DIB' can't both be true.")
    if(VCA) DIB = .false.
    if(any(atomtypes < 1)) &
         call exit_with_message("Atomtypes can't be negative.")
    if(T < 0.0_r64) &
         call exit_with_message("Negative temperatures not allowed here.")
    if(.not. DIB .and. .not. VCA .and. any(masses < 0)) then
       call exit_with_message("Provide valid basis atom masses for 'VCA' and 'DIB' false.")
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
    self%VCA = VCA
    self%DIB = DIB
    if(.not. VCA .and. .not. DIB) self%masses = masses
    self%gfactors_VCA = 0.0_r64
    self%gfactors_DIB = 0.0_r64
    self%twod = twod
    self%subs_masses = subs_masses
    self%subs_conc = subs_conc
    self%bound_length = bound_length
    self%thinfilm_height = thinfilm_height
    self%thinfilm_normal = thinfilm_normal
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

    !Generate periodic table
    call create_periodic_table(isotope_element, isotope_mass, isotope_abundance)

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
    if(self%VCA) &
         call calculate_mavg_and_g(self%elements, self%masses, self%gfactors_VCA)
    if(self%DIB) &
         call calculate_g_DIB(self%elements, self%masses, self%gfactors_DIB)

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
       
       if(self%VCA) write(*,"(A)") 'Isotopic average of masses (VCA) will be used.'
       if(self%DIB) write(*,"(A)") 'Dominant isotopic masses (DIB) will be used.'
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

    nelems = size(elements)
    niso = 287
        
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

  subroutine calculate_g_DIB(elements, m, g)
    !! Find the dominant isotopic masses from the periodic table
    !! and the average mass perturbation for use in the DIB-1st Born ph-iso scattering theory.
    
    character(len=3), intent(in) :: elements(:)
    real(r64), intent(out) :: m(:), g(:)

    !Local variables
    integer :: i, total_numiso, nelems, e, iso_count, numiso
    real(r64), allocatable :: these_iso_masses(:), these_iso_abundances(:)

    nelems = size(elements)
    total_numiso = 287

    do e = 1, nelems
       !1st pass to get the number of isotopes of this element
       numiso = 0
       do i = 1, total_numiso
          if(isotope_element(i) == elements(e)) numiso = numiso + 1
       end do

       allocate(these_iso_masses(numiso), these_iso_abundances(numiso))

       !Get as 1d array the isotopic masses and abundances of element e
       iso_count = 0
       do i = 1, total_numiso
          if(isotope_element(i) .eq. elements(e)) then
             iso_count = iso_count + 1
             these_iso_masses(iso_count) = isotope_mass(i)
             these_iso_abundances(iso_count) = isotope_abundance(i)
          end if
       end do

       !Set the dominant isotopic mass of element e
       m(e) = these_iso_masses(maxloc(these_iso_abundances, dim = 1))
       
       !Calculate g_2 for this element.
       g(e) = 0.0_r64
       do i = 1, numiso
          !This step is an overkill since the perturbation
          !will come out zero for one of the isotopes anyway.
          !Still, explicitly taking out the dominant isotope from
          !the list of defects.
          if(i == maxloc(these_iso_abundances, dim = 1)) cycle
          
          g(e) = g(e) + these_iso_abundances(i)* &
               (1.0_r64 - these_iso_masses(i)/m(e))**2
       end do
       g(e) = g(e)/100.0_r64
       
       deallocate(these_iso_masses, these_iso_abundances)
    end do
  end subroutine calculate_g_DIB

  subroutine create_periodic_table(element, mass, abundance)
      character(len = 3), allocatable, intent(out) :: element(:)
      real(r64), allocatable, intent(out) :: mass(:)
      real(r64), allocatable, intent(out) :: abundance(:)

      allocate(element(287), mass(287), abundance(287)) 
      
      ! Fill in isotope data.
      element(1) = "Ag"
      mass(1) = 106.905095_r64
      abundance(1) = 51.84_r64
      element(2) = "Ag"
      mass(2) = 108.904754_r64
      abundance(2) = 48.16_r64
      element(3) = "Al"
      mass(3) = 26.981541_r64
      abundance(3) = 100.0_r64
      element(4) = "Ar"
      mass(4) = 35.967546_r64
      abundance(4) = 0.34_r64
      element(5) = "Ar"
      mass(5) = 37.962732_r64
      abundance(5) = 0.063_r64
      element(6) = "Ar"
      mass(6) = 39.962383_r64
      abundance(6) = 99.6_r64
      element(7) = "As"
      mass(7) = 74.921596_r64
      abundance(7) = 100.0_r64
      element(8) = "Au"
      mass(8) = 196.96656_r64
      abundance(8) = 100.0_r64
      element(9) = "B"
      mass(9) = 10.012938_r64
      abundance(9) = 19.8_r64
      element(10) = "B"
      mass(10) = 11.009305_r64
      abundance(10) = 80.2_r64
      element(11) = "Ba"
      mass(11) = 129.906277_r64
      abundance(11) = 0.11_r64
      element(12) = "Ba"
      mass(12) = 131.905042_r64
      abundance(12) = 0.1_r64
      element(13) = "Ba"
      mass(13) = 133.90449_r64
      abundance(13) = 2.42_r64
      element(14) = "Ba"
      mass(14) = 134.905668_r64
      abundance(14) = 6.59_r64
      element(15) = "Ba"
      mass(15) = 135.904556_r64
      abundance(15) = 7.85_r64
      element(16) = "Ba"
      mass(16) = 136.905816_r64
      abundance(16) = 11.23_r64
      element(17) = "Ba"
      mass(17) = 137.905236_r64
      abundance(17) = 71.7_r64
      element(18) = "Be"
      mass(18) = 9.012183_r64
      abundance(18) = 100.0_r64
      element(19) = "Bi"
      mass(19) = 208.980388_r64
      abundance(19) = 100.0_r64
      element(20) = "Br"
      mass(20) = 78.918336_r64
      abundance(20) = 50.69_r64
      element(21) = "Br"
      mass(21) = 80.91629_r64
      abundance(21) = 49.31_r64
      element(22) = "C"
      mass(22) = 12.0_r64
      abundance(22) = 98.9_r64
      element(23) = "C"
      mass(23) = 13.003355_r64
      abundance(23) = 1.1_r64
      element(24) = "Ca"
      mass(24) = 39.962591_r64
      abundance(24) = 96.95_r64
      element(25) = "Ca"
      mass(25) = 41.958622_r64
      abundance(25) = 0.65_r64
      element(26) = "Ca"
      mass(26) = 42.95877_r64
      abundance(26) = 0.14_r64
      element(27) = "Ca"
      mass(27) = 43.955485_r64
      abundance(27) = 2.086_r64
      element(28) = "Ca"
      mass(28) = 45.953689_r64
      abundance(28) = 0.004_r64
      element(29) = "Ca"
      mass(29) = 47.952532_r64
      abundance(29) = 0.19_r64
      element(30) = "Cd"
      mass(30) = 105.906461_r64
      abundance(30) = 1.25_r64
      element(31) = "Cd"
      mass(31) = 107.904186_r64
      abundance(31) = 0.89_r64
      element(32) = "Cd"
      mass(32) = 109.903007_r64
      abundance(32) = 12.49_r64
      element(33) = "Cd"
      mass(33) = 110.904182_r64
      abundance(33) = 12.8_r64
      element(34) = "Cd"
      mass(34) = 111.902761_r64
      abundance(34) = 24.13_r64
      element(35) = "Cd"
      mass(35) = 112.904401_r64
      abundance(35) = 12.22_r64
      element(36) = "Cd"
      mass(36) = 113.903361_r64
      abundance(36) = 28.73_r64
      element(37) = "Cd"
      mass(37) = 115.904758_r64
      abundance(37) = 7.49_r64
      element(38) = "Ce"
      mass(38) = 135.90714_r64
      abundance(38) = 0.19_r64
      element(39) = "Ce"
      mass(39) = 137.905996_r64
      abundance(39) = 0.25_r64
      element(40) = "Ce"
      mass(40) = 139.905442_r64
      abundance(40) = 88.48_r64
      element(41) = "Ce"
      mass(41) = 141.909249_r64
      abundance(41) = 11.08_r64
      element(42) = "Cl"
      mass(42) = 34.968853_r64
      abundance(42) = 75.77_r64
      element(43) = "Cl"
      mass(43) = 36.965903_r64
      abundance(43) = 24.23_r64
      element(44) = "Co"
      mass(44) = 58.933198_r64
      abundance(44) = 100.0_r64
      element(45) = "Cr"
      mass(45) = 49.946046_r64
      abundance(45) = 4.35_r64
      element(46) = "Cr"
      mass(46) = 51.94051_r64
      abundance(46) = 83.79_r64
      element(47) = "Cr"
      mass(47) = 52.940651_r64
      abundance(47) = 9.5_r64
      element(48) = "Cr"
      mass(48) = 53.938882_r64
      abundance(48) = 2.36_r64
      element(49) = "Cs"
      mass(49) = 132.905433_r64
      abundance(49) = 100.0_r64
      element(50) = "Cu"
      mass(50) = 62.929599_r64
      abundance(50) = 69.17_r64
      element(51) = "Cu"
      mass(51) = 64.927792_r64
      abundance(51) = 30.83_r64
      element(52) = "Dy"
      mass(52) = 155.924287_r64
      abundance(52) = 0.06_r64
      element(53) = "Dy"
      mass(53) = 157.924412_r64
      abundance(53) = 0.1_r64
      element(54) = "Dy"
      mass(54) = 159.925203_r64
      abundance(54) = 2.34_r64
      element(55) = "Dy"
      mass(55) = 160.926939_r64
      abundance(55) = 18.9_r64
      element(56) = "Dy"
      mass(56) = 161.926805_r64
      abundance(56) = 25.5_r64
      element(57) = "Dy"
      mass(57) = 162.928737_r64
      abundance(57) = 24.9_r64
      element(58) = "Dy"
      mass(58) = 163.929183_r64
      abundance(58) = 28.2_r64
      element(59) = "Er"
      mass(59) = 161.928787_r64
      abundance(59) = 0.14_r64
      element(60) = "Er"
      mass(60) = 163.929211_r64
      abundance(60) = 1.61_r64
      element(61) = "Er"
      mass(61) = 165.930305_r64
      abundance(61) = 33.6_r64
      element(62) = "Er"
      mass(62) = 166.932061_r64
      abundance(62) = 22.95_r64
      element(63) = "Er"
      mass(63) = 167.932383_r64
      abundance(63) = 26.8_r64
      element(64) = "Er"
      mass(64) = 169.935476_r64
      abundance(64) = 14.9_r64
      element(65) = "Eu"
      mass(65) = 150.91986_r64
      abundance(65) = 47.8_r64
      element(66) = "Eu"
      mass(66) = 152.921243_r64
      abundance(66) = 52.2_r64
      element(67) = "F"
      mass(67) = 18.998403_r64
      abundance(67) = 100.0_r64
      element(68) = "Fe"
      mass(68) = 53.939612_r64
      abundance(68) = 5.8_r64
      element(69) = "Fe"
      mass(69) = 55.934939_r64
      abundance(69) = 91.72_r64
      element(70) = "Fe"
      mass(70) = 56.935396_r64
      abundance(70) = 2.2_r64
      element(71) = "Fe"
      mass(71) = 57.933278_r64
      abundance(71) = 0.28_r64
      element(72) = "Ga"
      mass(72) = 68.925581_r64
      abundance(72) = 60.1_r64
      element(73) = "Ga"
      mass(73) = 70.924701_r64
      abundance(73) = 39.9_r64
      element(74) = "Gd"
      mass(74) = 151.919803_r64
      abundance(74) = 0.2_r64
      element(75) = "Gd"
      mass(75) = 153.920876_r64
      abundance(75) = 2.18_r64
      element(76) = "Gd"
      mass(76) = 154.822629_r64
      abundance(76) = 14.8_r64
      element(77) = "Gd"
      mass(77) = 155.92213_r64
      abundance(77) = 20.47_r64
      element(78) = "Gd"
      mass(78) = 156.923967_r64
      abundance(78) = 15.65_r64
      element(79) = "Gd"
      mass(79) = 157.924111_r64
      abundance(79) = 24.84_r64
      element(80) = "Gd"
      mass(80) = 159.927061_r64
      abundance(80) = 21.86_r64
      element(81) = "Ge"
      mass(81) = 69.92425_r64
      abundance(81) = 20.5_r64
      element(82) = "Ge"
      mass(82) = 71.92208_r64
      abundance(82) = 27.4_r64
      element(83) = "Ge"
      mass(83) = 72.923464_r64
      abundance(83) = 7.8_r64
      element(84) = "Ge"
      mass(84) = 73.921179_r64
      abundance(84) = 36.5_r64
      element(85) = "Ge"
      mass(85) = 75.921403_r64
      abundance(85) = 7.8_r64
      element(86) = "H"
      mass(86) = 1.007825_r64
      abundance(86) = 99.99_r64
      element(87) = "H"
      mass(87) = 2.014102_r64
      abundance(87) = 0.015_r64
      element(88) = "He"
      mass(88) = 3.016029_r64
      abundance(88) = 0.0001_r64
      element(89) = "He"
      mass(89) = 4.002603_r64
      abundance(89) = 100.0_r64
      element(90) = "Hf"
      mass(90) = 173.940065_r64
      abundance(90) = 0.16_r64
      element(91) = "Hf"
      mass(91) = 175.94142_r64
      abundance(91) = 5.2_r64
      element(92) = "Hf"
      mass(92) = 176.943233_r64
      abundance(92) = 18.6_r64
      element(93) = "Hf"
      mass(93) = 177.94371_r64
      abundance(93) = 27.1_r64
      element(94) = "Hf"
      mass(94) = 178.945827_r64
      abundance(94) = 13.74_r64
      element(95) = "Hf"
      mass(95) = 179.946561_r64
      abundance(95) = 35.2_r64
      element(96) = "Hg"
      mass(96) = 195.965812_r64
      abundance(96) = 0.15_r64
      element(97) = "Hg"
      mass(97) = 197.96676_r64
      abundance(97) = 10.1_r64
      element(98) = "Hg"
      mass(98) = 198.968269_r64
      abundance(98) = 17.0_r64
      element(99) = "Hg"
      mass(99) = 199.968316_r64
      abundance(99) = 23.1_r64
      element(100) = "Hg"
      mass(100) = 200.970293_r64
      abundance(100) = 13.2_r64
      element(101) = "Hg"
      mass(101) = 201.970632_r64
      abundance(101) = 29.65_r64
      element(102) = "Hg"
      mass(102) = 203.973481_r64
      abundance(102) = 6.8_r64
      element(103) = "Ho"
      mass(103) = 164.930332_r64
      abundance(103) = 100.0_r64
      element(104) = "I"
      mass(104) = 126.904477_r64
      abundance(104) = 100.0_r64
      element(105) = "In"
      mass(105) = 112.904056_r64
      abundance(105) = 4.3_r64
      element(106) = "In"
      mass(106) = 114.903875_r64
      abundance(106) = 95.7_r64
      element(107) = "Ir"
      mass(107) = 190.960603_r64
      abundance(107) = 37.3_r64
      element(108) = "Ir"
      mass(108) = 192.962942_r64
      abundance(108) = 62.7_r64
      element(109) = "K"
      mass(109) = 38.963708_r64
      abundance(109) = 93.2_r64
      element(110) = "K"
      mass(110) = 39.963999_r64
      abundance(110) = 0.012_r64
      element(111) = "K"
      mass(111) = 40.961825_r64
      abundance(111) = 6.73_r64
      element(112) = "Kr"
      mass(112) = 77.920397_r64
      abundance(112) = 0.35_r64
      element(113) = "Kr"
      mass(113) = 79.916375_r64
      abundance(113) = 2.25_r64
      element(114) = "Kr"
      mass(114) = 81.913483_r64
      abundance(114) = 11.6_r64
      element(115) = "Kr"
      mass(115) = 82.914134_r64
      abundance(115) = 11.5_r64
      element(116) = "Kr"
      mass(116) = 83.911506_r64
      abundance(116) = 57.0_r64
      element(117) = "Kr"
      mass(117) = 85.910614_r64
      abundance(117) = 17.3_r64
      element(118) = "La"
      mass(118) = 137.907114_r64
      abundance(118) = 0.09_r64
      element(119) = "La"
      mass(119) = 138.906355_r64
      abundance(119) = 99.91_r64
      element(120) = "Li"
      mass(120) = 6.015123_r64
      abundance(120) = 7.42_r64
      element(121) = "Li"
      mass(121) = 7.016005_r64
      abundance(121) = 92.58_r64
      element(122) = "Lu"
      mass(122) = 174.940785_r64
      abundance(122) = 97.4_r64
      element(123) = "Lu"
      mass(123) = 175.942694_r64
      abundance(123) = 2.6_r64
      element(124) = "Mg"
      mass(124) = 23.985045_r64
      abundance(124) = 78.9_r64
      element(125) = "Mg"
      mass(125) = 24.985839_r64
      abundance(125) = 10.0_r64
      element(126) = "Mg"
      mass(126) = 25.982595_r64
      abundance(126) = 11.1_r64
      element(127) = "Mn"
      mass(127) = 54.938046_r64
      abundance(127) = 100.0_r64
      element(128) = "Mo"
      mass(128) = 91.906809_r64
      abundance(128) = 14.84_r64
      element(129) = "Mo"
      mass(129) = 93.905086_r64
      abundance(129) = 9.25_r64
      element(130) = "Mo"
      mass(130) = 94.905838_r64
      abundance(130) = 15.92_r64
      element(131) = "Mo"
      mass(131) = 95.904676_r64
      abundance(131) = 16.68_r64
      element(132) = "Mo"
      mass(132) = 96.906018_r64
      abundance(132) = 9.55_r64
      element(133) = "Mo"
      mass(133) = 97.905405_r64
      abundance(133) = 24.13_r64
      element(134) = "Mo"
      mass(134) = 99.907473_r64
      abundance(134) = 9.63_r64
      element(135) = "N"
      mass(135) = 14.003074_r64
      abundance(135) = 99.63_r64
      element(136) = "N"
      mass(136) = 15.000109_r64
      abundance(136) = 0.37_r64
      element(137) = "Na"
      mass(137) = 22.98977_r64
      abundance(137) = 100.0_r64
      element(138) = "Nb"
      mass(138) = 92.906378_r64
      abundance(138) = 100.0_r64
      element(139) = "Nd"
      mass(139) = 141.907731_r64
      abundance(139) = 27.13_r64
      element(140) = "Nd"
      mass(140) = 142.909823_r64
      abundance(140) = 12.18_r64
      element(141) = "Nd"
      mass(141) = 143.910096_r64
      abundance(141) = 23.8_r64
      element(142) = "Nd"
      mass(142) = 144.912582_r64
      abundance(142) = 8.3_r64
      element(143) = "Nd"
      mass(143) = 145.913126_r64
      abundance(143) = 17.19_r64
      element(144) = "Nd"
      mass(144) = 147.916901_r64
      abundance(144) = 5.76_r64
      element(145) = "Nd"
      mass(145) = 149.9209_r64
      abundance(145) = 5.64_r64
      element(146) = "Ne"
      mass(146) = 19.992439_r64
      abundance(146) = 90.6_r64
      element(147) = "Ne"
      mass(147) = 20.993845_r64
      abundance(147) = 0.26_r64
      element(148) = "Ne"
      mass(148) = 21.991384_r64
      abundance(148) = 9.2_r64
      element(149) = "Ni"
      mass(149) = 57.935347_r64
      abundance(149) = 68.27_r64
      element(150) = "Ni"
      mass(150) = 59.930789_r64
      abundance(150) = 26.1_r64
      element(151) = "Ni"
      mass(151) = 60.931059_r64
      abundance(151) = 1.13_r64
      element(152) = "Ni"
      mass(152) = 61.928346_r64
      abundance(152) = 3.59_r64
      element(153) = "Ni"
      mass(153) = 63.927968_r64
      abundance(153) = 0.91_r64
      element(154) = "O"
      mass(154) = 15.994915_r64
      abundance(154) = 99.76_r64
      element(155) = "O"
      mass(155) = 16.999131_r64
      abundance(155) = 0.038_r64
      element(156) = "O"
      mass(156) = 17.999159_r64
      abundance(156) = 0.2_r64
      element(157) = "Os"
      mass(157) = 183.952514_r64
      abundance(157) = 0.02_r64
      element(158) = "Os"
      mass(158) = 185.953852_r64
      abundance(158) = 1.58_r64
      element(159) = "Os"
      mass(159) = 186.955762_r64
      abundance(159) = 1.6_r64
      element(160) = "Os"
      mass(160) = 187.95585_r64
      abundance(160) = 13.3_r64
      element(161) = "Os"
      mass(161) = 188.958156_r64
      abundance(161) = 16.1_r64
      element(162) = "Os"
      mass(162) = 189.958455_r64
      abundance(162) = 26.4_r64
      element(163) = "Os"
      mass(163) = 191.961487_r64
      abundance(163) = 41.0_r64
      element(164) = "P"
      mass(164) = 30.973763_r64
      abundance(164) = 100.0_r64
      element(165) = "Pb"
      mass(165) = 203.973037_r64
      abundance(165) = 1.4_r64
      element(166) = "Pb"
      mass(166) = 205.974455_r64
      abundance(166) = 24.1_r64
      element(167) = "Pb"
      mass(167) = 206.975885_r64
      abundance(167) = 22.1_r64
      element(168) = "Pb"
      mass(168) = 207.976641_r64
      abundance(168) = 52.4_r64
      element(169) = "Pd"
      mass(169) = 101.905609_r64
      abundance(169) = 1.02_r64
      element(170) = "Pd"
      mass(170) = 103.904026_r64
      abundance(170) = 11.14_r64
      element(171) = "Pd"
      mass(171) = 104.905075_r64
      abundance(171) = 22.33_r64
      element(172) = "Pd"
      mass(172) = 105.903475_r64
      abundance(172) = 27.33_r64
      element(173) = "Pd"
      mass(173) = 107.903894_r64
      abundance(173) = 26.46_r64
      element(174) = "Pd"
      mass(174) = 109.905169_r64
      abundance(174) = 11.72_r64
      element(175) = "Pr"
      mass(175) = 140.907657_r64
      abundance(175) = 100.0_r64
      element(176) = "Pt"
      mass(176) = 189.959937_r64
      abundance(176) = 0.01_r64
      element(177) = "Pt"
      mass(177) = 191.961049_r64
      abundance(177) = 0.79_r64
      element(178) = "Pt"
      mass(178) = 193.962679_r64
      abundance(178) = 32.9_r64
      element(179) = "Pt"
      mass(179) = 194.964785_r64
      abundance(179) = 33.8_r64
      element(180) = "Pt"
      mass(180) = 195.964947_r64
      abundance(180) = 25.3_r64
      element(181) = "Pt"
      mass(181) = 197.967879_r64
      abundance(181) = 7.2_r64
      element(182) = "Rb"
      mass(182) = 84.9118_r64
      abundance(182) = 72.17_r64
      element(183) = "Rb"
      mass(183) = 86.909184_r64
      abundance(183) = 27.84_r64
      element(184) = "Re"
      mass(184) = 184.952977_r64
      abundance(184) = 37.4_r64
      element(185) = "Re"
      mass(185) = 186.955765_r64
      abundance(185) = 62.6_r64
      element(186) = "Rh"
      mass(186) = 102.905503_r64
      abundance(186) = 100.0_r64
      element(187) = "Ru"
      mass(187) = 95.907596_r64
      abundance(187) = 5.52_r64
      element(188) = "Ru"
      mass(188) = 97.905287_r64
      abundance(188) = 1.88_r64
      element(189) = "Ru"
      mass(189) = 98.905937_r64
      abundance(189) = 12.7_r64
      element(190) = "Ru"
      mass(190) = 99.904218_r64
      abundance(190) = 12.6_r64
      element(191) = "Ru"
      mass(191) = 100.905581_r64
      abundance(191) = 17.0_r64
      element(192) = "Ru"
      mass(192) = 101.904348_r64
      abundance(192) = 31.6_r64
      element(193) = "Ru"
      mass(193) = 103.905422_r64
      abundance(193) = 18.7_r64
      element(194) = "S"
      mass(194) = 31.972072_r64
      abundance(194) = 95.02_r64
      element(195) = "S"
      mass(195) = 32.971459_r64
      abundance(195) = 0.75_r64
      element(196) = "S"
      mass(196) = 33.967868_r64
      abundance(196) = 4.21_r64
      element(197) = "S"
      mass(197) = 35.967079_r64
      abundance(197) = 0.02_r64
      element(198) = "Sb"
      mass(198) = 120.903824_r64
      abundance(198) = 57.3_r64
      element(199) = "Sb"
      mass(199) = 122.904222_r64
      abundance(199) = 42.7_r64
      element(200) = "Sc"
      mass(200) = 44.955914_r64
      abundance(200) = 100.0_r64
      element(201) = "Se"
      mass(201) = 73.922477_r64
      abundance(201) = 0.9_r64
      element(202) = "Se"
      mass(202) = 75.919207_r64
      abundance(202) = 9.0_r64
      element(203) = "Se"
      mass(203) = 76.919908_r64
      abundance(203) = 7.6_r64
      element(204) = "Se"
      mass(204) = 77.917304_r64
      abundance(204) = 23.5_r64
      element(205) = "Se"
      mass(205) = 79.916521_r64
      abundance(205) = 49.6_r64
      element(206) = "Se"
      mass(206) = 81.916709_r64
      abundance(206) = 9.4_r64
      element(207) = "Si"
      mass(207) = 27.976928_r64
      abundance(207) = 92.23_r64
      element(208) = "Si"
      mass(208) = 28.976496_r64
      abundance(208) = 4.67_r64
      element(209) = "Si"
      mass(209) = 29.973772_r64
      abundance(209) = 3.1_r64
      element(210) = "Sm"
      mass(210) = 143.912009_r64
      abundance(210) = 3.1_r64
      element(211) = "Sm"
      mass(211) = 146.914907_r64
      abundance(211) = 15.0_r64
      element(212) = "Sm"
      mass(212) = 147.914832_r64
      abundance(212) = 11.3_r64
      element(213) = "Sm"
      mass(213) = 148.917193_r64
      abundance(213) = 13.8_r64
      element(214) = "Sm"
      mass(214) = 149.917285_r64
      abundance(214) = 7.4_r64
      element(215) = "Sm"
      mass(215) = 151.919741_r64
      abundance(215) = 26.7_r64
      element(216) = "Sm"
      mass(216) = 153.922218_r64
      abundance(216) = 22.7_r64
      element(217) = "Sn"
      mass(217) = 111.904826_r64
      abundance(217) = 0.97_r64
      element(218) = "Sn"
      mass(218) = 113.902784_r64
      abundance(218) = 0.65_r64
      element(219) = "Sn"
      mass(219) = 114.903348_r64
      abundance(219) = 0.36_r64
      element(220) = "Sn"
      mass(220) = 115.901744_r64
      abundance(220) = 14.7_r64
      element(221) = "Sn"
      mass(221) = 116.902954_r64
      abundance(221) = 7.7_r64
      element(222) = "Sn"
      mass(222) = 117.901607_r64
      abundance(222) = 24.3_r64
      element(223) = "Sn"
      mass(223) = 118.90331_r64
      abundance(223) = 8.6_r64
      element(224) = "Sn"
      mass(224) = 119.902199_r64
      abundance(224) = 32.4_r64
      element(225) = "Sn"
      mass(225) = 121.90344_r64
      abundance(225) = 4.6_r64
      element(226) = "Sn"
      mass(226) = 123.905271_r64
      abundance(226) = 5.6_r64
      element(227) = "Sr"
      mass(227) = 83.913428_r64
      abundance(227) = 0.56_r64
      element(228) = "Sr"
      mass(228) = 85.909273_r64
      abundance(228) = 9.86_r64
      element(229) = "Sr"
      mass(229) = 86.908902_r64
      abundance(229) = 7.0_r64
      element(230) = "Sr"
      mass(230) = 87.905625_r64
      abundance(230) = 82.58_r64
      element(231) = "Ta"
      mass(231) = 179.947489_r64
      abundance(231) = 0.012_r64
      element(232) = "Ta"
      mass(232) = 180.948014_r64
      abundance(232) = 99.99_r64
      element(233) = "Tb"
      mass(233) = 158.92535_r64
      abundance(233) = 100.0_r64
      element(234) = "Te"
      mass(234) = 119.904021_r64
      abundance(234) = 0.096_r64
      element(235) = "Te"
      mass(235) = 121.903055_r64
      abundance(235) = 2.6_r64
      element(236) = "Te"
      mass(236) = 122.904278_r64
      abundance(236) = 0.91_r64
      element(237) = "Te"
      mass(237) = 123.902825_r64
      abundance(237) = 4.82_r64
      element(238) = "Te"
      mass(238) = 124.904435_r64
      abundance(238) = 7.14_r64
      element(239) = "Te"
      mass(239) = 125.90331_r64
      abundance(239) = 18.95_r64
      element(240) = "Te"
      mass(240) = 127.904464_r64
      abundance(240) = 31.69_r64
      element(241) = "Te"
      mass(241) = 129.906229_r64
      abundance(241) = 33.8_r64
      element(242) = "Th"
      mass(242) = 232.038054_r64
      abundance(242) = 100.0_r64
      element(243) = "Ti"
      mass(243) = 45.952633_r64
      abundance(243) = 8.0_r64
      element(244) = "Ti"
      mass(244) = 46.951765_r64
      abundance(244) = 7.3_r64
      element(245) = "Ti"
      mass(245) = 47.947947_r64
      abundance(245) = 73.8_r64
      element(246) = "Ti"
      mass(246) = 48.947871_r64
      abundance(246) = 5.5_r64
      element(247) = "Ti"
      mass(247) = 49.944786_r64
      abundance(247) = 5.4_r64
      element(248) = "Tl"
      mass(248) = 202.972336_r64
      abundance(248) = 29.52_r64
      element(249) = "Tl"
      mass(249) = 204.97441_r64
      abundance(249) = 70.48_r64
      element(250) = "Tm"
      mass(250) = 168.934225_r64
      abundance(250) = 100.0_r64
      element(251) = "U"
      mass(251) = 234.040947_r64
      abundance(251) = 0.006_r64
      element(252) = "U"
      mass(252) = 235.043925_r64
      abundance(252) = 0.72_r64
      element(253) = "U"
      mass(253) = 238.050786_r64
      abundance(253) = 99.27_r64
      element(254) = "V"
      mass(254) = 49.947161_r64
      abundance(254) = 0.25_r64
      element(255) = "V"
      mass(255) = 50.943963_r64
      abundance(255) = 99.75_r64
      element(256) = "W"
      mass(256) = 179.946727_r64
      abundance(256) = 0.13_r64
      element(257) = "W"
      mass(257) = 181.948225_r64
      abundance(257) = 26.3_r64
      element(258) = "W"
      mass(258) = 182.950245_r64
      abundance(258) = 14.3_r64
      element(259) = "W"
      mass(259) = 183.950953_r64
      abundance(259) = 30.67_r64
      element(260) = "W"
      mass(260) = 185.954377_r64
      abundance(260) = 28.6_r64
      element(261) = "Xe"
      mass(261) = 123.905894_r64
      abundance(261) = 0.1_r64
      element(262) = "Xe"
      mass(262) = 125.904281_r64
      abundance(262) = 0.09_r64
      element(263) = "Xe"
      mass(263) = 127.903531_r64
      abundance(263) = 1.91_r64
      element(264) = "Xe"
      mass(264) = 128.90478_r64
      abundance(264) = 26.4_r64
      element(265) = "Xe"
      mass(265) = 129.90351_r64
      abundance(265) = 4.1_r64
      element(266) = "Xe"
      mass(266) = 130.905076_r64
      abundance(266) = 21.2_r64
      element(267) = "Xe"
      mass(267) = 131.904148_r64
      abundance(267) = 26.9_r64
      element(268) = "Xe"
      mass(268) = 133.905395_r64
      abundance(268) = 10.4_r64
      element(269) = "Xe"
      mass(269) = 135.907219_r64
      abundance(269) = 8.9_r64
      element(270) = "Y"
      mass(270) = 88.905856_r64
      abundance(270) = 100.0_r64
      element(271) = "Yb"
      mass(271) = 167.933908_r64
      abundance(271) = 0.13_r64
      element(272) = "Yb"
      mass(272) = 169.934774_r64
      abundance(272) = 3.05_r64
      element(273) = "Yb"
      mass(273) = 170.936338_r64
      abundance(273) = 14.3_r64
      element(274) = "Yb"
      mass(274) = 171.936393_r64
      abundance(274) = 21.9_r64
      element(275) = "Yb"
      mass(275) = 172.938222_r64
      abundance(275) = 16.12_r64
      element(276) = "Yb"
      mass(276) = 173.938873_r64
      abundance(276) = 31.8_r64
      element(277) = "Yb"
      mass(277) = 175.942576_r64
      abundance(277) = 12.7_r64
      element(278) = "Zn"
      mass(278) = 63.929145_r64
      abundance(278) = 48.6_r64
      element(279) = "Zn"
      mass(279) = 65.926035_r64
      abundance(279) = 27.9_r64
      element(280) = "Zn"
      mass(280) = 66.927129_r64
      abundance(280) = 4.1_r64
      element(281) = "Zn"
      mass(281) = 67.924846_r64
      abundance(281) = 18.8_r64
      element(282) = "Zn"
      mass(282) = 69.925325_r64
      abundance(282) = 0.6_r64
      element(283) = "Zr"
      mass(283) = 89.904708_r64
      abundance(283) = 51.45_r64
      element(284) = "Zr"
      mass(284) = 90.905644_r64
      abundance(284) = 11.27_r64
      element(285) = "Zr"
      mass(285) = 91.905039_r64
      abundance(285) = 17.17_r64
      element(286) = "Zr"
      mass(286) = 93.906319_r64
      abundance(286) = 17.33_r64
      element(287) = "Zr"
      mass(287) = 95.908272_r64
      abundance(287) = 2.78_r64
  end subroutine create_periodic_table
end module crystal_module

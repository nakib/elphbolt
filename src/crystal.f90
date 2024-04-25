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

  use precision, only: r64, i64
  use params, only: twopi
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
     integer :: dim
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
       self%dim = 2
    else
       self%dim = 3
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
    use params, only: lookup_periodic_table 
    use isotopes_module, only: isotopes
    
    character(len=3), intent(in) :: elements(:)
    real(r64), intent(out) :: m(:), g(:)

    !Local variables
    integer(i64) :: i, nelems, e
    type(isotopes) :: iso_data

    nelems = size(elements)
        
    do e = 1, nelems
       iso_data = lookup_periodic_table(elements(e))
       m(e) = 0.0_r64
       do i = 1, iso_data%numisotopes
             m(e) = m(e) + iso_data%masses(i)*iso_data%abundances(i)
       end do
       m(e) = m(e)/100.0_r64
       g(e) = 0.0_r64
       do i = 1, iso_data%numisotopes
             g(e) = g(e) + iso_data%abundances(i)*(1.0_r64 - iso_data%masses(i)/m(e))**2
       end do
       g(e) = g(e)/100.0_r64
    end do
  end subroutine calculate_mavg_and_g

  subroutine calculate_g_DIB(elements, m, g)
    !! Find the dominant isotopic masses from the periodic table
    !! and the average mass perturbation for use in the DIB-1st Born ph-iso scattering theory.
    use params, only: lookup_periodic_table
    use isotopes_module, only: isotopes
    
    character(len=3), intent(in) :: elements(:)
    real(r64), intent(out) :: m(:), g(:)

    !Local variables
    integer :: i, nelems, e
    type(isotopes) :: iso_data 

    nelems = size(elements)

    do e = 1, nelems
       iso_data = lookup_periodic_table(elements(e))

       !Set the dominant isotopic mass of element e
       m(e) = iso_data%masses(maxloc(iso_data%abundances, dim = 1))
       
       !Calculate g_2 for this element.
       g(e) = 0.0_r64
       do i = 1, iso_data%numisotopes
          !This step is an overkill since the perturbation
          !will come out zero for one of the isotopes anyway.
          !Still, explicitly taking out the dominant isotope from
          !the list of defects.
          if(i == maxloc(iso_data%abundances, dim = 1)) cycle
          
          g(e) = g(e) + iso_data%abundances(i)* &
               (1.0_r64 - iso_data%masses(i)/m(e))**2
       end do
       g(e) = g(e)/100.0_r64
       
    end do
  end subroutine calculate_g_DIB

end module crystal_module

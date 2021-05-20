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
     !! Elements in the in the basis.
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
    real(dp), allocatable :: masses(:), born(:,:,:), basis(:,:), basis_cart(:,:)
    real(dp) :: epsilon(3,3), lattvecs(3,3), volume, reclattvecs(3,3), volume_bz, T
    character(len=3), allocatable :: elements(:)
    character(len=100) :: name
    logical :: polar
    
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ name, elements, atomtypes, basis, lattvecs, &
         polar, born, epsilon, masses, T

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
         masses(numatoms), basis(3,numatoms), basis_cart(3,numatoms))
    allocate(c%elements(c%numelements), c%atomtypes(c%numatoms), c%born(3,3,c%numatoms), &
         c%masses(c%numatoms), c%basis(3,c%numatoms), c%basis_cart(3,c%numatoms))
    
    !Read crystal_info
    c%polar = .false.
    c%epsilon = 0.0_dp
    c%born = 0.0_dp
    c%T = -1.0_dp
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
    c%masses = masses
    c%basis = basis
    c%polar = polar
    c%lattvecs = lattvecs
    c%T = T

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
       print*, 'Material: ', c%name
       print*, 'Lattice vectors [nm]:'
       print*, c%lattvecs(:,1)
       print*, c%lattvecs(:,2)
       print*, c%lattvecs(:,3)
       print*, 'Primitive cell volume =', c%volume, 'nm^3'

       print*, 'Reciprocal lattice vectors [1/nm]:'
       print*, c%reclattvecs(:,1)
       print*, c%reclattvecs(:,2)
       print*, c%reclattvecs(:,3)
       print*, 'Brillouin zone volume =', c%volume_bz, '1/nm^3'

       if(c%polar) then
          print*, 'System is polar.'
          print*, 'Dielectric tensor:'
          print*, c%epsilon(:,1)
          print*, c%epsilon(:,2)
          print*, c%epsilon(:,3)
       end if
       
       print*, 'Crystal temperature =', c%T, 'K'
    end if
  end subroutine read_input_and_setup_crystal

  subroutine calculate_wavevectors_full(mesh,wavevecs,blocks,indexlist)
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
end module crystal_module

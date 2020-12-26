module phonon_module
  !! Module containing type and procedures related to the phononic properties.

  use params, only: dp, k4
  use misc, only: print_message
  use numerics_module, only: numerics
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal, calculate_wavevectors_full
  use symmetry_module, only: symmetry, find_irred_wedge, create_fbz2ibz_map
  
  implicit none

  private
  public phonon
  
  type phonon
     !! Data and procedures related to phonons.

     integer(k4) :: numbranches
     !! Total number of phonon branches.
     integer(k4) :: nq
     !! Number of phonon wave vectors in the full Brillouin zone (FBZ).
     integer(k4) :: nq_irred
     !! Number of phonon wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     integer(k4) :: qmesh(3) 
     !! Phonon wave vector mesh.
     real(dp), allocatable :: wavevecs(:,:)
     !! List of all phonon wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible phonon wave vectors (crystal coordinates).
     integer(k4), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(k4), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k4), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k4), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ phonon point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k4), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ phonon point to its IBZ wedge image.
     real(dp), allocatable :: ens(:,:)
     !! List of phonon energies on FBZ.
     real(dp), allocatable :: ens_irred(:,:)
     !! List of phonon energies on IBZ.
     real(dp), allocatable :: vels(:,:,:)
     !! List of phonon velocities on FBZ.
     real(dp), allocatable :: vels_irred(:,:,:)
     !! List of phonon velocites on IBZ.
     complex(dp), allocatable :: evecs(:,:,:)
     !! List of all phonon eigenvectors.
     complex(dp), allocatable :: evecs_irred(:,:,:)
     !! List of IBZ wedge phonon eigenvectors.

   contains

     procedure :: initialize
     
  end type phonon

contains

  subroutine initialize(ph, wann, crys, sym, num)
    !! Initialize the phonon data type and calculate ground state phonon properties.

    class(phonon), intent(out) :: ph
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num

    !Set phonon branches
    ph%numbranches = crys%numatoms*3
    !Set wave vector mesh
    ph%qmesh = num%qmesh
    !Set number of phonon wave vectors
    ph%nq = product(ph%qmesh(:))

    call calculate_phonons(ph, wann, crys, sym, num)
    
  end subroutine initialize
  
  subroutine calculate_phonons(ph, wann, crys, sym, num)
    !! Calculate phonon quantities on the FBZ and IBZ meshes.

    class(phonon), intent(inout) :: ph
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(k4) :: iq
    !Switch for mesh utilites with or without energy restriction
    logical :: blocks

    blocks = .false.

    call print_message("Calculating phonon FBZ quantities...")

    allocate(ph%indexlist(ph%nq))
    do iq = 1, ph%nq
       ph%indexlist(iq) = iq
    end do

    !Calculate FBZ mesh
    call calculate_wavevectors_full(ph%qmesh, ph%wavevecs, blocks)

    !Calculate FBZ phonon quantities
    allocate(ph%ens(ph%nq, ph%numbranches))
    allocate(ph%vels(ph%nq, ph%numbranches, 3))
    allocate(ph%evecs(ph%nq, ph%numbranches, ph%numbranches))    
    call wann%ph_wann_epw(crys, ph%nq, ph%wavevecs, ph%ens, ph%vels, ph%evecs)

    !Calculate IBZ mesh
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(ph%qmesh, ph%nq_irred, ph%wavevecs_irred, &
         ph%indexlist_irred, ph%nequiv, sym%nsymm_rot, sym%qrotations, ph%ibz2fbz_map, blocks)

    !Create FBZ to IBZ map
    call print_message("Calculating FBZ -> IBZ mappings...")
    call create_fbz2ibz_map(ph%fbz2ibz_map, ph%nq, ph%nq_irred, ph%indexlist, &
         ph%nequiv, ph%ibz2fbz_map)
  end subroutine calculate_phonons
end module phonon_module

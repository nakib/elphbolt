module derived_types
  !! Module containing the derived data types.

  use params, only: dp, k4, twopi
  
  implicit none

  public

  type crystal_data
     !! Container for crystal related data.

     integer(k4) :: numelements
     !! Number of types of basis atoms.
     integer(k4) :: numatoms
     !! Number of basis atoms.
     character(len=3), allocatable :: elements(:)
     !! Elements in the in the basis.
     integer(k4), allocatable :: atomtypes(:)
     !! Integer tagging unique elements in the basis.
     real(dp), allocatable :: masses(:)
     !! Masses of the basis atoms.
     real(dp), allocatable :: basis(:,:)
     !! Basis vectors (crystal coordinates).
     real(dp), allocatable :: basis_cart(:,:)
     !! Basis vectors (Cartesian coordinates).
     real(dp) :: lattvecs(3,3)
     !! Lattice vectors (nm).
     real(dp) :: volume
     !! Volume of primitive cell (nm^3).
  end type crystal_data

  type reciprocal_lattice_data
     !! Container for reciprocal lattice related data.

     real(dp) :: reclattvecs(3,3)
     !! Reciprocal lattice vectors.
     real(dp) :: volume_bz
     !! Brillouin zone volume (nm^-3).
     integer(k4) :: qmesh_coarse(3) 
     !! Coarse phonon wave vector mesh.
     integer(k4) :: qmesh_fine(3) 
     !! Fine phonon wave vector mesh.
     integer(k4) :: kmesh_coarse(3)
     !! Electron wave vector mesh.
     integer(k4) :: kmesh_fine(3)
     !! Fine electron wave vector mesh.
  end type reciprocal_lattice_data

  type symmetry_data
     !! Container for symmetry related data.
     
     integer(k4) :: nsymm
     !! Number of spacegroup symmetries.
     integer(k4) :: nsymm_rot
     !! Number of rotations.
     !Rotation matrices without time-reversal:
     integer(k4), allocatable :: rotations_orig(:,:,:) !Real space, crystal coordinates
     real(dp), allocatable :: crotations_orig(:,:,:) !Real space, cartesian coordinates
     real(dp), allocatable :: qrotations_orig(:,:,:) !Reciprocal space, crystal coordinates
     !And with time-reversal (time reversed sector is the 2nd half of last axis):
     integer(k4), allocatable :: rotations(:,:,:) 
     !! Rotation matrices: real space, crystal coordinates.
     real(dp), allocatable :: crotations(:,:,:)
     !! Rotation matrices: real space, cartesian coordinates.
     real(dp), allocatable :: qrotations(:,:,:)
     !! Rotation matrices:: reciprocal space, crystal coordinates.
     integer(k4), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors (full Brillouin zone).
     integer(k4), allocatable :: equiv_map_blocks(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors (energy windowed blocks of full Brillouin zone).
     character(len=10) :: international
     !! Spacegroup in Hermann–Mauguin (or international) notation.
  end type symmetry_data

  type electron_data
     !! Container for electronic data.
     
     integer(k4) :: numbands
     !! Total number of electronic bands.
     integer(kind = 4) :: coarse_mesh(3)
     !! Coarse wave vector mesh.
     integer(kind = 4) :: fine_mesh(3)
     !! Fine wave vector mesh.
     integer(k4) :: nq
     !! Number of electron wave vectors in the full Brillouin zone (FBZ).
     integer(k4) :: nq_irred
     !! Number of electron wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     integer(k4) :: nstates_inwindow
     !! Number of electron wave vectors within transport window.
     integer(k4) :: nstates_irred_inwindow
     !! Number of IBZ wedge electron wave vectors within transport window.
     integer(k4), allocatable :: IBZ_inwindow_states(:,:)
     !! List of irreducible wedge states within transport window.
     real(dp), allocatable :: wavevecs(:, :)
     !! List of all electron wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:, :)
     !! List of irreducible electron wave vectors (crystal coordinates).
     integer(k4), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(k4), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k4), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k4), allocatable :: ibz2fbz_map(:, :, :)
     !! Map from an IBZ electron point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k4), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ electron point to its IBZ wedge image.
     real(dp), allocatable :: ens(:, :)
     !! List of electron energies on FBZ.
     real(dp), allocatable :: ens_irred(:, :)
     !! List of electron energies on IBZ.
     real(dp), allocatable :: vels(:, :, :)
     !! List of electron velocities on FBZ
     real(dp), allocatable :: vels_irred(:, :, :)
     !! List of electron velicites on IBZ.
     complex(dp), allocatable :: evecs(:, :, :)
     !! List of all electron eigenvectors.
     complex(dp), allocatable :: evecs_irred(:, :, :)
     !! List of IBZ wedge electron eigenvectors.
  end type electron_data

  !type phonon_data
     !! Container for electronic data.

     
  !end type phonon_data

  !type directories
     !! Container for i/o directories
  !end type directories
end module derived_types



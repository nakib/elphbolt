module derived_types
  !! Module containing the derived data types.

  use params, only: dp, k4
  
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
     integer(k4), allocatable :: rotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations_orig(:,:,:)
     !! Rotations without time-reversal, reciprocal space, crystal coordinates.
     !And with time-reversal (time reversed sector is the 2nd half of last axis):
     integer(k4), allocatable :: rotations(:,:,:) 
     !! Rotations with time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations(:,:,:)
     !! Rotations with time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations(:,:,:)
     !! Rotations with time-reversal, reciprocal space, crystal coordinates.
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
     integer(k4) :: nk
     !! Number of fine electron wave vectors in the full Brillouin zone (FBZ).
     integer(k4) :: nk_irred
     !! Number of fine electron wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     integer(k4) :: nstates_inwindow
     !! Number of electron wave vectors within transport window.
     integer(k4) :: nstates_irred_inwindow
     !! Number of IBZ wedge electron wave vectors within transport window.
     integer(k4), allocatable :: IBZ_inwindow_states(:,:)
     !! List of irreducible wedge states within transport window.
     real(dp), allocatable :: wavevecs(:,:)
     !! List of all electron wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible electron wave vectors (crystal coordinates).
     integer(k4), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(k4), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k4), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k4), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ electron point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k4), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ electron point to its IBZ wedge image.
     real(dp), allocatable :: ens(:,:)
     !! List of electron energies on FBZ.
     real(dp), allocatable :: ens_irred(:,:)
     !! List of electron energies on IBZ.
     real(dp), allocatable :: vels(:,:,:)
     !! List of electron velocities on FBZ.
     real(dp), allocatable :: vels_irred(:,:,:)
     !! List of electron velocites on IBZ.
     complex(dp), allocatable :: evecs(:,:,:)
     !! List of all electron eigenvectors.
     complex(dp), allocatable :: evecs_irred(:,:,:)
     !! List of IBZ wedge electron eigenvectors.
  end type electron_data

  type phonon_data
     !! Container for phononic data.

     integer(k4) :: numbranches
     !! Total number of phonon branches.
     integer(kind = 4) :: mesh(3)
     !! Wave vector mesh.
     integer(k4) :: nq
     !! Number of phonon wave vectors in the full Brillouin zone (FBZ).
     integer(k4) :: nq_irred
     !! Number of phonon wave vectors in the irreducible wedge of Brillouin zone (IBZ).
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
  end type phonon_data

  type EPW_data
     !! Container for the EPW Wannier representation data.

     integer(k4) :: numwannbands
     !! Number of Wannier bands.
     integer(k4) :: numbranches
     !! Number of phonon branches.
     integer(k4) :: nwsk
     !! Number of real space cells for electrons.
     integer(k4) :: nwsq
     !! Number of real space cells for phonons.
     integer(k4) :: nwsg
     !! Number of real space cells for electron-phonon vertex.
     integer(k4), allocatable :: rcells_k(:, :)
     !! Real space cell locations for electrons.
     integer(k4), allocatable :: rcells_q(:, :)
     !! Real space cell locations for phonons.
     integer(k4), allocatable :: rcells_g(:, :)
     !! Real space cell locations for electron-phonon vertex.     
     integer(k4), allocatable :: elwsdeg(:)
     !! Real space cell multiplicity for electrons.
     integer(k4), allocatable :: phwsdeg(:)
     !! Real space cell multiplicity for phonons.
     integer(k4), allocatable :: gwsdeg(:)
     !! Real space cell multiplicity for electron-phonon vertex.
     complex(dp), allocatable :: Hwann(:, :, :)
     !! Hamiltonian in Wannier representation.
     complex(dp), allocatable :: Dphwann(:, :, :)
     !! Dynamical matrix in Wannier representation.
     complex(dp), allocatable :: gwann(:, :, :, :, :)
     !! e-ph vertex in Wannier representation.
  end type EPW_data

  type control_data
     !! Container for i/o directories and system control flags.

     character(:), allocatable :: cwd
     !! Current working directory.
     character(:), allocatable :: datadumpdir
     !! Runtime data dump repository.
     character(:), allocatable :: g2dir
     !! Directory for e-ph vertex.
     character(:), allocatable :: Vdir
     !! Directory for ph-ph vertex.

     logical :: read_g2
     !! Choose if earlier e-ph vertex is to be used.
     logical :: read_V
     !! Choose if earlier p-ph vertex is to be used.
  end type control_data
end module derived_types



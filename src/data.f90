module data
  !! Module containing the derived data types.

  use params, only: dp, k4
  
  implicit none

  public

  !Crystal related data:
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

  !Reciprocal lattice related data:
  real(dp) :: reclattvecs(3,3)
  !! Reciprocal lattice vectors.
  real(dp) :: volume_bz
  !! Brillouin zone volume (nm^-3).

  !Symmetry related data.
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

  !Electronic data.
  integer(k4) :: numbands
  !! Total number of electronic bands.
  integer(k4) :: mesh_ref
  !! Electron mesh refinement factor compared to the phonon mesh.
  integer(k4) :: kmesh(3)
  !! Electron wave vector mesh.
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
  real(dp) :: enref
  !! Electron reference energy (eV).
  real(dp) :: fsthick
  !! Fermi surface thickness in (eV).
  real(dp), allocatable :: el_wavevecs(:,:)
  !! List of all electron wave vectors (crystal coordinates).
  real(dp), allocatable :: el_wavevecs_irred(:,:)
  !! List of irreducible electron wave vectors (crystal coordinates).
  integer(k4), allocatable :: el_indexlist(:)
  !! List of muxed indices of the FBZ wave vectors.
  integer(k4), allocatable :: el_indexlist_irred(:)
  !! List of muxed indices of the IBZ wedge.
  integer(k4), allocatable :: el_nequiv(:)
  !! List of the number of equivalent points for each IBZ point.
  integer(k4), allocatable :: el_ibz2fbz_map(:,:,:)
  !! Map from an IBZ electron point to its images.
  !! The third axis contains the pair (symmetry index, image).
  integer(k4), allocatable :: el_fbz2ibz_map(:)
  !! Map from an FBZ electron point to its IBZ wedge image.
  real(dp), allocatable :: el_ens(:,:)
  !! List of electron energies on FBZ.
  real(dp), allocatable :: el_ens_irred(:,:)
  !! List of electron energies on IBZ.
  real(dp), allocatable :: el_vels(:,:,:)
  !! List of electron velocities on FBZ.
  real(dp), allocatable :: el_vels_irred(:,:,:)
  !! List of electron velocites on IBZ.
  complex(dp), allocatable :: el_evecs(:,:,:)
  !! List of all electron eigenvectors.
  complex(dp), allocatable :: el_evecs_irred(:,:,:)
  !! List of IBZ wedge electron eigenvectors.

  !Phononic data.
  integer(k4) :: numbranches
  !! Total number of phonon branches.
  integer(k4) :: nq
  !! Number of phonon wave vectors in the full Brillouin zone (FBZ).
  integer(k4) :: nq_irred
  !! Number of phonon wave vectors in the irreducible wedge of Brillouin zone (IBZ).
  integer(k4) :: qmesh(3) 
  !! Coarse phonon wave vector mesh.
  integer(k4) :: qmesh_fine(3) 
  !! Fine phonon wave vector mesh.
  real(dp), allocatable :: ph_wavevecs(:,:)
  !! List of all phonon wave vectors (crystal coordinates).
  real(dp), allocatable :: ph_wavevecs_irred(:,:)
  !! List of irreducible phonon wave vectors (crystal coordinates).
  integer(k4), allocatable :: ph_indexlist(:)
  !! List of muxed indices of the FBZ wave vectors.
  integer(k4), allocatable :: ph_indexlist_irred(:)
  !! List of muxed indices of the IBZ wedge.
  integer(k4), allocatable :: ph_nequiv(:)
  !! List of the number of equivalent points for each IBZ point.
  integer(k4), allocatable :: ph_ibz2fbz_map(:,:,:)
  !! Map from an IBZ phonon point to its images.
  !! The third axis contains the pair (symmetry index, image).
  integer(k4), allocatable :: ph_fbz2ibz_map(:)
  !! Map from an FBZ phonon point to its IBZ wedge image.
  real(dp), allocatable :: ph_ens(:,:)
  !! List of phonon energies on FBZ.
  real(dp), allocatable :: ph_ens_irred(:,:)
  !! List of phonon energies on IBZ.
  real(dp), allocatable :: ph_vels(:,:,:)
  !! List of phonon velocities on FBZ.
  real(dp), allocatable :: ph_vels_irred(:,:,:)
  !! List of phonon velocites on IBZ.
  complex(dp), allocatable :: ph_evecs(:,:,:)
  !! List of all phonon eigenvectors.
  complex(dp), allocatable :: ph_evecs_irred(:,:,:)
  !! List of IBZ wedge phonon eigenvectors.     

  !EPW Wannier representation data.
  integer(k4) :: numwannbands
  !! Number of Wannier bands.
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

  !I/O directories and system control flags.
  character(len=500) :: cwd
  !! Current working directory.
  character(len=500) ::datadumpdir
  !! Runtime data dump repository.
  character(len=500) :: g2dir
  !! Directory for e-ph vertex.
  character(len=500), allocatable :: Vdir
  !! Directory for ph-ph vertex.
  logical :: read_g2
  !! Choose if earlier e-ph vertex is to be used.
  logical :: read_V
  !! Choose if earlier p-ph vertex is to be used.
end module data



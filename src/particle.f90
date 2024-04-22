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

module particle_module
  !! Module containing the particle abstract data type.

  use precision, only: r64, i64
  
  implicit none

  private
  public particle

  type particle
     !! Data related to generic particle properties.

     integer(i64) :: numbands
     !! Total number of energy dispersion bands.
     integer(i64) :: wvmesh(3)
     !! Particle wave vector mesh.
     integer(i64) :: nwv
     !! Number of particle wave vectors in the full Brillouin zone (FBZ).
     integer(i64) :: nwv_irred
     !! Number of particle wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     real(r64), allocatable :: wavevecs(:,:)
     !! List of all particle wave vectors (crystal coordinates).
     real(r64), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible particle wave vectors (crystal coordinates).
     integer(i64), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(i64), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(i64), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(i64), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ particle wave vector point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(i64), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ particle wave vector point to its IBZ wedge image.
     integer(i64), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors.
     real(r64), allocatable :: symmetrizers(:,:,:)
     !! Symmetrizers of wave vector dependent vectors.
     integer(i64), allocatable :: simplicial_complex(:,:)
     !! List of all the wave vector mesh tetrahedra vertices.
     !! First axis list simplices(tetrahedra/triangles) and the second axis list the vertices.
     integer(i64), allocatable :: simplex_count(:)
     !! The number of simplicies in which a wave vector belongs.
     integer(i64), allocatable :: simplex_map(:,:,:)
     !! Mapping from a wave vector to the (simplex, vertex) where it belongs.
     real(r64), allocatable :: simplex_evals(:,:,:)
     !! Simplex vertices filled with eigenvalues.     
     real(r64), allocatable :: ens(:,:)
     !! List of particle energies on FBZ.
     real(r64), allocatable :: ens_irred(:,:)
     !! List of particle energies on IBZ.
     real(r64), allocatable :: vels(:,:,:)
     !! List of particle velocities on FBZ.
     real(r64), allocatable :: vels_irred(:,:,:)
     !! List of particle velocites on IBZ.
     complex(r64), allocatable :: evecs(:,:,:)
     !! List of all particle eigenvectors.
     complex(r64), allocatable :: evecs_irred(:,:,:)
     !! List of IBZ wedge particle eigenvectors.
     real(r64), allocatable :: dos(:,:)
     !! Band resolved density of states.
  end type particle
end module particle_module

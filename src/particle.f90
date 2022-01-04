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

module particle_module
  !! Module containing the particle abstract data type.

  use params, only: dp, k8
  
  implicit none

  private
  public particle

  type particle
     !! Data related to generic particle properties.

     integer(k8) :: numbands
     !! Total number of energy dispersion bands.
     integer(k8) :: wvmesh(3)
     !! Particle wave vector mesh.
     integer(k8) :: nwv
     !! Number of particle wave vectors in the full Brillouin zone (FBZ).
     integer(k8) :: nwv_irred
     !! Number of particle wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     real(dp), allocatable :: wavevecs(:,:)
     !! List of all particle wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible particle wave vectors (crystal coordinates).
     integer(k8), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(k8), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k8), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k8), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ particle wave vector point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k8), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ particle wave vector point to its IBZ wedge image.
     integer(k8), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors.
     real(dp), allocatable :: symmetrizers(:,:,:)
     !! Symmetrizers of wave vector dependent vectors.
     integer(k8), allocatable :: tetra(:,:)
     !! List of all the wave vector mesh tetrahedra vertices.
     !! First axis list tetraheda and the second axis list the vertices.
     integer(k8), allocatable :: tetracount(:)
     !! The number of tetrahedra in which a wave vector belongs.
     integer(k8), allocatable :: tetramap(:,:,:)
     !! Mapping from a wave vector to the (tetrahedron, vertex) where it belongs.
     real(dp), allocatable :: tetra_evals(:,:,:)
     !! Tetrahedra vertices filled with eigenvalues.
     integer(k8), allocatable :: triang(:,:)
     !! List of all the wave vector mesh triangles vertices.
     !! First axis lists triangles and the second axis lists the vertices.
     integer(k8), allocatable :: triangcount(:)
     !! The number of triangles in which a wave vector belongs.
     integer(k8), allocatable :: triangmap(:,:,:)
     !! Mapping from a wave vector to the (triangle, vertex) where it belongs.
     real(dp), allocatable :: triang_evals(:,:,:)
     !! Triangles vertices filled with eigenvalues.
     real(dp), allocatable :: ens(:,:)
     !! List of particle energies on FBZ.
     real(dp), allocatable :: ens_irred(:,:)
     !! List of particle energies on IBZ.
     real(dp), allocatable :: vels(:,:,:)
     !! List of particle velocities on FBZ.
     real(dp), allocatable :: vels_irred(:,:,:)
     !! List of particle velocites on IBZ.
     complex(dp), allocatable :: evecs(:,:,:)
     !! List of all particle eigenvectors.
     complex(dp), allocatable :: evecs_irred(:,:,:)
     !! List of IBZ wedge particle eigenvectors.
     real(dp), allocatable :: dos(:,:)
     !! Band resolved density of states.
  end type particle
end module particle_module

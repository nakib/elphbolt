! Copyright (C) 2022- Nakib Haider Protik <nakib.haider.protik@gmail.com>
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

module phonon_defect_module
  !! Module containing phonon defect related data type and procedures.

  use params, only: k8, dp
  use misc, only: exit_with_message, subtitle, demux_vector, twonorm
  use crystal_module, only: crystal
  use phonon_module, only: phonon

  implicit none

  private
  public phonon_defect

  type :: phonon_defect
     !! Data and procedures related to phonon defects.

     real(dp) :: range
     !! Radius of the defect in nm. This defines a block of cells in the defective supercell.
     integer(k8) :: numcells
     !! Number of cells in the defective supercell block.
     integer(k8), allocatable :: cell_pos_intvec(:, :)
     !! Unitcell positions as 0-based integer triplets in the defective supercell block.  
     real(dp), allocatable :: atom_pos(:, :, :)
     !! Cartesian positions of atoms in the defective supercell block.
     integer(k8), allocatable :: pcell_atom_label(:, :)
     !! Primitive cell equivalence (integer label) of atoms in the defective supercell block.
     real(dp), allocatable :: V_onsite_mass(:)
     !! On-site mass defect potential.
     real(dp), allocatable :: V_extended(:, :)
     !! General space-dependent, pairwise defect potential.
     logical mass_defect
     !! Choose if mass defect is going to be used.

   contains
     
     procedure, private :: generate_V_onsite_mass
  end type phonon_defect

contains

  subroutine initialize(self, ph, crys)
    !! Initialize the phonon defect data type.
    
    class(phonon_defect), intent(out) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    
    !Local variables
    real(dp) :: range, supercell_lattvecs(3, 3)
    logical :: mass_defect
    integer(k8) :: cell, atom, cell_count, cell_intvec(3), i, &
         supercell_numatoms, supercell_cell_pos_intvec(3, product(ph%scell))

    namelist /phonon_defect/ mass_defect, range
    
    call subtitle("Setting up phonon defect...")

    mass_defect = .false.
    range = 0.0_dp
    
    !Read defect input
    open(1, file = 'input.nml', status = 'old')

    read(1, nml = phonon_defect)

    if(mass_defect .and. range < 0) then
       call exit_with_message("Must provide non-zero range (in nm) of defect.")
    end if

    self%mass_defect = mass_defect
    self%range = range !nm
    
    !Apply defect radius
    cell_count = 0
    do cell = 1, product(ph%scell)
       !Demultiplex cell index into an 0-based integer triplet
       !giving the coordinates of a cell in the supercell.
       call demux_vector(cell_count, cell_intvec, ph%scell, 0_k8)
       
       !If position of cell is within range, then keep it.
       if(distance_from_origin(cell_intvec, ph%scell, crys%lattvecs) <= range) then
          cell_count = cell_count + 1

          supercell_cell_pos_intvec(:, cell_count) = cell_intvec 
       end if
    end do
    self%numcells = cell_count
    allocate(self%cell_pos_intvec(3, self%numcells))
    self%cell_pos_intvec(:, 1:self%numcells) = &
         supercell_cell_pos_intvec(:, 1:self%numcells)

    !Calculate supercell lattice vectors
    do i = 1, 3
       supercell_lattvecs(:, i) = ph%scell(i)*crys%lattvecs(:, i)
    end do

    !Number of atoms in the defective block of the supercell
    supercell_numatoms = self%numcells*crys%numatoms
    
    !Calculate defective supercell block atomic positions and
    !primitive cell equivalent atom labels.
    allocate(self%atom_pos(3, crys%numatoms, self%numcells))
    allocate(self%pcell_atom_label(crys%numatoms, self%numcells))
    do cell = 1, self%numcells
       do atom = 1, crys%numatoms
          self%atom_pos(:, atom, cell) = &
               matmul(supercell_lattvecs, &
               self%cell_pos_intvec(:, cell) + crys%basis(:, atom))

          self%pcell_atom_label(atom, cell) = crys%atomtypes(atom)
       end do
    end do
    
    !Create on-site mass perturbation
    if(self%mass_defect) call generate_V_onsite_mass(self, crys)
  end subroutine initialize
  
  subroutine generate_V_onsite_mass(self, crys)
    !! Create a single, on-site defect potential.
    !! The defect is always taken to be in the central unit cell.
    !
    ! V_M = -(M' - M)/M.
    ! Note that the expression evaluated here is defined purely on the real space.
    ! The phonon frequency squared scaling will be put back in later in the calculation.

    class(phonon_defect), intent(inout) :: self
    type(crystal), intent(in) :: crys

    !Local variables
    integer(k8) :: elem, cell

    allocate(self%V_onsite_mass(crys%numelements))
    self%V_onsite_mass = (1.0_dp - crys%subs_masses/crys%masses)
  end subroutine generate_V_onsite_mass

  pure real(dp) function distance_from_origin(cell_intvec, scell, lattvecs)
    !! Function to calculate the minimum distance (nm) of a vector measured from
    !! the origin of a supercell.

    integer(k8), intent(in) :: cell_intvec(3), scell(3)
    real(dp), intent(in) :: lattvecs(3, 3)

    !Local variables
    real(dp) :: distance_from_origins(3**3), supercell_lattvecs(3, 3)
    integer(k8) :: i, j, k, count
    
    !Calculate supercell lattice vectors
    do i = 1, 3
       supercell_lattvecs(:, i) = scell(i)*lattvecs(:, i)
    end do

    !Calculate the distance in the images of the central supercell
    count = 1
    do i = -1, 1
       do j = -1, 1
          do k = -1, 1
             distance_from_origins(count) = twonorm(&
                  matmul(supercell_lattvecs, cell_intvec - [i, j, k]*scell))
             count = count + 1
          end do
       end do
    end do

    !Return the minimum distance
    distance_from_origin = minval(distance_from_origins)
  end function distance_from_origin
end module phonon_defect_module

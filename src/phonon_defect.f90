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

  use params, only: k8, dp, hbar_eVps
  use misc, only: exit_with_message, subtitle, demux_vector, twonorm, write2file_rank2_real, &
       kronecker, expi, demux_state, invert, distribute_points, mux_state
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
     integer(k8) :: numhosts
     !! Number of host sites in the unit cell. This can't exceed the number of unique elements.
     integer(k8), allocatable :: cell_pos_intvec(:, :)
     !! Unitcell positions as 0-based integer triplets in the defective supercell block.  
     real(dp), allocatable :: atom_pos(:, :)
     !! Cartesian positions of atoms in the defective supercell block.
     integer(k8), allocatable :: pcell_atom_label(:)
     !! Primitive cell equivalence (integer label) of atoms in the defective supercell block.
     real(dp), allocatable :: V_mass(:)
     !! On-site mass defect potential.
     real(dp), allocatable :: V_bond(:, :)
     !! General space-dependent, pairwise defect potential.
     complex(dp), allocatable :: D0(:, :, :)
     !! Retarded, bare Green's function defined on the defect space.
     !complex(dp), allocatable :: irred_diagT(:, :, :)
     !!! Diagonal of T-matrix for irreducible phonons.
     logical mass_defect
     !! Choose if mass defect is going to be used.

   contains

     procedure, public :: initialize, calculate_phonon_Tmatrix
     procedure, private :: generate_V_mass, calculate_phonon_Tmatrix_host
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
         supercell_numatoms, supercell_cell_pos_intvec(3, product(ph%scell)), &
         atom_count

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
       call demux_vector(cell, cell_intvec, ph%scell, 0_k8)
       
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

!!$    !TODO Check if these are needed
!!$    !Calculate supercell lattice vectors
!!$    do i = 1, 3
!!$       supercell_lattvecs(:, i) = ph%scell(i)*crys%lattvecs(:, i)
!!$    end do

    !Number of atoms in the defective block of the supercell
    supercell_numatoms = self%numcells*crys%numatoms

    !TODO Check if atom_pos is even needed anywhere!
    !
    !Calculate defective supercell block atomic positions and
    !primitive cell equivalent atom labels.
    allocate(self%atom_pos(3, supercell_numatoms))
    allocate(self%pcell_atom_label(supercell_numatoms))
    atom_count = 0
    do cell = 1, self%numcells
       do atom = 1, crys%numatoms
          atom_count = atom_count + 1

          !DBG
          self%atom_pos(:, atom_count) = &
               matmul(crys%lattvecs, self%cell_pos_intvec(:, cell)) + &
               crys%basis_cart(:, atom)

          self%pcell_atom_label(atom_count) = atom
       end do
    end do
    
    print*, 'Number of cells in defective supercell:', self%numcells
    print*, self%atom_pos(:, 1)
    print*, self%atom_pos(:, 2)
    print*, self%atom_pos(:, 3)
    print*, self%atom_pos(:, 4)
    
    !Create on-site mass perturbation
    if(self%mass_defect) call generate_V_mass(self, crys)
  end subroutine initialize
  
  subroutine generate_V_mass(self, crys)
    !! Create a single, on-site defect potential.
    !! The defect is always taken to be in the central unit cell.
    !
    ! V_M = -(M' - M)/M.
    ! Note that the expression evaluated here is defined purely on the real space.
    ! The phonon frequency squared scaling will be put back in later in the calculation.

    class(phonon_defect), intent(inout) :: self
    type(crystal), intent(in) :: crys

    allocate(self%V_mass(crys%numelements))
    self%V_mass = (1.0_dp - crys%subs_masses/crys%masses)
  end subroutine generate_V_mass

  pure real(dp) function distance_from_origin(cell_intvec, scell, lattvecs)
    !! Function to calculate the minimum distance (nm) of a vector measured from
    !! the origin of a supercell.

    integer(k8), intent(in) :: cell_intvec(3), scell(3)
    real(dp), intent(in) :: lattvecs(3, 3)

    !Local variables
    real(dp) :: distance_from_origins(27), supercell_lattvecs(3, 3)
    integer(k8) :: i, j, k, count
    
    !Calculate supercell lattice vectors
    do i = 1, 3
       supercell_lattvecs(:, i) = scell(i)*lattvecs(:, i)
    end do

    !Calculate the distance in the images of the central supercell
    count = 0
    do i = -1, 1
       do j = -1, 1
          do k = -1, 1
             count = count + 1
             !TODO Check this calculation
             distance_from_origins(count) = twonorm(&
                  matmul(supercell_lattvecs, cell_intvec/dble(scell) - [i, j, k]))
          end do
       end do
    end do
    
    !Return the minimum distance
    distance_from_origin = minval(distance_from_origins)
  end function distance_from_origin

  subroutine calculate_phonon_Tmatrix(self, ph, crys, approx)

    class(phonon_defect), intent(inout) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    character(len=*), intent(in) :: approx

    integer(k8) :: host, numhosts, ik
    real(dp) :: def_frac
    real(dp), allocatable :: scatt_rates(:, :)
    complex(dp), allocatable :: irred_diagT(:, :, :)
    
    !DBG
    numhosts = 2
    
    allocate(irred_diagT(ph%nwv_irred, ph%numbands, numhosts))

    allocate(scatt_rates(ph%nwv_irred, ph%numbands))
    scatt_rates = 0.0_dp
    
    do host = 1, numhosts
       call self%calculate_phonon_Tmatrix_host(ph, crys, approx, crys%defect_hosts(host), &
            irred_diagT(:, :, host))

       def_frac = crys%subs_conc(crys%atomtypes(crys%defect_hosts(host)))*(1.0e-21_dp*crys%volume)
       
       do ik = 1, ph%nwv_irred
          scatt_rates(ik, :) = scatt_rates(ik, :) + &
               def_frac*imag(irred_diagT(ik, :, host))/ph%ens(ph%indexlist_irred(ik), :)
       end do

       if(this_image() == 1) then
          print*, 'Calculating ph-defect scattering rates at host site', crys%defect_hosts(host)
          print*, 'Host atom type', crys%atomtypes(crys%defect_hosts(host))
          print*, 'Defect conc', crys%subs_conc(crys%atomtypes(crys%defect_hosts(host)))
       end if
    end do
    scatt_rates = -scatt_rates/hbar_eVps

    !Deal with Gamma point acoustic phonons!
    scatt_rates(1, 1:3) = 0.0_dp
    
    !Write to file
    call write2file_rank2_real(ph%prefix // '.W_rta_'//ph%prefix//'defect', scatt_rates)
  end subroutine calculate_phonon_Tmatrix
  
  subroutine calculate_phonon_Tmatrix_host(self, ph, crys, approx, host, diagT)
    !! Parallel calculator of the scattering T-matrix for phonons for a given approximation.
    !!
    !! D0 Retarded, bare Green's function in real space
    !! V_mass On-site scattering potential in real space
    !! diagT Diagonoal of the scattering T-matrix in reciprocal space
    !! approx Approximation for the scattering theory

    class(phonon_defect), intent(in) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    character(len=*), intent(in) :: approx
    complex(dp), intent(out) :: diagT(:, :)
    integer(k8), intent(in) :: host

    !Local variables
    integer(k8) :: num_dof_def, numstates_irred, istate, &
         chunk, start, end, num_active_images, i, j, a, def_numatoms, &
         dof_counter, iq, s, tau_sc, tau_uc, cell, atom
    complex(dp), allocatable :: inv_one_minus_VD0(:, :), T(:, :, :), phi(:)
    real(dp), allocatable :: V(:, :), identity(:, :)
    complex(dp) :: phase, ev(ph%numbands, ph%numbands)
    real(dp) :: q_cart(3), en_sq, val

    !Displacement degrees of freedom in the defective supercell 
    num_dof_def = size(self%D0, 1)

    !Number of atoms in the defective block of the supercell
    def_numatoms = num_dof_def/3

    !Number of irreducible phonon states
    numstates_irred = size(self%D0, 3)

    allocate(T(num_dof_def, num_dof_def, numstates_irred))
    allocate(V(num_dof_def, num_dof_def))
    T = 0.0_dp
    V = 0.0_dp

    allocate(identity(num_dof_def, num_dof_def))
    identity = 0.0_dp
    do i = 1, num_dof_def
       identity(i, i) = 1.0_dp
    end do

    if(approx == 'full Born') allocate(inv_one_minus_VD0(num_dof_def, num_dof_def))

    !Divide phonon states among images
    call distribute_points(numstates_irred, chunk, start, end, num_active_images)

    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)

       !Squared energy of phonon 1
       en_sq = ph%ens(ph%indexlist_irred(iq), s)**2

       !Scale the mass perturbation with squared energy.
       !Since on-site mass perturbation is forced to be in the unit cell,
       !only these elements of V get a contribution.
       dof_counter = 0
       do a = 1, crys%numatoms !Number of atoms in the unit cell
          val = en_sq*self%V_mass(crys%atomtypes(a))*kronecker(a, host)
          do i = 1, 3 !Cartesian directions
             dof_counter = dof_counter + 1
             V(dof_counter, dof_counter) = val
          end do
       end do

       select case(approx)
       case('lowest order')
          ! Lowest order:
          ! T = V
          !                            
          !    *                     
          !    |                        
          !  V |                         
          !    |
          !                
          T(:, :, istate) = V
       case('1st Born')
          ! 1st Born approximation:
          ! T = V + V.D0.V
          !                            
          !    *             *            
          !    |            / \              
          !  V |     +     /   \              
          !    |          /_____\
          !                 D0
          !
          !T(:, :, istate) = V + matmul(matmul(V, self%D0(:, :, istate)), V)
          T(:, :, istate) = matmul(V, matmul(self%D0(:, :, istate), V))
       case('full Born')
          ! Full Born approximation:
          ! T = V + V.D0.T = [I - VD0]^-1 . V
          !                              _                                      _
          !    *             *          |    *         *            *            |
          !    |            /           |    |        / \          /|\           |
          !  V |     +     /       x    |    |   +   /   \   +    / | \  +  ...  |
          !    |          /_____        |_   |      /_____\      /__|__\        _|
          !                 D0           
          !
          inv_one_minus_VD0 = identity - matmul(V, self%D0(:, :, istate))
          call invert(inv_one_minus_VD0)
          T(:, :, istate) = matmul(inv_one_minus_VD0, V)
       case default
          call exit_with_message("T-matrix approximation not recognized.")
       end select
    end do

    !Reduce T
    sync all
    call co_sum(T)
    sync all
    
    !Release some memory
    deallocate(V, identity)
    if(allocated(inv_one_minus_VD0)) deallocate(inv_one_minus_VD0)
    
    !Calculate diagonal T in reciprocal space.
    allocate(phi(num_dof_def))
    diagT = 0.0_dp
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)

       !This phonon eigenvector
       ev = ph%evecs(ph%indexlist_irred(iq), :, :)
       !ev = ph%evecs(iq, :, :)
       !ev = conjg(ph%evecs(ph%indexlist_irred(iq), :, :))
       !ev = transpose(conjg(ph%evecs(ph%indexlist_irred(iq), :, :)))
       !ev = transpose(ph%evecs(ph%indexlist_irred(iq), :, :))
       !ev = ph%evecs(ph%indexlist_irred(iq), :, :)/ph%nequiv(iq)
       !ev = ph%evecs(ph%indexlist_irred(iq), :, :)*ph%nequiv(iq)
       
       !Calculate wave vector in Cartesian coordinates
       q_cart = matmul(crys%reclattvecs, ph%wavevecs_irred(iq, :))

       !TODO this should be done in a separate subroutine.
       !Precompute the eigenfunction at (q, s)
       dof_counter = 0
       do cell = 1, self%numcells
          !Phase factor
          !DBG
          !TODO Check units below. def_supercell_cell_pos is in integer triplet form.
          phase = 1.0_dp !expi( &
               !dot_product(q_cart, def_supercell_cell_pos(:, cell)) )
          
          do atom = 1, crys%numatoms
             tau_sc = mux_state(crys%numatoms, atom, cell) 

             !phase = expi( &
             !     dot_product(q_cart, self%atom_pos(:, tau_sc)) )
             
             !Get primitive cell equivalent atom of supercell atom
             tau_uc = self%pcell_atom_label(tau_sc)
             
             !Run over Cartesian directions
             do a = 1, 3
                dof_counter = dof_counter + 1
                phi(dof_counter) = phase*ev(s, mux_state(3_k8, a, tau_uc))
             end do
          end do
       end do

       !Fourier transform
       diagT(iq, s) = dot_product(phi, matmul(T(:, :, istate), phi))
    end do

    !Reduce T
    sync all
    call co_sum(diagT)
    sync all
  end subroutine calculate_phonon_Tmatrix_host
end module phonon_defect_module

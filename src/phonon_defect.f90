! Copyright 2022 elphbolt contributors.
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

  use params, only: i64, r64, hbar_eVps, twopi, pi
  use misc, only: exit_with_message, subtitle, demux_vector, twonorm, write2file_rank2_real, &
       kronecker, expi, demux_state, invert, distribute_points, mux_state
  use crystal_module, only: crystal
  use phonon_module, only: phonon

  implicit none

  private
  public phonon_defect

  type :: phonon_defect
     !! Data and procedures related to phonon defects.

     real(r64) :: range
     !! Radius of the defect in nm. This defines a block of cells in the defective supercell.
     integer(i64) :: numcells
     !! Number of cells in the defective supercell block.
     integer(i64) :: numhosts
     !! Number of host sites in the unit cell. This can't exceed the number of unique elements.
     integer(i64), allocatable :: cell_pos_intvec(:, :), dimp_cell_pos_intvec(:, :)
     !! Unitcell positions as 0-based integer triplets in the defective supercell block.  
     integer(i64), allocatable :: pcell_atom_label(:)
     !! Primitive cell equivalence (integer label) of atoms in the defective supercell block.
     integer(i64), allocatable :: pcell_atom_dof(:)
     !! Primitive cell equivalent atomic degree of freedom.
     real(r64), allocatable :: V_mass(:)
     !! On-site mass defect potential.
     real(r64), allocatable :: V_bond(:, :)
     !! General space-dependent, pairwise defect potential.
     complex(r64), allocatable :: D0(:, :, :)
     !! Retarded, bare Green's function defined on the defect space.
     logical mass_defect
     !! Choose if mass defect is going to be used.
     character(len=100) :: approx
     !! Approximation of scattering T-matrix.
     
   contains

     procedure, public :: initialize, calculate_phonon_Tmatrix
     procedure, private :: calculate_phonon_Tmatrix_host
  end type phonon_defect

contains

  subroutine initialize(self, ph, crys)
    !! Initialize the phonon defect data type.
    
    class(phonon_defect), intent(out) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    
    !Local variables
    real(r64) :: range
    logical :: mass_defect
    integer(i64) :: cell, atom, cell_count, cell_intvec(3), &
         supercell_numatoms, atom_count, supercell_cell_pos_intvec(3, product(ph%scell)), &
         uc_dof_count, sc_dof_count, a
    character(len=100) :: approx

    namelist /phonon_defect/ mass_defect, range, approx
    
    call subtitle("Setting up phonon defect...")

    mass_defect = .false.
    range = 0.0_r64
    approx = 'full Born'
    
    !Read defect input
    open(1, file = 'input.nml', status = 'old')

    read(1, nml = phonon_defect)

    if(mass_defect .and. range < 0) then
       call exit_with_message("Must provide non-zero range (in nm) of defect.")
    end if

    if(.not. (approx .eq. 'lowest order' .or. &
         approx .eq. '1st Born' .or. &
         approx .eq. 'full Born')  ) then
       call exit_with_message(&
            "T-matrix approximation must be either 'lowest order', '1st Born', or 'full Born.'")
    end if

    self%mass_defect = mass_defect
    self%range = range !nm
    self%approx = approx
    
    !Apply defect radius
    cell_count = 0
    do cell = 1, product(ph%scell)
       !Demultiplex cell index into an 0-based integer triplet
       !giving the coordinates of a cell in the supercell.
       call demux_vector(cell, cell_intvec, ph%scell, 0_i64)
       !cell_intvec = cell_intvec - ph%scell/2
       
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

    !Number of atoms in the defective block of the supercell
    supercell_numatoms = self%numcells*crys%numatoms

    !Calculate defective supercell block atomic positions and
    !primitive cell equivalent atom labels.
    allocate(self%pcell_atom_label(supercell_numatoms))
    allocate(self%pcell_atom_dof(supercell_numatoms*3))

    allocate(self%dimp_cell_pos_intvec(3, supercell_numatoms*3))
    
    atom_count = 0
    sc_dof_count = 0
    do cell = 1, self%numcells
       uc_dof_count = 0
       do atom = 1, crys%numatoms
          atom_count = atom_count + 1

          self%pcell_atom_label(atom_count) = atom
          
          do a = 1, 3
             uc_dof_count = uc_dof_count + 1
             sc_dof_count = sc_dof_count + 1
             
             self%pcell_atom_dof(sc_dof_count) = uc_dof_count

             self%dimp_cell_pos_intvec(:, sc_dof_count) = self%cell_pos_intvec(:, cell)
          end do
       end do
    end do

    if(this_image() == 1) then
       write(*, "(A, A)") 'T-matrix approximation level: ', self%approx
       write(*,"(A, F7.2, A)") 'Defect range = ', self%range, ' nm'
       if(self%mass_defect) write(*,"(A)") 'On-site mass-defect scattering will be used.'
    end if
  end subroutine initialize
  
  pure real(r64) function distance_from_origin(cell_intvec, scell, lattvecs)
    !! Function to calculate the minimum distance (nm) of a vector measured from
    !! the origin of a supercell.

    integer(i64), intent(in) :: cell_intvec(3), scell(3)
    real(r64), intent(in) :: lattvecs(3, 3)

    !Local variables
    real(r64) :: distance_from_origins(5**3), supercell_lattvecs(3, 3)
    integer(i64) :: i, j, k, count
    
    !Calculate supercell lattice vectors
    do i = 1, 3
       supercell_lattvecs(:, i) = scell(i)*lattvecs(:, i)
    end do

    !Calculate the distance in the images of the central supercell
    count = 0
    do i = -2, 2
       do j = -2, 2
          do k = -2, 2
             count = count + 1
             distance_from_origins(count) = twonorm(&
                  matmul(supercell_lattvecs, cell_intvec/dble(scell) - [i, j, k]))
          end do
       end do
    end do
    
    !Return the minimum distance
    distance_from_origin = minval(distance_from_origins)
  end function distance_from_origin
  
  subroutine calculate_phonon_Tmatrix(self, ph, crys)

    class(phonon_defect), intent(inout) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
  
    integer(i64) :: host, ik, i, j, dopant
    real(r64) :: def_frac
    real(r64), allocatable :: scatt_rates(:, :), renorm_ens(:, :), lineshifts(:, :)
    complex(r64), allocatable :: irred_diagT(:, :, :)

    real(r64) :: V_mass_iso(crys%numelements)
    integer(i64) :: num_atomtypes(crys%numelements)
    
    num_atomtypes(:) = 0_i64
    do i = 1, crys%numelements
       do j = 1, crys%numatoms
          if(crys%atomtypes(j) == i) num_atomtypes(i) = num_atomtypes(i) + 1 
       end do
    end do

    allocate(irred_diagT(ph%nwv_irred, ph%numbands, crys%numelements))

    allocate(scatt_rates(ph%nwv_irred, ph%numbands), &
         lineshifts(ph%nwv_irred, ph%numbands), renorm_ens(ph%nwv_irred, ph%numbands))
    
    scatt_rates = 0.0_r64
    lineshifts = 0.0_r64
    
    if(self%mass_defect) then
       do host = 1, crys%numelements
          do dopant = 1, crys%numdopants_types(host) !dopants of this host atom
             V_mass_iso = 0.0_r64
             V_mass_iso(host) = 1.0_r64 - crys%dopant_masses(dopant, host)/crys%masses(host)

             call self%calculate_phonon_Tmatrix_host(ph, crys, host, &
                  irred_diagT(:, :, host), V_mass_iso)

             def_frac = crys%dopant_conc(dopant, host)*(1.0e-21_r64*crys%volume)/num_atomtypes(host)

             do ik = 1, ph%nwv_irred
                scatt_rates(ik, :) = scatt_rates(ik, :) + &
                     def_frac*imag(irred_diagT(ik, :, host))/ph%ens(ph%indexlist_irred(ik), :)

                lineshifts(ik, :) = lineshifts(ik, :) + &
                     def_frac*real(irred_diagT(ik, :, host))
             end do
          end do
       end do
    end if

    scatt_rates = -scatt_rates/hbar_eVps

    do ik = 1, ph%nwv_irred
       renorm_ens(ik, :) = sqrt(ph%ens(ph%indexlist_irred(ik), :)**2 + &
            lineshifts(ik, :))
    end do
    
    !Deal with Gamma point acoustic phonons
    scatt_rates(1, 1:3) = 0.0_r64

    !Write to file
    call write2file_rank2_real(ph%prefix // '.W_rta_'//ph%prefix//'defect', scatt_rates)
    !call write2file_rank2_real(ph%prefix // '.lineshifts_ibz_'//ph%prefix//'defect', lineshifts)
    call write2file_rank2_real(ph%prefix // '.ens_renorm_ibz_'//ph%prefix//'defect', renorm_ens)
  end subroutine calculate_phonon_Tmatrix
  
  subroutine calculate_phonon_Tmatrix_host(self, ph, crys, host_atom_type, diagT, V_mass)
    !! Parallel calculator of the scattering T-matrix for phonons for a given approximation.
    !!
    !! D0 Retarded, bare Green's function in real space
    !! V_mass On-site scattering potential in real space
    !! diagT Diagonoal of the scattering T-matrix in reciprocal space
    !! approx Approximation for the scattering theory

    class(phonon_defect), intent(in) :: self
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    complex(r64), intent(out) :: diagT(:, :)
    real(r64), intent(in) :: V_mass(crys%numelements)
    integer(i64), intent(in) :: host_atom_type

    !Local variables
    integer(i64) :: num_dof_def, numstates_irred, istate, &
         chunk, start, end, num_active_images, i, a, def_numatoms, &
         dof_counter, iq, s, cell, atom
    complex(r64), allocatable :: inv_one_minus_VD0(:, :), T(:, :, :), phi(:)
    real(r64), allocatable :: V(:, :), identity(:, :)
    complex(r64) :: phase, ev(ph%numbands)
    real(r64) :: en_sq, val

    !Displacement degrees of freedom in the defective supercell 
    num_dof_def = size(self%D0, 1)

    !Number of atoms in the defective block of the supercell
    def_numatoms = num_dof_def/3

    !Number of irreducible phonon states
    numstates_irred = ph%nwv_irred*ph%numbands
    
    allocate(T(num_dof_def, num_dof_def, numstates_irred))
    allocate(V(num_dof_def, num_dof_def))
    T = 0.0_r64
    V = 0.0_r64

    allocate(identity(num_dof_def, num_dof_def))
    identity = 0.0_r64
    do i = 1, num_dof_def
       identity(i, i) = 1.0_r64
    end do

    if(self%approx == 'full Born') allocate(inv_one_minus_VD0(num_dof_def, num_dof_def))

    !Divide phonon states among images
    call distribute_points(numstates_irred, chunk, start, end, num_active_images)
    
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)

       !Squared energy of phonon 1
       en_sq = ph%ens(ph%indexlist_irred(iq), s)**2
       
       !Scale the mass perturbation with squared energy.
       dof_counter = 0
       do cell = 1, self%numcells
          do atom = 1, crys%numatoms !Number of atoms in the unit cell
             val = en_sq*V_mass(crys%atomtypes(atom))
             do a = 1, 3 !Cartesian directions
                dof_counter = dof_counter + 1
                !Since on-site mass perturbation is forced to be in the central unit cell,
                !only these elements of V get a non-zero contribution.
                if(all(self%cell_pos_intvec(:, cell) == 0)) then
                   V(dof_counter, dof_counter) = val
                end if
             end do
          end do
       end do
       
       select case(self%approx)
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
          T(:, :, istate) = V + matmul(V, matmul(self%D0(:, :, istate), V))
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
    diagT = 0.0_r64
    
    !Optical theorem
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)
       
       !This phonon eigenvector
       ev = ph%evecs(ph%indexlist_irred(iq), s, :)
       
       do dof_counter = 1, num_dof_def
          phase = expi( &
               twopi*dot_product(ph%wavevecs_irred(iq, :), &
               self%dimp_cell_pos_intvec(:, dof_counter)) )

          phi(dof_counter) = phase*ev(self%pcell_atom_dof(dof_counter))
       end do
       
       !<i|T|j> --> <sq|T|sq>
       diagT(iq, s) = dot_product(phi, matmul(T(:, :, istate), phi))
    end do
    
    !Reduce T
    sync all
    call co_sum(diagT)
    sync all
  end subroutine calculate_phonon_Tmatrix_host
end module phonon_defect_module

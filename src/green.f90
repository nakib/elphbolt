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

module Green_function
  !! Module containing Green's function related procedures.

  use params, only: i64, r64, pi, oneI, twopi, hbar_eVps
  use electron_module, only: electron
  use phonon_module, only: phonon
  use crystal_module, only: crystal
  use delta, only: real_tetra, delta_fn, get_delta_fn_pointer
  use misc, only: exit_with_message, distribute_points, expi, demux_state, invert, &
       write2file_rank2_real, kronecker, mux_state
  
  implicit none

  private
  public calculate_retarded_phonon_D0, resolvent
  
contains
  
  complex(r64) function resolvent(species, ib, iwv, sampling_point)
    !! Calculate the resolvant
    !!   electron: 1/[z - E(k)], lim z -> E + i0^{+}.
    !!   phonon:   1/[z - omega^2(q)], lim z -> omega^2 + i0^{+}.
    !!
    !! species Object of particle type
    !! ib, iwv Band, wave vector indices
    !! sampling_point Sampling energy (or energy squared) in eV (or eV^2), depending on particle type.

    class(*), intent(in) :: species
    integer(i64), intent(in) :: ib, iwv
    real(r64), intent(in) :: sampling_point

    !Local variables
    real(r64) :: Im_resolvent, Re_resolvent
    procedure(delta_fn), pointer :: delta_fn_ptr => null()

    !TODO Pass delta_fn_ptr to this function.
    
    !Associate delta function procedure pointer
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra = .true.)
    
    select type(species)
    class is(phonon)
       !Imaginary part of resolvent
       if(sampling_point < 1.0e-3) then !=>omega < 1e-3 eV^2
          !For very small energies, use the tetrahedra filled with omega.
          !This works better, numerically.
          Im_resolvent = -pi*delta_fn_ptr(sqrt(sampling_point), iwv, ib, species%wvmesh, species%simplex_map, &
               species%simplex_count, species%simplex_evals)/2.0_r64/sqrt(sampling_point)
       else
          !Otherwise, use the omega^2 tetrahedra
          Im_resolvent = -pi*delta_fn_ptr(sampling_point, iwv, ib, species%wvmesh, &
               species%simplex_map, &
               species%simplex_count, species%simplex_squared_evals)
       end if

       !Real part of resolvent
       Re_resolvent = real_tetra(sampling_point, iwv, ib, species%wvmesh, &
            species%simplex_map, &
            species%simplex_count, species%simplex_squared_evals)
    class is(electron)
       !Imaginary part of resolvent
       Im_resolvent = -pi*delta_fn_ptr(sampling_point, iwv, ib, species%wvmesh, &
            species%simplex_map, &
            species%simplex_count, species%simplex_evals)

       !Real part of resolvent
       Re_resolvent = real_tetra(sampling_point, iwv, ib, species%wvmesh, &
            species%simplex_map, &
            species%simplex_count, species%simplex_evals)
    class default
       call exit_with_message(&
            "Unknown particle species in resolvent. Exiting.")
    end select

    resolvent = Re_resolvent + oneI*Im_resolvent

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end function resolvent

  subroutine calculate_retarded_phonon_D0(ph, crys, def_supercell_cell_pos_intvec, &
       pcell_atom_label, D0, dimp_cell_pos_intvec, pcell_atom_dof)
    !! Parallel driver of the retarded, bare phonon Green's function, D0, over
    !! the IBZ states.
    !!
    !! ph Phonon object
    !! crys Crystal object
    !! def_supercell_cell_pos_intvec Positions of unitcells (integer 3 vector) in the defective supercell
    !! pcell_atom_label Primitive cell equivalence of atom labels in the defective supercell
    !! D0 Green's function

    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    integer(i64), intent(in) :: def_supercell_cell_pos_intvec(:, :)
    integer(i64), intent(in) :: pcell_atom_label(:)

    integer(i64), intent(in) :: dimp_cell_pos_intvec(:, :), pcell_atom_dof(:)
    
    complex(r64), allocatable, intent(out) :: D0(:, :, :)
    
    !Local variables
    integer(i64) :: nstates_irred, chunk, start, end, num_active_images, &
         istate1, s1, iq1_ibz, iq1, s2, iq2, j, num_dof_def, dof_counter, iq, &
         def_numatoms, def_numcells
    real(r64) :: en1_sq
    complex(r64) :: d0_istate, phase, ev(ph%numbands, ph%numbands)!, dos(ph%nwv_irred, ph%numbands)
    complex(r64), allocatable :: phi(:, :, :), phi_internal(:)

    !Total number of atoms in the defective block of the supercell
    def_numatoms = size(pcell_atom_label)
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    !Total number of degrees of freedom in defective supercell
    num_dof_def = def_numatoms*3
    
    !Number of primitive unit cells in the defective supercell    
    def_numcells = num_dof_def/crys%numatoms/3
    
    allocate(phi(num_dof_def, ph%numbands, ph%nwv), phi_internal(num_dof_def), &
         D0(num_dof_def, num_dof_def, nstates_irred))

    !Precompute the defective supercell eigenfunctions
    do iq = 1, ph%nwv
       !This phonon eigenvector
       ev = ph%evecs(iq, :, :)

       do dof_counter = 1, num_dof_def
          phase = expi( &
               twopi*dot_product(ph%wavevecs(iq, :), &
               dimp_cell_pos_intvec(:, dof_counter)) )
          
          phi(dof_counter, :, iq) = &
               phase*ev(:, pcell_atom_dof(dof_counter))
       end do
    end do

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    !Initialize D0
    D0 = 0.0_r64
    
    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over first (IBZ) phonon states
       do istate1 = start, end
          !Demux state index into branch (s) and wave vector (iq) indices
          call demux_state(istate1, ph%numbands, s1, iq1_ibz)

          !Muxed index of wave vector from the IBZ index list.
          !This will be used to access IBZ information from the FBZ quantities.
          !Recall that for the phonons, unlike for the electrons, we don't save these IBZ info separately.
          iq1 = ph%indexlist_irred(iq1_ibz)

          !Squared energy of phonon 1
          en1_sq = ph%ens(iq1, s1)**2

          !Sum over internal (FBZ) phonon wave vectors
          do iq2 = 1, ph%nwv
             !Sum over internal phonon bands
             do s2 = 1, ph%numbands
                phi_internal(:) = phi(:, s2, iq2)

                d0_istate = resolvent(ph, s2, iq2, en1_sq)

                do j = 1, num_dof_def
                   D0(:, j, istate1) = D0(:, j, istate1) + &
                        d0_istate*phi_internal(:)*conjg(phi_internal(j))
                end do
             end do
          end do
       end do
    end if

    !Reduce D0
    sync all
    call co_sum(D0)
    sync all

!!$    !Sanity check: print DOS
!!$    dos = 0.0_r64
!!$    !Only work with the active images
!!$    if(this_image() <= num_active_images) then
!!$       do istate1 = start, end
!!$          !Demux state index into branch (s) and wave vector (iq) indices
!!$          call demux_state(istate1, ph%numbands, s1, iq1_ibz)
!!$
!!$          iq1 = ph%indexlist_irred(iq1_ibz)
!!$
!!$          do i = 1, num_dof_def
!!$             dos(iq1_ibz, s1) = dos(iq1_ibz, s1) + &
!!$                  D0(i, i, istate1)
!!$          end do
!!$
!!$          dos(iq1_ibz, s1) = dos(iq1_ibz, s1)*ph%ens(iq1, s1)
!!$       end do
!!$    end if
!!$
!!$    !Reduce dos
!!$    sync all
!!$    call co_sum(dos)
!!$    sync all
!!$
!!$    call write2file_rank2_real(ph%prefix // '.D0test_'//ph%prefix//'dos', imag(-2.0/pi*dos))
!!$    sync all
!!$    !!
  end subroutine calculate_retarded_phonon_D0
end module Green_function

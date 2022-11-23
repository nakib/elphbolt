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

module Green_function
  !! Module containing Green's function related procedures.

  use params, only: k8, dp, pi, oneI
  use electron_module, only: electron
  use phonon_module, only: phonon
  use crystal_module, only: crystal
  use delta, only: delta_fn_tetra, real_tetra
  use misc, only: exit_with_message, distribute_points, expi, demux_state, invert, &
       write2file_rank2_real, kronecker, mux_state
  
  implicit none

  private
  public calculate_retarded_phonon_D0
  
contains
  
  complex(dp) function resolvent(species, ib, iwv, sampling_point)
    !! Calculate the resolvant
    !!   electron: 1/[z - E(k)], lim z -> E + i0^{+}.
    !!   phonon:   1/[z - omega^2(q)], lim z -> omega^2 + i0^{+}.
    !!
    !! species Object of particle type
    !! ib, iwv Band, wave vector indices
    !! sampling_point Sampling energy (or energy squared) in eV (or eV^2), depending on particle type.

    class(*), intent(in) :: species
    integer(k8), intent(in) :: ib, iwv
    real(dp), intent(in) :: sampling_point

    !Local variables
    real(dp) :: Im_resolvent, Re_resolvent

    select type(species)
    class is(phonon)
       !Imaginary part of resolvent
       Im_resolvent = -pi*delta_fn_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_squared_evals)

       !Real part of resolvent
       Re_resolvent = real_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_squared_evals)
    class is(electron)
       !Imaginary part of resolvent
       Im_resolvent = -pi*delta_fn_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_evals)

       !Real part of resolvent
       Re_resolvent = real_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_evals)
    class default
       call exit_with_message(&
            "Unknown particle species in resolvent. Exiting.")
    end select

    resolvent = Re_resolvent + oneI*Im_resolvent
  end function resolvent

!!$  subroutine calculate_retarded_phonon_D0(ph, crys, def_supercell_cell_pos_intvec, &
!!$       pcell_atom_label, D0)
  subroutine calculate_retarded_phonon_D0(ph, crys, def_atom_pos, &
       pcell_atom_label, D0)
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
    !integer(k8), intent(in) :: def_supercell_cell_pos_intvec(:, :) !Cartesian
    real(dp), intent(in) :: def_atom_pos(:, :)
    integer(k8), intent(in) :: pcell_atom_label(:)
    complex(dp), allocatable, intent(out) :: D0(:, :, :)
    
    !Local variables
    integer(k8) :: nstates_irred, chunk, start, end, num_active_images, &
         istate1, s1, iq1_ibz, iq1, s2, iq2, i, j, num_dof_def, a, dof_counter, iq, &
         tau_sc, tau_uc, def_numatoms, def_numcells, atom, cell
    real(dp) :: en1_sq, q_cart(3)
    complex(dp) :: d0_istate, phase, ev(ph%numbands, ph%numbands)
    complex(dp), allocatable :: phi(:, :, :), phi_internal(:)

    !Total number of atoms in the defective block of the supercell
    def_numatoms = size(pcell_atom_label)
    
    !Total number of IBZ and FBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    !Total number of degrees of freedom in defective supercell
    num_dof_def = def_numatoms*3

    !Number of primitive unit cells in the defective supercell    
    !def_numcells = size(def_supercell_cell_pos_intvec, 2)
    def_numcells = num_dof_def/def_numatoms/3
    
    allocate(phi(num_dof_def, ph%numbands, ph%nwv), phi_internal(num_dof_def), &
         D0(num_dof_def, num_dof_def, nstates_irred))

    !Precompute the defective supercell eigenfunctions
    do iq = 1, ph%nwv
       !Calculate wave vector in Cartesian coordinates
       q_cart = matmul(crys%reclattvecs, ph%wavevecs(iq, :))

       !This phonon eigenvector
       ev = ph%evecs(iq, :, :)
       
       !TODO this should be done in a separate subroutine.
       dof_counter = 0
       do cell = 1, def_numcells
          !TODO Check units below. def_supercell_cell_pos_intvec is in integer triplet form.
          !DBG
          phase = 1.0_dp !expi( &
          !dot_product(q_cart, def_supercell_cell_pos_intvec(:, cell)) )

          do atom = 1, crys%numatoms
             !Index of basis atom in supercell
             tau_sc = mux_state(crys%numatoms, atom, cell)

             !phase = expi( &
             !     dot_product(q_cart, def_atom_pos(:, tau_sc)) )
             
             !Get primitive cell equivalent atom of supercell atom
             tau_uc = pcell_atom_label(tau_sc)

             !print*, tau_uc, (tau_uc - 1)*3 + 1, (tau_uc - 1)*3 + 3
             
             !Run over Cartesian directions
             do a = 1, 3
                dof_counter = dof_counter + 1
                phi(dof_counter, :, iq) = &
                     phase*ev(:, mux_state(3_k8, a, tau_uc))
             end do
          end do
       end do
    end do

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    !Initialize D0
    D0 = 0.0_dp
    
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
                     d0_istate*conjg(phi_internal(j))*phi_internal(:)
             end do
          end do
       end do
    end do

    !Reduce D0
    sync all
    call co_sum(D0)
    sync all
  end subroutine calculate_retarded_phonon_D0
end module Green_function

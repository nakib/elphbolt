module Green_function
  !! Module containing Green's function related procedures.

  use params, only: k8, dp, twopi
  use electron_module, only: electron
  use phonon_module, only: phonon
  use misc, only: exit_with_message, distribute_points, expi
  
  implicit none

  private
  public calculate_retarded_G0
  
contains
  
  pure complex(dp) function resolvent(species, ib, iwv, sampling_point)
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
    real(dp) :: ImG0, ReG0

    select type(species)
    class is(phonon)
       !Imaginary part of G0
       ImG0 = -pi*delta_fn_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_squared_evals)

       !Real part of G0
       ReG0 = real_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_squared_evals)
    class is(electron)
       !Imaginary part of G0
       ImG0 = -pi*delta_fn_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_evals)

       !Real part of G0
       ReG0 = real_tetra(sampling_point, iwv, ib, species%wvmesh, species%tetramap, &
            species%tetracount, species%tetra_evals)
    class default
       call exit_with_message(&
            "Unknown particle species in resolvent. Exiting.")
    end select

    resolvent = complex(ReG0, ImG0, type = dp)
  end function resolvent

  subroutine calculate_retarded_phonon_D0(ph, def)
    !! Parallel driver of the retarded, bare phonon Green's function, D0, over
    !! the IBZ states.
    !!
    !! TODO Description of the subroutine
    !!

    type(phonon), intent(in) :: ph
    type(defect), intent(inout) :: def !TODO have to define this

    !Local variables
    integer(k8) :: nstates, nstates_irred, chunk, start, end, num_active_images, &
         istate1, s1, iq1_ibz, iq1, s2, iq2, i, j, pcell_i, pcell_j, &
         i_disp, j_disp, num_dof_def
    real(dp) :: en1_sq, q2(3)
    complex(dp), allocatable :: phi(:, :)
   
    !Total number of IBZ and FBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands
    nstates = ph%nwv*ph%numbands

    !Total number of degrees of freedom in defective supercell
    num_dof_def = def%num_def_cells*3
    
    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    allocate(phi(num_dof_def, num_dof_def))

    !Precompute the defective supercell eigenvectors
    do iq = 1, ph%nwv
       dof_counter = 0
       do tau_sc = 1, def%num_def_supercell !Atoms in defective supercell
          !Phase factor
          phase = expi( &
               twopi*dot_product(ph%wavevecs(iq, :), def%r(tau_sc)) )

          !Get primitive cell equivalent atom of supercell atom
          tau_uc = def%pcell_equiv_atom_label(tau_sc)

          !Run over Cartesian directions
          do a = 1, 3
             dof_counter = dof_counter + 1
             phi(iq, :, dof_counter) = phase*ph%evecs(iq, :, (tau_uc - 1)*3 + a)
          end do
       end do
    end do
        
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
             phi_internal(:) = phi(iq2, s2, :)

             d0_istate = resolvent(ph, s2, iq2, en1_sq) 

             do j = 1, num_dof_def
                do i = 1, num_dof_def
                   D0(i, j, istate1) = D0(i, j, istate1) + &
                        *d0_istate*conjg(phi_internal(i))*phi_internal(j)
                end do
             end do
          end do
       end do
    end do
  end subroutine calculate_retarded_phonon_D0
end module Green_function

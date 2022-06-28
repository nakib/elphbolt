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

module eliashberg
  !! Module containing the procedures related to the computation of the Eliashberg
  !! spectral function a2F and the e-ph coupling factor lambda.
  
  use params, only: k8, dp
  use misc, only: exit_with_message, print_message, distribute_points, &
       demux_state, mux_vector, write2file_rank1_real, write2file_rank2_real, &
       compsimps
  use wannier_module, only: epw_wannier
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics
  use delta, only: delta_fn_tetra, delta_fn_triang

  implicit none

  private
  public calculate_a2F, calculate_iso_Matsubara_lambda, calculate_aniso_Matsubara_lambda

  external chdir

contains

  subroutine calculate_a2F(wann, el, ph, num, omegas, iso_lambda0, omegalog, external_eps_switch)
    !! Parallel driver of a2F_mk(nk'|sq) over IBZ electron states.
    !! Other zero temperature quantities such as the isotropic electron-phonon
    !! coupling (lambda0) and the log-averaged boson energy (omegalog) are
    !! also calculated here.
    !!
    !! This subroutine will calculate the anisotropic Eliashberg spectral function for
    !! all the energy window restricted electron-phonon processes for a given
    !! irreducible initial electron state = (band, wave vector). 
    !! This list will be written to disk in files tagged with the muxed state index.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.
    
    type(epw_wannier), intent(in) :: wann
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    real(dp), intent(in) :: omegas(:)
    real(dp), intent(out) :: iso_lambda0, omegalog
    logical, intent(in), optional :: external_eps_switch

    !Local variables
    integer(k8) :: nstates_irred, istate, m, ik, n, ikp, s, &
         iq, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, nprocs, num_active_images, &
         iomega, numomega
    real(dp) :: k(3), kp(3), en_el, WWp, aux, domega
    real(dp), allocatable :: ph_deltas(:, :, :), g2_istate(:), a2F_istate(:, :), &
         iso_a2F_branches(:, :), cum_iso_lambda_branches(:, :), eps_squared(:)
    character(len = 1024) :: filename
    logical :: use_external_eps = .false.

    !Number of equidistant Boson energy points
    numomega = size(omegas)

    !Boson energy difference
    domega = omegas(2) - omegas(1)
    
    !Precalculate the phonon delta functions.
    call print_message("Precalculating FBZ phonon delta functions...")
    
    ! Allocate and initialize
    allocate(ph_deltas(wann%numbranches, ph%nwv, numomega))
    ph_deltas = 0.0_dp

    call distribute_points(numomega, chunk, start, end, num_active_images)

    ! Handle tetrahedron vs triangle choice.
    if(num%tetrahedra) then !tetrahedron method
       do iomega = start, end
          do iq = 1, ph%nwv
             do s = 1, wann%numbranches
                ph_deltas(s, iq, iomega) = &
                     delta_fn_tetra(omegas(iomega), iq, s, ph%wvmesh, ph%tetramap, &
                     ph%tetracount, ph%tetra_evals)
             end do
          end do
       end do
    else !triangle method
       do iomega = start, end
          do iq = 1, ph%nwv
             do s = 1, wann%numbranches
                ph_deltas(s, iq, iomega) = &
                     delta_fn_triang(omegas(iomega), iq, s, ph%wvmesh, ph%triangmap, &
                     ph%triangcount, ph%triang_evals)
             end do
          end do
       end do
    end if

    ! The delta weights above were supercell number normalized.
    ! Taking this fact into account here:
    if(num%tetrahedra) then
       ph_deltas = ph_deltas*product(ph%wvmesh)
    else
       ph_deltas = ph_deltas*product(ph%wvmesh)
    end if

    ! Absorb DOS(Ef) in the definition of ph_deltas
    ph_deltas = ph_deltas*el%spinnormed_dos_fermi
    
    ! Reduce ph_deltas
    call co_sum(ph_deltas)

    call print_message("Calculating a2F for all IBZ electrons...")

    if(present(external_eps_switch) .and. external_eps_switch) then
       use_external_eps = .true.
       call print_message(" Externally generated screening will be used...")

       allocate(eps_squared(numomega))

       !Read |epsilon|^2 from file
       if(this_image() == 1) then
          call chdir(num%cwd)
          filename = 'eps_squared'
          open(1, file=trim(filename), status='old')
          do iomega = 1, numomega
             read(1, *) eps_squared(iomega)
          end do
          close(1)
       end if
       call co_broadcast(eps_squared, 1)
    end if
    
    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*wann%numwannbands

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

    !Allocate and initialize isotropic a2F
    allocate(iso_a2F_branches(numomega, wann%numbranches))
    iso_a2F_branches = 0.0_dp
    
    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, wann%numwannbands, m, ik)

       !Electron energy
       en_el = el%ens_irred(ik, m)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(en_el - el%enref) > el%fsthick) cycle

       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el%wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*el%wvmesh)

       !Load g2_istate from disk for a2F_istate calculation
       ! Change to data output directory
       call chdir(trim(adjustl(num%g2dir)))

       ! Read data in binary format
       write (filename, '(I9)') istate
       filename = 'gk2.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'old', access = 'stream')
       read(1) nprocs
       if(allocated(g2_istate)) deallocate(g2_istate, a2F_istate)
       allocate(g2_istate(nprocs), a2F_istate(nprocs, numomega))
       if(nprocs > 0) read(1) g2_istate
       close(1)

       ! Change back to working directory
       call chdir(num%cwd)

       !Initialize eligible process counter for this state
       count = 0

       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nwv
          !Final wave vector (crystal coords.)
          kp = el%wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*el%wvmesh)
          
          !Find interacting phonon wave vector
          ! Note that q, k, and k' are all on the same mesh
          q_indvec = modulo(kp_indvec - k_indvec, el%wvmesh)

          ! Muxed index of q
          iq = mux_vector(q_indvec, el%wvmesh, 0_k8)
          
          !Run over final electron bands
          do n = 1, wann%numwannbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle

             !Delta function contributions and k-point weight
             WWp = el%Ws_irred(ik, m)*el%Ws(ikp, n)*el%nequiv(ik)
             
             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1

                !Note that the phonon branch index iterates last for a2F_istate
                a2F_istate(count, :) = g2_istate(count)*ph_deltas(s, iq, :)

                !Sum contribuion to the isotropic a2F
                iso_a2F_branches(:, s) = iso_a2F_branches(:, s) + &
                     a2F_istate(count, :)*WWp
             end do !s
          end do !n
       end do !ikp

       !Screen, if desired
       if(use_external_eps) then
          do s = 1, wann%numbranches
             a2F_istate(:, s) = a2F_istate(:, s)/eps_squared(:)
          end do
       end if
       
       !Change to data output directory
       call chdir(trim(adjustl(num%scdir)))

       !Write data in binary format
       !Note: this will overwrite existing data!
       write (filename, '(I9)') istate
       filename = 'a2F.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'replace', access = 'stream')
       write(1) count
       write(1) a2F_istate
       close(1)

       !Change back to working directory
       call chdir(num%cwd)

       deallocate(g2_istate, a2F_istate)
    end do

    !Screen, if desired.
    if(use_external_eps) then
       do s = 1, wann%numbranches
          iso_a2F_branches(:, s) = iso_a2F_branches(:, s)/eps_squared(:)
       end do
    end if
    
    !Reduce iso_a2F_branches
    call co_sum(iso_a2F_branches)

    !Write isotropic a2F to file
    call chdir(num%cwd)
    call write2file_rank2_real('a2F_iso_branch_resolved', iso_a2F_branches)
    call write2file_rank1_real('a2F_iso', sum(iso_a2F_branches, dim = 2))

    ! Calculate cumulative lambda
    allocate(cum_iso_lambda_branches(numomega, wann%numbranches))
    do s = 1, wann%numbranches
       do iomega = 1, numomega
          call compsimps(&
               iso_a2F_branches(1:iomega, s)*2.0_dp/omegas(1:iomega), &
               domega, aux)
          cum_iso_lambda_branches(iomega, s) = aux
       end do
    end do

    ! Calculate and print out standard isotropic lambda
    iso_lambda0 = sum(cum_iso_lambda_branches(numomega, :))
    if(this_image() == 1) then
       write(*,"(A, (1E16.8))") ' Standard, isotropic e-ph coupling =', &
            iso_lambda0
    end if

    ! Calculate and print out log-averaged phonon energy
    call compsimps(&
         log(omegas(:))*sum(iso_a2F_branches, dim = 2)*2.0_dp/omegas(:), &
         domega, aux)
    omegalog = exp(aux/iso_lambda0)
    if(this_image() == 1) then
       write(*,"(A, (1E16.8, x), A)") ' Log-averaged phonon energy =', &
            omegalog*1.0e3_dp, ' meV'
    end if

    ! Print cumulative lambda to file
    call write2file_rank2_real('cum_lambda_iso_branch_resolved', cum_iso_lambda_branches)
    
    sync all
  end subroutine calculate_a2F

  subroutine calculate_iso_Matsubara_lambda(wann, num, omegas, bose_matsubara_ens, &
       iso_matsubara_lambda)
    !! Calculate the isotropic Matsubara electron-phonon coupling, lambda(l).
    !! Here l is the Bosonic Matsubara energy index. 
    
    type(epw_wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    real(dp), intent(in) :: omegas(:), bose_matsubara_ens(:)
    real(dp), intent(out) :: iso_matsubara_lambda(:)

    !Local variables
    integer(k8) :: s, l, iomega, numomega, nummatsubara
    real(dp) :: aux, domega
    real(dp), allocatable :: iso_a2F_branches(:, :)
    character(len = 1024) :: filename
    
    !Number of equidistant Boson energy points
    numomega = size(omegas)

    !Boson energy difference
    domega = omegas(2) - omegas(1)

    !Number of Bosonic Matsubara energy points
    nummatsubara = size(bose_matsubara_ens)

    call print_message("   Calculating isotropic Matsubara lambda...")

    !Read iso_a2F_branches from file
    allocate(iso_a2F_branches(numomega, wann%numbranches))
    if(this_image() == 1) then
       call chdir(num%cwd)
       filename = 'a2F_iso_branch_resolved'
       open(1, file=trim(filename), status='old')
       do iomega = 1, numomega
          read(1, *) iso_a2F_branches(iomega, :)
       end do
       close(1)
    end if
    call co_broadcast(iso_a2F_branches, 1)
    
    !Isotropic theory
    iso_matsubara_lambda = 0.0_dp

    !Sum over phonon branches
    do s = 1, wann%numbranches
       !Run over Matsubara frequencies
       do l = 1, nummatsubara
          call compsimps(&
               iso_a2F_branches(:, s)*2.0_dp*omegas(:)/ &
               (omegas(:)**2 + bose_matsubara_ens(l)**2), domega, aux)
          iso_matsubara_lambda(l) = iso_matsubara_lambda(l) + aux
       end do
    end do
  end subroutine calculate_iso_Matsubara_lambda

  subroutine calculate_aniso_Matsubara_lambda(wann, num, el, omegas, bose_matsubara_ens)
    !! Parallel driver of lambda_mk(nk'|sq, l) over IBZ electron states.
    !! Here l is the Bosonic Matsubara energy index. 
    !!
    !! This subroutine will calculate the anisotropic e-ph coupling function for
    !! all the energy window restricted electron-phonon processes for a given
    !! irreducible initial electron state = (band, wave vector). 
    !! This list will be written to disk in files tagged with the muxed state index.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.
    
    type(epw_wannier), intent(in) :: wann
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    real(dp), intent(in) :: omegas(:), bose_matsubara_ens(:)

    !Local variables
    integer(k8) :: nstates_irred, istate, m, ik, n, ikp, s, l, &
         start, end, chunk, count, nprocs, &
         num_active_images, numomega, nummatsubara
    real(dp) :: aux, domega
    real(dp), allocatable :: a2F_istate(:, :), matsubara_lambda_istate(:, :)
    character(len = 1024) :: filename
    
    call print_message("   Calculating lambda for all IBZ electrons...")
    
    !Number of equidistant Boson energy points
    numomega = size(omegas)

    !Boson energy difference
    domega = omegas(2) - omegas(1)

    !Number of Bosonic Matsubara energy points
    nummatsubara = size(bose_matsubara_ens)

    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*wann%numwannbands

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") "    #states = ", nstates_irred
       write(*, "(A, I10)") "    #states/image = ", chunk
    end if

    !Anisotropic theory
    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, wann%numwannbands, m, ik)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle

       !Load a2F_istate from disk for matsubara_lambda_istate calculation
       ! Change to data output directory
       call chdir(trim(adjustl(num%scdir)))

       ! Read data in binary format
       write (filename, '(I9)') istate
       filename = 'a2F.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'old', access = 'stream')
       read(1) nprocs
       if(allocated(a2F_istate)) deallocate(a2F_istate, matsubara_lambda_istate)
       allocate(a2F_istate(nprocs, numomega), matsubara_lambda_istate(nprocs, nummatsubara))
       if(nprocs > 0) read(1) a2F_istate
       close(1)

       ! Change back to working directory
       call chdir(num%cwd)

       !Initialize eligible process counter for this state
       count = 0

       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nwv
          !Run over final electron bands
          do n = 1, wann%numwannbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle

             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1

                !Note that the phonon branch index iterates last for a2F_istate
                !and matsubara_lambda_istate is similarly phonon branch resolved.
                !
                !Calculate lambda for each Matsubara energy
                do l = 1, nummatsubara
                   call compsimps(&
                        a2F_istate(count, :)*2.0_dp*omegas(:)/ &
                        (omegas(:)**2 + bose_matsubara_ens(l)**2), domega, aux)
                   matsubara_lambda_istate(count, l) = aux
                end do
             end do !s
          end do !n
       end do !ikp

       !Change to data output directory
       call chdir(trim(adjustl(num%scdir)))

       !Write data in binary format
       !Note: this will overwrite existing data!
       write (filename, '(I9)') istate
       filename = 'lambda.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'replace', access = 'stream')
       write(1) count
       write(1) matsubara_lambda_istate
       close(1)

       !Change back to working directory
       call chdir(num%cwd)

       deallocate(a2F_istate, matsubara_lambda_istate)
    end do
    
    sync all
  end subroutine calculate_aniso_Matsubara_lambda
end module eliashberg

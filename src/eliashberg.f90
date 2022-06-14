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

module eliashberg
  !! Module containing the procedures related to the computation of the Eliashberg
  !! spectral function a2F and the e-ph coupling factor lambda.
  
  use params, only: k8, dp
  use misc, only: exit_with_message, print_message, distribute_points, &
       demux_state, mux_vector
  use wannier_module, only: epw_wannier
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics
  use delta, only: delta_fn_tetra, delta_fn_triang

  implicit none

  private
  public calculate_a2F

contains

  subroutine calculate_a2F(wann, el, ph, num, omegas, dos_ef)
    !! Parallel driver of a2F_mk(nk'|sq) over IBZ electron states.
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
    real(dp), intent(in) :: omegas(:), dos_ef

    !Local variables
    integer(k8) :: nstates_irred, istate, m, ik, n, ikp, s, &
         iq, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, nprocs, num_active_images, &
         iomega, numomega
    real(dp) :: k(3), kp(3), en_el
    real(dp), allocatable :: ph_deltas(:, :, :), g2_istate(:), a2F_istate(:, :)
    character(len = 1024) :: filename
    
    call print_message("Calculating a2F for all IBZ electrons...")

    !Number of equidistant Boson energy points
    numomega = size(omegas)
    
    !Precalculate the phonon delta functions.
    allocate(ph_deltas(wann%numbranches, ph%nwv, numomega))
    do iomega = 1, numomega
       do iq = 1, ph%nwv
          do s = 1, wann%numbranches
             if(num%tetrahedra) then
                ph_deltas(s, iq, iomega) = &
                     delta_fn_tetra(omegas(iomega), iq, s, ph%wvmesh, ph%tetramap, &
                     ph%tetracount, ph%tetra_evals)
             else
                ph_deltas(s, iq, iomega) = &
                     delta_fn_triang(omegas(iomega), iq, s, ph%wvmesh, ph%triangmap, &
                     ph%triangcount, ph%triang_evals)
             end if
          end do
       end do
    end do

    !Number of scattering processes within the transport window
    nprocs = el%nstates_inwindow*wann%numbranches

    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*wann%numwannbands

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

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
          q_indvec = kp_indvec - k_indvec

          ! Muxed index of q
          iq = mux_vector(q_indvec, ph%wvmesh, 0_k8)
          
          !Run over final electron bands
          do n = 1, wann%numwannbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle
             
             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1

                !Note that the phonon branch index iterates last for a2F_istate
                a2F_istate(count, :) = g2_istate(count)*ph_deltas(s, iq, :)
             end do !s
          end do !n
       end do !ikp

       !Multiply with electron DOS at Fermi level
       a2F_istate = a2F_istate*dos_ef

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
    sync all
  end subroutine calculate_a2F
  
end module eliashberg

module interactions
  !! Module containing the procedures related to the computation of interactions.

  use params, only: k4, dp
  use misc, only: print_message, distribute_points, demux_state, mux_vector
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics

  implicit none

  private
  public calculate_g_mixed, calculate_g_bloch
  
contains

  subroutine calculate_g_mixed(wann, el, num)
    !! Parallel driver of gmixed_epw over IBZ electron wave vectors.

    type(epw_wannier), intent(in) :: wann
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num

    !Local variables
    integer(k4) :: ik, ikstart, ikend, chunk, ierr

    call print_message("Doing g(Re,Rp) -> g(k,Rp) for all IBZ k...")

    call distribute_points(el%nk_irred, chunk, ikstart, ikend)

    if(this_image() == 1) then
       print*, "   #k = ", el%nk_irred
       print*, "   #k/image = ", chunk
    end if

    do ik = ikstart, ikend
       call wann%gmixed_epw(num, ik, el%wavevecs_irred(ik, :))
    end do

    sync all
  end subroutine calculate_g_mixed

  subroutine calculate_g_bloch(wann, crys, el, ph, num)
    !! Parallel driver of g2_epw over IBZ electron states within the Fermi window.
    !! This subroutine will calculate the full Bloch rep. matrix elements for
    !! all the energy window restricted electron-phonon processes for a given
    !! irreducible initial electron state = (band, wave vector). 
    !! This list will be written to disk in files tagged with the muxed state index.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside energy window
    !will be skipped in the calculation as they are irrelevant for transport.

    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(k4) :: nstates_irred, istate, m, ik, ik_muxed, n, ikp, s, &
         iq, start, end, chunk, ierr, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, g2size
    real(dp) :: k(3), kp(3), q(3)
    real(dp), allocatable :: g2_istate(:)
    complex(dp) :: gmixed_ik(wann%numwannbands,wann%numwannbands,wann%numbranches,wann%nwsq)
    character(len = 1024) :: filename

    call print_message("Doing g(k,Rp) -> |g(k,q)|^2 for all IBZ states...")

    !Length of g2_istate
    g2size = el%nstates_inwindow*wann%numbranches
    allocate(g2_istate(g2size))

    !Total number of IBZ blocks states
    nstates_irred = el%nk_irred*wann%numwannbands

    call distribute_points(nstates_irred, chunk, start, end)

    if(this_image() == 1) then
       print*, "   #states = ", nstates_irred
       print*, "   #states/image = ", chunk
    end if

    count = 0
    do istate = start, end !over IBZ blocks states
       !Initialize eligible process counter for this state
       count = 0

       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate,wann%numwannbands, m, ik)

       !Get the muxed index of wave vector from the IBZ blocks index list
       ik_muxed = el%indexlist_irred(ik)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle

       !Load gmixed(ik) here for use inside the loops below
       call chdir(trim(adjustl(num%g2dir)))
       write (filename, '(I6)') ik
       filename = 'gmixed.ik'//trim(adjustl(filename))
       open(1,file=filename,status="old",access='stream')
       read(1) gmixed_ik
       close(1)
       call chdir(num%cwd)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el%wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*el%kmesh)

       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nk
          !Final wave vector (crystal coords.)
          kp = el%wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*el%kmesh)

          !Run over final electron bands
          do n = 1, wann%numwannbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle

             !Find interacting phonon wave vector
             !Note that q, k, and k' are all on the same mesh
             q_indvec = modulo(kp_indvec - k_indvec, el%kmesh) !0-based index vector
             q = q_indvec/dble(el%kmesh) !crystal coords.

             !Muxed index of q
             iq = mux_vector(q_indvec, el%kmesh, 0_dp)

             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1

                !Calculate |g_mns(<k>,q)|^2
                g2_istate(count) = wann%g2_epw(crys, q, el%evecs_irred(ik, m, :), &
                     el%evecs(ikp, n, :), ph%evecs(iq, s, :), &
                     ph%ens(iq, s), gmixed_ik)
             end do !s
          end do !n
       end do !ikp

       !Change to data output directory
       call chdir(trim(adjustl(num%g2dir)))

       !Write data in binary format
       !Note: this will overwrite existing data!
       write (filename, '(I9)') istate
       filename = 'g2.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'replace', access = 'stream')
       write(1) g2_istate
       close(1)
       
       !Change back to working directory
       call chdir(num%cwd)
    end do

    sync all
    
    !Delete the gmixed disk data
    if(this_image() == 1) then
       call chdir(trim(adjustl(num%g2dir)))
       call system('rm gmixed.*')
       call chdir(num%cwd)
    endif
  end subroutine calculate_g_bloch
end module interactions

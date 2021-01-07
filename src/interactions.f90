module interactions
  !! Module containing the procedures related to the computation of interactions.

  use params, only: k4, dp
  use misc, only: print_message, distribute_points, demux_state, mux_vector
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon, Vm2_3ph
  use numerics_module, only: numerics

  implicit none

  private
  public calculate_g_mixed, calculate_g2_bloch, calculate_V_3ph
  
contains

  subroutine calculate_V_3ph(ph, num)
    !! Parallel driver of the 3-ph vertex calculator for all IBZ phonon wave vectors.
    !! This subroutine calculates a subset of the the upper triangular (phonon 2, phonon 3)
    !! subspace of |V-(<q1>|q2,q3)|^2 for each irreducible phonon 1 and save the result to disk.

    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(k4) :: start, end, chunk, istate, nstates_irred, nstates_triangle, &
         s, iq, iq_muxed

    call print_message("Calculating minimal set of |V-(<q1>|q2,q3)|^2...")

    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    !Size of upper triangle of phonon 2 and phonon 3 wave vectors
    nstates_triangle = ph%nq*ph%numbranches**2

    allocate(Vm2(nstates_triangle))

    call distribute_points(nstates_irred, chunk, start, end)

    if(this_image() == 1) then
       print*, "   #states = ", nstates_irred
       print*, "   #states/image = ", chunk
    end if

    !Run over first phonon IBZ states
    do istate1 = start, end
       !Initialize process counter for this state
       count = 0

       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate1, ph%numbranches, s1, iq1)

       !Energy of phonon 1
       en1 = ph%ens_irred(iq, s)
       
       !Get the muxed index of wave vector from the IBZ index list
       iq1_muxed = el%indexlist_irred(iq1)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs_irred(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%qmesh)

       !Run over second (FBZ) phonon wave vectors
       do iq2 = 1, ph%nq
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%qmesh)

          !Folded final phonon wave vector
          q3_indvec = modulo(q1_indvec - q2_indvec, ph%qmesh) !0-based index vector
          q3 = q3_indvec/dble(ph%qmesh) !crystal coords.

          !Use q2 <-> q3 exchange symmetry of V-(q1|q2,q3)
          if(all(q3_indvec(:) < q2_indvec(:))) cycle
          
          !Run over branches of second phonon
          do s2 = 1, ph%numbranches
             !Energy of phonon 2
             en2 = ph%ens_irred(iq2, s2)

             !Run over branches of third phonon
             do s3 = 1, ph%numbranches
                !Energy of phonon 3
                en3 = ph%ens_irred(iq3, s3)

                !Count process
                count = count + 1

                !Calculate the minus process vertex
                Vm2(count) = Vm2_3ph(ph, s1, s2, s3, iq1, iq2, iq3)
             end do !s3
          end do !s2
       end do !iq2
    end do !istate1

    sync all
  end subroutine calculate_V_3ph

  subroutine calculate_g_mixed(wann, el, num)
    !! Parallel driver of gmixed_epw over IBZ electron wave vectors.

    type(epw_wannier), intent(in) :: wann
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num

    !Local variables
    integer(k4) :: ik, ikstart, ikend, chunk

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

  subroutine calculate_g2_bloch(wann, crys, el, ph, num)
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
         iq, start, end, chunk, k_indvec(3), kp_indvec(3), &
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
  end subroutine calculate_g2_bloch
end module interactions

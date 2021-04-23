module interactions
  !! Module containing the procedures related to the computation of interactions.

  use params, only: k4, dp, pi, amu, qe, hbar_eVps
  use misc, only: exit_with_message, print_message, distribute_points, &
       demux_state, mux_vector, mux_state, expi, Bose
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics
  use delta, only: delta_fn_tetra

  implicit none

  private
  public calculate_g_mixed, calculate_g2_bloch, calculate_3ph_interaction, &
       calculate_ph_rta_rates, read_transition_probabilities

contains
  
  pure real(dp) function Vm2_3ph(s1, s2, s3, ev1, ev2, ev3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3)
    !! Function to calculate the 3-ph interaction vertex |V-|^2.
    
    integer(k4), intent(in) :: s1, s2, s3, Index_i(:), Index_j(:), Index_k(:)
    complex(dp), intent(in) :: phases_q2q3(:), ev1(:,:), ev2(:,:), ev3(:,:)
    real(dp), intent(in) :: ifc3(:,:,:,:)

    !Local variables
    integer(k4) :: it, a, b, c, aind, bind, cind, ntrip
    complex(dp) :: aux1, aux2, V0

    ntrip = size(Index_k(:))
    
    aux1 = (0.0_dp, 0.0_dp)
    do it = 1, ntrip
       V0 = (0.0_dp, 0.0_dp)
       do a = 1, 3
          aind = a + 3*(Index_k(it) - 1)
          do b = 1, 3
             bind = b + 3*(Index_j(it) - 1)
             aux2 = conjg(ev2(s2, bind)*ev3(s3, aind))
             do c = 1, 3
                cind = c + 3*(Index_i(it) - 1)
                V0 = V0 + ifc3(c, b, a, it)*ev1(s1, cind)*aux2
             end do
          end do
       end do
       aux1 = aux1 + V0*phases_q2q3(it)
    end do

    Vm2_3ph = abs(aux1)**2
  end function Vm2_3ph
  
  subroutine calculate_3ph_interaction(ph, crys, num, key)
    !! Parallel driver of the 3-ph vertex calculator for all IBZ phonon wave vectors.
    !! This subroutine calculates |V-(s1<q1>|s2q2,s3q3)|^2, W-(s1<q1>|s2q2,s3q3),
    !! and W+(s1<q1>|s2q2,s3q3) for each irreducible phonon and saves the results to disk.
    !!
    !! key = 'V', 'W' for vertex, scattering rate calculation, respectively.

    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    character(len = 1), intent(in) :: key
    
    !Local variables
    integer(k4) :: start, end, chunk, istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), index_minus, index_plus, &
         neg_iq2, neg_q2_indvec(3)
    real(dp) :: en1, en2, en3, massfac, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         delta, occup_fac, Vp2_index_plus, const, bose2, bose3
    real(dp), allocatable :: Vm2(:), Wm(:), Wp(:)
    integer(k4), allocatable :: istate2_plus(:), istate3_plus(:), istate2_minus(:), istate3_minus(:)
    complex(dp) :: phases_q2q3(ph%numtriplets)
    character(len = 1024) :: filename, filename_Wm, filename_Wp, Wdir, tag

    if(key /= 'V' .and. key /= 'W') then
       call exit_with_message("Invalid value of key in call to calculate_3ph_interaction. Exiting.")
    end if

    if(key == 'V') then
       call print_message("Calculating 3-ph vertices for all IBZ q...")
    else
       call print_message("Calculating 3-ph scattering rates for all IBZ q...")
    end if

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Wdir = trim(adjustl(num%datadumpdir))//'W_T'//trim(adjustl(tag))
    if(this_image() == 1 .and. key == 'W') then
       !Create directory
       call system('mkdir '//trim(adjustl(Wdir)))
    end if
    sync all
   
    !Conversion factor in transition probability expression
    const = pi/4.0_dp*hbar_eVps**5*(qe/amu)**3*1.0d-12

    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    !Total number of 3-phonon processes for a given initial phonon state
    nprocs = ph%nq*ph%numbranches**2

    !Allocate |V^-|^2 and, if needed, W- and W+
    allocate(Vm2(nprocs))
    if(key == 'W') then
       allocate(Wp(nprocs), Wm(nprocs))
       allocate(istate2_plus(nprocs), istate3_plus(nprocs),&
            istate2_minus(nprocs),istate3_minus(nprocs))
    end if

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end)

    if(this_image() == 1) then
       print*, "   #states = ", nstates_irred
       print*, "   #states/image = ", chunk
    end if

    !Run over first phonon IBZ states
    do istate1 = start, end
       !Initialize transition probabilities
       if(key == 'W') then
          Wp(:) = 0.0_dp
          Wm(:) = 0.0_dp
          istate2_plus(:) = 0_k4
          istate3_plus(:) = 0_k4
          istate2_minus(:) = 0_k4
          istate3_minus(:) = 0_k4
       end if

       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate1, ph%numbranches, s1, iq1_ibz)

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Energy of phonon 1
       en1 = ph%ens(iq1, s1)
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%qmesh)

       !Load |V^-|^2 from disk for scattering rates calculation
       if(key == 'W') then
          !Change to data output directory
          call chdir(trim(adjustl(num%Vdir)))

          !Read data in binary format
          write (filename, '(I9)') istate1
          filename = 'Vm2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')
          read(1) Vm2
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
       end if
       
       !Run over second (FBZ) phonon wave vectors
       do iq2 = 1, ph%nq
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%qmesh)

          !Folded final phonon wave vector
          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%qmesh) !0-based index vector
          q3_minus = q3_minus_indvec/dble(ph%qmesh) !crystal coords.

          !Muxed index of q3_minus
          iq3_minus = mux_vector(q3_minus_indvec, ph%qmesh, 0_k4)
          
          if(key == 'V') then
             !Calculate the numtriplet number of mass-normalized phases for this (q2,q3) pair
             do it = 1, ph%numtriplets
                massfac = 1.0_dp/sqrt(&
                     crys%masses(crys%atomtypes(ph%Index_i(it)))*&
                     crys%masses(crys%atomtypes(ph%Index_j(it)))*&
                     crys%masses(crys%atomtypes(ph%Index_k(it))))
                q2_cart = matmul(crys%reclattvecs, q2)
                q3_minus_cart = matmul(crys%reclattvecs, q3_minus)
                phases_q2q3(it) = massfac*&
                     expi(-dot_product(q2_cart, ph%R_j(:,it)) -&
                     dot_product(q3_minus_cart, ph%R_k(:,it)))
             end do
          end if
          
          !Run over branches of second phonon
          do s2 = 1, ph%numbranches
             if(key == 'W') then
                !Energy of phonon 2
                en2 = ph%ens(iq2, s2)

                !Bose factor for phonon 2
                bose2 = Bose(en2, crys%T)

                !Get index of -q2
                neg_q2_indvec = modulo(-q2_indvec, ph%qmesh)
                neg_iq2 = mux_vector(neg_q2_indvec, ph%qmesh, 0_k4)
             end if
             
             !Run over branches of third phonon
             do s3 = 1, ph%numbranches                
                !Minus process index
                index_minus = ((iq2 - 1)*ph%numbranches + (s2 - 1))*ph%numbranches + s3

                if(key == 'V') then
                   !Calculate the minus process vertex
                   Vm2(index_minus) = Vm2_3ph(s1, s2, s3, ph%evecs(iq1,:,:), &
                        ph%evecs(iq2,:,:), ph%evecs(iq3_minus,:,:), ph%Index_i(:), &
                        ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), phases_q2q3)
                end if

                if(key == 'W') then
                   !Energy of phonon 3
                   en3 = ph%ens(iq3_minus, s3)

                   !Bose factor for phonon 3
                   bose3 = Bose(en3, crys%T)

                   !Calculate W-:
                   
                   !Evaluate delta function
                   delta = delta_fn_tetra(en1 - en3, iq2, s2, ph%qmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals)
                   
                   !Temperature dependent occupation factor
                   !(bose1 + 1)*bose2*bose3/(bose1*(bose1 + 1))
                   ! = (bose2 + bose3 + 1)
                   occup_fac = (bose2 + bose3 + 1.0_dp)

                   !Save W-
                   if(en1*en2*en3 /= 0.0_dp) then
                      Wm(index_minus) = Vm2(index_minus)*occup_fac*delta/en1/en2/en3
                   end if

                   !Calculate W+:

                   !Grab corresponding plus process using
                   !V-(s1q1|s2q2,s3q3) = V+(s1q1|s2-q2,s3q3)
                   Vp2_index_plus = Vm2(index_minus)
                   index_plus = ((neg_iq2 - 1)*ph%numbranches + (s2 - 1))*ph%numbranches + s3
                   
                   !Evaluate delta function
                   delta = delta_fn_tetra(en3 - en1, neg_iq2, s2, ph%qmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals)

                   !Temperature dependent occupation factor
                   !(bose1 + 1)*(bose2 + 1)*bose3/(bose1*(bose1 + 1))
                   ! = bose2 - bose3.
                   occup_fac = (bose2 - bose3)
                   
                   !Save W+
                   if(en1*en2*en3 /= 0.0_dp) then
                      Wp(index_plus) = Vp2_index_plus*occup_fac*delta/en1/en2/en3
                   end if
                   
                   !Save 2nd and 3rd phonon states
                   istate2_minus(index_minus) = mux_state(ph%numbranches, s2, iq2)
                   istate2_plus(index_plus) = mux_state(ph%numbranches, s2, neg_iq2)
                   istate3_minus(index_minus) = mux_state(ph%numbranches, s3, iq3_minus)
                   istate3_plus(index_plus) = istate3_minus(index_minus)
                end if
             end do !s3
          end do !s2
       end do !iq2

       if(key == 'V') then
          !Change to data output directory
          call chdir(trim(adjustl(num%Vdir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate1
          filename = 'Vm2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) Vm2
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
       end if

       if(key == 'W') then
          !Multiply constant factor, unit factor, etc.
          Wm(:) = const*Wm(:) !THz
          Wp(:) = const*Wp(:) !THz

          !Write W+ and W- to disk
          !Change to data output directory
          call chdir(trim(adjustl(Wdir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate1

          filename_Wm = 'Wm.istate'//trim(adjustl(filename))
          open(1, file = trim(filename_Wm), status = 'replace', access = 'stream')
          write(1) Wm
          write(1) istate2_minus
          write(1) istate3_minus
          close(1)

          filename_Wp = 'Wp.istate'//trim(adjustl(filename))
          open(1, file = trim(filename_Wp), status = 'replace', access = 'stream')
          write(1) Wp
          write(1) istate2_plus
          write(1) istate3_plus
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
       end if
    end do !istate1
    sync all
  end subroutine calculate_3ph_interaction

  subroutine calculate_ph_rta_rates(ph, num, crys, rta_rates)
    !! Subroutine for parallel reading of the 3-ph transition probabilities
    !! from disk and calculating the relaxation time approximation (RTA)
    !! scattering rates for phonons.

    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    real(dp), allocatable, intent(out) :: rta_rates(:,:)

    !Local variables
    integer(k4) :: nstates_irred, procs, istate, nprocs, iproc, chunk, s, &
         iq, im
    integer(k4), allocatable :: start[:], end[:]
    real(dp), allocatable :: rta_rates_psum(:,:)[:], W(:)
    character(len = 1024) :: filepath_Wm, filepath_Wp, Wdir, tag
    
    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Wdir = trim(adjustl(num%datadumpdir))//'W_T'//trim(adjustl(tag))
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    !Total number of 3-phonon processes for a given phonon state
    nprocs = ph%nq*ph%numbranches**2

    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end)

    !Allocate and initialize scattering rates coarrays
    allocate(rta_rates_psum(ph%nq_irred, ph%numbranches)[*])
    rta_rates_psum(:, :) = 0.0_dp
    
    !Allocate and initialize scattering rates
    allocate(rta_rates(ph%nq_irred, ph%numbranches))
    rta_rates(:, :) = 0.0_dp
    
    !Allocate transition probabilities variable
    allocate(W(nprocs))
    
    !Run over first phonon IBZ states
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbranches, s, iq)

       !Set W+ filename
       write(tag, '(I9)') istate
       filepath_Wp = trim(adjustl(Wdir))//'/Wp.istate'//trim(adjustl(tag))
       
       !Read W+ from file
       call read_transition_probabilities(trim(adjustl(filepath_Wp)), W)

       do iproc = 1, nprocs
          rta_rates_psum(iq, s) = rta_rates_psum(iq, s) + W(iproc) 
       end do

       !Set W- filename
       filepath_Wm = trim(adjustl(Wdir))//'/Wm.istate'//trim(adjustl(tag))
       
       !Read W- from file
       call read_transition_probabilities(trim(adjustl(filepath_Wm)), W)

       do iproc = 1, nprocs
          rta_rates_psum(iq, s) = rta_rates_psum(iq, s) + 0.5_dp*W(iproc)
       end do
    end do

    sync all

    !Reduce coarray partial sums
    do im = 1, num_images()
       rta_rates(:,:) = rta_rates(:,:) + rta_rates_psum(:,:)[im]
    end do
    sync all
  end subroutine calculate_ph_rta_rates

  subroutine read_transition_probabilities(filepath, W, istate2, istate3)
    !! Subroutine to read the phonon transition probabilities from disk.

    character(len = *), intent(in) :: filepath
    real(dp), intent(out) :: W(:)
    integer(k4), optional, intent(out) :: istate2(:), istate3(:)

    !Read data
    open(1, file = trim(adjustl(filepath)), status = 'old', access = 'stream')
    read(1) W
    if(present(istate2) .and. present(istate3)) then
       read(1) istate2
       read(1) istate3
    end if
    close(1)
  end subroutine read_transition_probabilities
  
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
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.

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

             !TODO Here need to 
             
             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1
                
                !TODO Calculate |g_mns(<k>,q)|^2
!!$                g2_istate(count) = wann%g2_epw(crys, q, el%evecs_irred(ik, m, :), &
!!$                     el%evecs(ikp, n, :), ph%evecs(iq, s, :), &
!!$                     ph%ens(iq, s), gmixed_ik)
                !NOTE: The above will crash for unequal k and q meshes.
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
    sync all
  end subroutine calculate_g2_bloch
end module interactions

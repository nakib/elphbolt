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

module interactions
  !! Module containing the procedures related to the computation of interactions.

  use params, only: k8, dp, pi, twopi, amu, qe, hbar_eVps
  use misc, only: exit_with_message, print_message, distribute_points, &
       demux_state, mux_vector, mux_state, expi, Bose, binsearch, Fermi
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics
  use delta, only: delta_fn_tetra

  implicit none

  private
  public calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_ph_rta_rates, read_transition_probs_3ph, read_transition_probs_eph, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk, calculate_el_rta_rates

contains
  
  pure real(dp) function Vm2_3ph(ev1_s1, conjg_ev2_s2, conjg_ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
    !! Function to calculate the 3-ph interaction vertex |V-|^2.
    
    integer(k8), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(dp), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), conjg_ev2_s2(nb), conjg_ev3_s3(nb)
    real(dp), intent(in) :: ifc3(3, 3, 3, ntrip)

    !Local variables
    integer(k8) :: it, a, b, c, aind, bind, cind
    complex(dp) :: aux1, aux2, aux3, V0
    
    aux1 = (0.0_dp, 0.0_dp)
    do it = 1, ntrip
       aind = 3*(Index_k(it) - 1)
       bind = 3*(Index_j(it) - 1)
       cind = 3*(Index_i(it) - 1)
       V0 = (0.0_dp, 0.0_dp)
       do a = 1, 3
          aux2 = conjg_ev3_s3(a + aind)
          do b = 1, 3
             aux3 = aux2*conjg_ev2_s2(b + bind)
             do c = 1, 3
                V0 = V0 + ifc3(c, b, a, it)*ev1_s1(c + cind)*aux3
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
    !! key = 'V', 'W' for vertex, transition probabilitiy calculation, respectively.

    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    character(len = 1), intent(in) :: key
    
    !Local variables
    integer(k8) :: start, end, chunk, istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), index_minus, index_plus, &
         neg_iq2, neg_q2_indvec(3), num_active_images
    real(dp) :: en1, en2, en3, massfac, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         occup_fac, Vp2_index_plus, const, bose2, bose3, delta_minus, delta_plus
    real(dp), allocatable :: Vm2(:), Wm(:), Wp(:)
    integer(k8), allocatable :: istate2_plus(:), istate3_plus(:), istate2_minus(:), istate3_minus(:)
    complex(dp) :: phases_q2q3(ph%numtriplets)
    character(len = 1024) :: filename, filename_Wm, filename_Wp, Wdir, tag

    if(key /= 'V' .and. key /= 'W') then
       call exit_with_message("Invalid value of key in call to calculate_3ph_interaction. Exiting.")
    end if

    if(key == 'V') then
       call print_message("Calculating 3-ph vertices for all IBZ phonons...")
    else
       call print_message("Calculating 3-ph transition probabilities for all IBZ phonons...")
    end if

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Wdir = trim(adjustl(num%datadumpdir))//'W_T'//trim(adjustl(tag))
    if(this_image() == 1 .and. key == 'W') then
       !Create directory
       call system('mkdir -p '//trim(adjustl(Wdir)))
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
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

    !Run over first phonon IBZ states
    do istate1 = start, end
       !Initialize transition probabilities
       if(key == 'W') then
          Wp(:) = 0.0_dp
          Wm(:) = 0.0_dp
          istate2_plus(:) = 0_k8
          istate3_plus(:) = 0_k8
          istate2_minus(:) = 0_k8
          istate3_minus(:) = 0_k8
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
          iq3_minus = mux_vector(q3_minus_indvec, ph%qmesh, 0_k8)
          
          if(key == 'V') then
             if(en1 /= 0.0_dp) then
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
          end if
          
          !Run over branches of second phonon
          do s2 = 1, ph%numbranches
             !Energy of phonon 2
             en2 = ph%ens(iq2, s2)

             !Get index of -q2
             neg_q2_indvec = modulo(-q2_indvec, ph%qmesh)
             neg_iq2 = mux_vector(neg_q2_indvec, ph%qmesh, 0_k8)

             if(key == 'W') then
                !Bose factor for phonon 2
                bose2 = Bose(en2, crys%T)
             end if
             
             !Run over branches of third phonon
             do s3 = 1, ph%numbranches                
                !Minus process index
                index_minus = ((iq2 - 1)*ph%numbranches + (s2 - 1))*ph%numbranches + s3

                !Energy of phonon 3
                en3 = ph%ens(iq3_minus, s3)

                !Evaluate delta functions
                delta_minus = delta_fn_tetra(en1 - en3, iq2, s2, ph%qmesh, ph%tetramap, &
                     ph%tetracount, ph%tetra_evals)

                delta_plus = delta_fn_tetra(en3 - en1, neg_iq2, s2, ph%qmesh, ph%tetramap, &
                     ph%tetracount, ph%tetra_evals)
                
                if(key == 'V') then
                   if(en1*en2*en3 == 0.0_dp .or. &
                        (delta_minus == 0.0_dp .and. delta_plus == 0.0_dp)) then
                      Vm2(index_minus) = 0.0_dp
                   else
                      !Calculate the minus process vertex
                      Vm2(index_minus) = Vm2_3ph(ph%evecs(iq1, s1, :), &
                           conjg(ph%evecs(iq2, s2, :)), conjg(ph%evecs(iq3_minus, s3, :)), &
                           ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
                           phases_q2q3, ph%numtriplets, ph%numbranches)
                   end if
                end if

                if(key == 'W') then
                   !Bose factor for phonon 3
                   bose3 = Bose(en3, crys%T)

                   !Calculate W-:
                   
                   !Temperature dependent occupation factor
                   !(bose1 + 1)*bose2*bose3/(bose1*(bose1 + 1))
                   ! = (bose2 + bose3 + 1)
                   occup_fac = (bose2 + bose3 + 1.0_dp)

                   !Save W-
                   if(en1*en2*en3 /= 0.0_dp) then
                      Wm(index_minus) = Vm2(index_minus)*occup_fac*delta_minus/en1/en2/en3
                   end if

                   !Calculate W+:

                   !Grab corresponding plus process using
                   !V-(s1q1|s2q2,s3q3) = V+(s1q1|s2-q2,s3q3)
                   Vp2_index_plus = Vm2(index_minus)
                   index_plus = ((neg_iq2 - 1)*ph%numbranches + (s2 - 1))*ph%numbranches + s3

                   !Temperature dependent occupation factor
                   !(bose1 + 1)*(bose2 + 1)*bose3/(bose1*(bose1 + 1))
                   ! = bose2 - bose3.
                   occup_fac = (bose2 - bose3)
                   
                   !Save W+
                   if(en1*en2*en3 /= 0.0_dp) then
                      Wp(index_plus) = Vp2_index_plus*occup_fac*delta_plus/en1/en2/en3
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
       end if

       !Change back to working directory
       call chdir(num%cwd)
    end do !istate1
    sync all
  end subroutine calculate_3ph_interaction

  subroutine calculate_gReq(wann, ph, num)
    !! Parallel driver of gReq_epw over IBZ phonon wave vectors.

    type(epw_wannier), intent(in) :: wann
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num

    !Local variables
    integer(k8) :: iq, iqstart, iqend, chunk, num_active_images

    call print_message("Calculating g(Re,Rp) -> g(Re,q) for all IBZ q...")

    call distribute_points(ph%nq_irred, chunk, iqstart, iqend, num_active_images)

    if(this_image() == 1) then
       print*, "   #q = ", ph%nq_irred
       print*, "   #q/image = ", chunk
    end if

    do iq = iqstart, iqend
       call wann%gReq_epw(num, iq, ph%wavevecs_irred(iq, :))
    end do

    sync all
  end subroutine calculate_gReq
  
  subroutine calculate_gkRp(wann, el, num)
    !! Parallel driver of gkRp_epw over IBZ electron wave vectors.

    type(epw_wannier), intent(in) :: wann
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num

    !Local variables
    integer(k8) :: ik, ikstart, ikend, chunk, num_active_images

    call print_message("Calculating g(Re,Rp) -> g(k,Rp) for all IBZ k...")

    call distribute_points(el%nk_irred, chunk, ikstart, ikend, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #k = ", el%nk_irred
       write(*, "(A, I10)") " #k/image = ", chunk
    end if

    do ik = ikstart, ikend
       call wann%gkRp_epw(num, ik, el%wavevecs_irred(ik, :))
    end do

    sync all
  end subroutine calculate_gkRp

  subroutine calculate_eph_interaction_ibzq(wann, crys, el, ph, num, key)
    !! Parallel driver of g2(q,k) over IBZ phonon states.
    !!
    !! This subroutine will calculate the full Bloch rep. matrix elements for
    !! all the energy window restricted electron-phonon processes for a given
    !! irreducible initial phonon state = (branch, wave vector). 
    !! This list will be written to disk in files tagged with the muxed state index.
    !!
    !! key = 'g', 'Y' for vertex, transition probability calculation, respectively.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.

    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    character(len = 1), intent(in) :: key
    
    !Local variables
    integer(k8) :: nstates_irred, istate, m, iq, iq_fbz, n, ik, ikp, s, &
         ikp_window, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), nprocs, count, num_active_images
    integer(k8), allocatable :: istate1(:), istate2(:)
    real(dp) :: k(3), kp(3), q(3), en_ph, en_el, en_elp, const, delta, &
         invboseplus1, fermi1, fermi2, occup_fac
    real(dp), allocatable :: g2_istate(:), Y_istate(:)
    complex(dp) :: gReq_iq(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsk)
    character(len = 1024) :: filename, filename_Y, Ydir, tag

    if(key /= 'g' .and. key /= 'Y') then
       call exit_with_message(&
            "Invalid value of key in call to calculate_eph_interaction_ibzq. Exiting.")
    end if

    if(key == 'g') then
       call print_message("Calculating g(Re,q) -> |g(k,q)|^2 for all IBZ phonons...")
    else
       call print_message("Calculating ph-e transition probabilities for all IBZ phonons...")
    end if
    
    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Ydir = trim(adjustl(num%datadumpdir))//'Y_T'//trim(adjustl(tag))
    if(this_image() == 1 .and. key == 'Y') then
       !Create directory
       call system('mkdir -p '//trim(adjustl(Ydir)))
    end if
    sync all
    
    !Allocate and initialize g2_istate
    if(key == 'g') then
       !Maximum length of g2_istate
       nprocs = el%nstates_inwindow*ph%numbranches
       allocate(g2_istate(nprocs))
       g2_istate(:) = 0.0_dp
    end if

    !Conversion factor in transition probability expression
    const = twopi/hbar_eVps
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

    do istate = start, end !over IBZ blocks states
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbranches, s, iq)

       if(key == 'g') then
          !Load gReq(iq) here for use inside the loops below
          call chdir(trim(adjustl(num%g2dir)))
          write (filename, '(I6)') iq
          filename = 'gReq.iq'//trim(adjustl(filename))
          open(1,file=filename,status="old",access='stream')
          read(1) gReq_iq
          close(1)
          call chdir(num%cwd)
       end if

       !Get the muxed index of FBZ wave vector from the IBZ blocks index list
       iq_fbz = ph%indexlist_irred(iq)

       !Energy of phonon
       en_ph = ph%ens(iq_fbz, s)

       !1/(1 + Bose factor) for phonon
       if(key == 'Y') then
          if(en_ph /= 0.0_dp) then
             invboseplus1 = 1.0_dp/(1.0_dp + Bose(en_ph, crys%T))
          else
             invboseplus1 = 0.0_dp
          end if
       end if
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       q = ph%wavevecs(iq_fbz, :)

       !Convert from crystal to 0-based index vector
       q_indvec = nint(q*ph%qmesh)

       !Load g2_istate from disk for scattering rates calculation
       if(key == 'Y') then
          !Change to data output directory
          call chdir(trim(adjustl(num%g2dir)))

          !Read data in binary format
          write (filename, '(I9)') istate
          filename = 'gq2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')
          read(1) nprocs
          if(allocated(g2_istate)) deallocate(g2_istate, Y_istate, istate1, istate2)
          allocate(g2_istate(nprocs))
          read(1) g2_istate
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
          
          !Allocate and initialize quantities related to transition probabilities
          allocate(Y_istate(nprocs))
          allocate(istate1(nprocs), istate2(nprocs))
          istate1(:) = -1_k8
          istate2(:) = -1_k8
          Y_istate(:) = 0.0_dp
       end if

       !Initialize process counter
       count = 0
       
       !Run over initial (in-window, FBZ blocks) electron wave vectors
       do ik = 1, el%nk
          !Initial wave vector (crystal coords.)
          k = el%wavevecs(ik, :)

          !Convert from crystal to 0-based index vector
          k_indvec = nint(k*el%kmesh)
          
          !Find final electron wave vector
          kp_indvec = modulo(k_indvec + el%mesh_ref*q_indvec, el%kmesh) !0-based index vector

          !Muxed index of kp
          ikp = mux_vector(kp_indvec, el%kmesh, 0_k8)

          !Check if final electron wave vector is within energy window
          call binsearch(el%indexlist, ikp, ikp_window)
          if(ikp_window < 0) cycle
          
          !Run over initial electron bands
          do m = 1, el%numbands
             !Energy of initial electron
             en_el = el%ens(ik, m)
             
             !Apply energy window to initial electron
             if(abs(en_el - el%enref) > el%fsthick) cycle

             !Fermi factor for initial and final electrons
             if(key == 'Y') then
                fermi1 = Fermi(en_el, el%chempot, crys%T)
                fermi2 = Fermi(en_el + en_ph, el%chempot, crys%T)
             end if
             
             !Run over final electron bands
             do n = 1, el%numbands
                !Energy of final electron
                en_elp = el%ens(ikp_window, n)
                
                !Apply energy window to final electron
                if(abs(en_elp - el%enref) > el%fsthick) cycle
                
                !Increment g2 process counter
                count = count + 1

                if(key == 'g') then
                   !Calculate |g_mns(k,<q>)|^2
                   g2_istate(count) = wann%g2_epw(crys, k, q, el%evecs(ik, m, :), &
                        el%evecs(ikp_window, n, :), ph%evecs(iq_fbz, s, :), &
                        ph%ens(iq_fbz, s), gReq_iq, 'el')
                end if
                
                if(key == 'Y') then                   
                   !Calculate Y:

                   !Evaluate delta function
                   delta = delta_fn_tetra(en_elp - en_ph, ik, m, el%kmesh, el%tetramap, &
                        el%tetracount, el%tetra_evals)

                   !Temperature dependent occupation factor
                   occup_fac = fermi1*(1.0_dp - fermi2)*invboseplus1
                   
                   !Save Y
                   if(en_ph >= 0.5e-3) then !Use a small phonon energy cut-off
                      Y_istate(count) = g2_istate(count)*occup_fac*delta
                   end if

                   !Save initial and final electron states
                   istate1(count) = mux_state(el%numbands, m, ik)
                   istate2(count) = mux_state(el%numbands, n, ikp_window)
                end if
             end do !n
          end do !m
       end do !ik

       if(key == 'g') then
          !Change to data output directory
          call chdir(trim(adjustl(num%g2dir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate
          filename = 'gq2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) count
          write(1) g2_istate(1:count)
          close(1)
       end if

       if(key == 'Y') then
          !Multiply constant factor, unit factor, etc.
          Y_istate(1:count) = const*Y_istate(1:count) !THz

          !Change to data output directory
          call chdir(trim(adjustl(Ydir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate
          filename = 'Y.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) count
          write(1) Y_istate(1:count)
          write(1) istate1(1:count)
          write(1) istate2(1:count)
          close(1)
       end if

       !Change back to working directory
       call chdir(num%cwd)
       
       if(key == 'Y') deallocate(g2_istate, Y_istate, istate1, istate2)
    end do
    sync all

    if(key == 'g') then
       !Delete the gReq disk data
       if(this_image() == 1) then
          call chdir(trim(adjustl(num%g2dir)))
          call system('rm gReq.*')
          call chdir(num%cwd)
       end if
    end if
    sync all
  end subroutine calculate_eph_interaction_ibzq
  
  subroutine calculate_eph_interaction_ibzk(wann, crys, el, ph, num, key)
    !! Parallel driver of g2(k,q) over IBZ electron states.
    !!
    !! This subroutine will calculate the full Bloch rep. matrix elements for
    !! all the energy window restricted electron-phonon processes for a given
    !! irreducible initial electron state = (band, wave vector). 
    !! This list will be written to disk in files tagged with the muxed state index.
    !!
    !! key = 'g', 'X' for vertex, transition probability calculation, respectively.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.

    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    character(len = 1), intent(in) :: key
    
    !Local variables
    integer(k8) :: nstates_irred, istate, m, ik, n, ikp, s, &
         iq, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, nprocs, num_active_images
    real(dp) :: k(3), kp(3), q(3), ph_ens_iq(1, ph%numbranches), const, bosefac, &
         fermi_minus_fac, fermi_plus_fac, en_ph, en_el, delta, occup_fac
    real(dp), allocatable :: g2_istate(:), Xplus_istate(:), Xminus_istate(:)
    integer(k8), allocatable :: istate_el(:), istate_ph(:)
    complex(dp) :: gkRp_ik(wann%numwannbands,wann%numwannbands,wann%numbranches,wann%nwsq), &
         ph_evecs_iq(1, ph%numbranches,ph%numbranches)
    character(len = 1024) :: filename, filename_X, Xdir, tag
    logical :: needfinephon

    if(key /= 'g' .and. key /= 'X') then
       call exit_with_message(&
            "Invalid value of key in call to calculate_eph_interaction_ibzk. Exiting.")
    end if

    if(key == 'g') then
       call print_message("Calculating g(k,Rp) -> |g(k,q)|^2 for all IBZ electrons...")
    else
       call print_message("Calculating e-ph transition probabilities for all IBZ electrons...")
    end if

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Xdir = trim(adjustl(num%datadumpdir))//'X_T'//trim(adjustl(tag))
    if(this_image() == 1 .and. key == 'X') then
       !Create directory
       call system('mkdir -p '//trim(adjustl(Xdir)))
    end if
    sync all
    
    !Allocate and initialize g2_istate
    if(key == 'g') then
       !Length of g2_istate
       nprocs = el%nstates_inwindow*wann%numbranches
       allocate(g2_istate(nprocs))
       g2_istate(:) = 0.0_dp
    end if

    !Conversion factor in transition probability expression
    const = twopi/hbar_eVps

    !Total number of IBZ blocks states
    nstates_irred = el%nk_irred*wann%numwannbands
    
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, wann%numwannbands, m, ik)

       if(key == 'g') then
          !Load gkRp(ik) here for use inside the loops below
          call chdir(trim(adjustl(num%g2dir)))
          write (filename, '(I6)') ik
          filename = 'gkRp.ik'//trim(adjustl(filename))
          open(1,file=filename,status="old",access='stream')
          read(1) gkRp_ik
          close(1)
          call chdir(num%cwd)
       end if

       !Get the muxed index of FBZ wave vector from the IBZ blocks index list
       !ik_fbz = el%indexlist_irred(ik)

       !Electron energy
       en_el = el%ens_irred(ik, m)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(en_el - el%enref) > el%fsthick) cycle

       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el%wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*el%kmesh)

       !Load g2_istate from disk for scattering rates calculation
       if(key == 'X') then
          !Change to data output directory
          call chdir(trim(adjustl(num%g2dir)))

          !Read data in binary format
          write (filename, '(I9)') istate
          filename = 'gk2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')
          read(1) nprocs
          if(allocated(g2_istate)) deallocate(g2_istate, Xplus_istate, Xminus_istate, &
               istate_el, istate_ph)
          allocate(g2_istate(nprocs))
          read(1) g2_istate
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
          
          !Allocate and initialize quantities related to transition probabilities
          allocate(Xplus_istate(nprocs), Xminus_istate(nprocs))
          allocate(istate_el(nprocs), istate_ph(nprocs))
          istate_el(:) = 0_k8
          istate_ph(:) = 0_k8
          Xplus_istate(:) = 0.0_dp
          Xminus_istate(:) = 0.0_dp
       end if

       !Initialize eligible process counter for this state
       count = 0
       
       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nk
          !Final wave vector (crystal coords.)
          kp = el%wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*el%kmesh)

          !Find interacting phonon wave vector
          !Note that q, k, and k' are all on the same mesh
          q_indvec = kp_indvec - k_indvec
          needfinephon = .false.
          if(any(mod(q_indvec(:), el%mesh_ref) /= 0_k8)) then
             needfinephon = .true.
             q_indvec = modulo(q_indvec, el%kmesh) !0-based index vector
             q = q_indvec/dble(el%kmesh) !crystal coords.
             !Muxed index of q
             iq = mux_vector(q_indvec, el%kmesh, 0_k8)

             !Calculate the fine mesh phonon.
             call wann%ph_wann_epw(crys, 1_k8, q, ph_ens_iq, ph_evecs_iq)
             !TODO: for consistency, above we should also use ifcs
          else !Original (coarser) mesh phonon
             q_indvec = modulo(q_indvec/el%mesh_ref, ph%qmesh) !0-based index vector
             q = q_indvec/dble(ph%qmesh) !crystal coords.
             !Muxed index of q
             iq = mux_vector(q_indvec, ph%qmesh, 0_k8)
          end if
          
          !Run over final electron bands
          do n = 1, wann%numwannbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle
             
             !Run over phonon branches
             do s = 1, wann%numbranches
                !Increment g2 processes counter
                count = count + 1

                if(key == 'g') then
                   !Calculate |g_mns(<k>,q)|^2
                   if(needfinephon) then
                      g2_istate(count) = wann%g2_epw(crys, k, q, el%evecs_irred(ik, m, :), &
                           el%evecs(ikp, n, :), ph_evecs_iq(1, s, :), &
                           ph_ens_iq(1, s), gkRp_ik, 'ph')
                   else
                      g2_istate(count) = wann%g2_epw(crys, k, q, el%evecs_irred(ik, m, :), &
                           el%evecs(ikp, n, :), ph%evecs(iq, s, :), &
                           ph%ens(iq, s), gkRp_ik, 'ph')
                   end if
                end if

                if(key == 'X') then
                   !Phonon energy
                   if(needfinephon) then
                      en_ph = ph_ens_iq(1, s)
                   else
                      en_ph = ph%ens(iq, s)
                   end if

                   !Bose and Fermi factors
                   if(en_ph /= 0.0_dp) then
                      bosefac = Bose(en_ph, crys%T)
                   else
                      bosefac = 0.0_dp
                   end if
                   fermi_plus_fac = Fermi(en_el + en_ph, el%chempot, crys%T)
                   fermi_minus_fac = Fermi(en_el - en_ph, el%chempot, crys%T)
                   
                   !Calculate X+:

                   !Evaulate delta function
                   delta = delta_fn_tetra(en_el + en_ph, ikp, n, el%kmesh, el%tetramap, &
                        el%tetracount, el%tetra_evals)

                   !Temperature dependent occupation factor
                   occup_fac = bosefac + fermi_plus_fac
                   
                   !Save X+
                   if(en_ph >= 0.5e-3) then !Use a small phonon energy cut-off
                      Xplus_istate(count) = g2_istate(count)*occup_fac*delta
                   end if
                   
                   !Calculate X-:

                   !Evaulate delta function
                   delta = delta_fn_tetra(en_el - en_ph, ikp, n, el%kmesh, el%tetramap, &
                        el%tetracount, el%tetra_evals)

                   !Temperature dependent occupation factor
                   occup_fac = 1.0_dp + bosefac - fermi_minus_fac
                   
                   !Save X-
                   if(en_ph >= 0.5e-3) then !Use a small phonon energy cut-off
                      Xminus_istate(count) = g2_istate(count)*occup_fac*delta
                   end if

                   !Save final electron and interacting phonon states (same for + and -)
                   istate_el(count) = mux_state(el%numbands, n, ikp)
                   if(needfinephon) then
                      !Write fine phonon index as negative so that the iterator
                      !knows to interpolate phonon quantities at this wave vector.
                      istate_ph(count) = -mux_state(ph%numbranches, s, iq)
                   else
                      istate_ph(count) = mux_state(ph%numbranches, s, iq)
                   end if
                end if
             end do !s
          end do !n
       end do !ikp

       if(key == 'g') then
          !Change to data output directory
          call chdir(trim(adjustl(num%g2dir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate
          filename = 'gk2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) count
          write(1) g2_istate
          close(1)
       end if

       if(key == 'X') then
          !Multiply constant factor, unit factor, etc.
          Xplus_istate(1:count) = const*Xplus_istate(1:count) !THz
          Xminus_istate(1:count) = const*Xminus_istate(1:count) !THz

          !Change to data output directory
          call chdir(trim(adjustl(Xdir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate
          filename = 'Xplus.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) count
          write(1) Xplus_istate(1:count)
          write(1) istate_el(1:count)
          write(1) istate_ph(1:count)
          close(1)

          write (filename, '(I9)') istate
          filename = 'Xminus.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'replace', access = 'stream')
          write(1) count
          write(1) Xminus_istate(1:count)
          write(1) istate_el(1:count)
          write(1) istate_ph(1:count)
          close(1)
       end if
       
       !Change back to working directory
       call chdir(num%cwd)

       if(key == 'X') deallocate(g2_istate, Xplus_istate, Xminus_istate, &
            istate_el, istate_ph)
    end do
    sync all
    
    if(key == 'g') then
       !Delete the gkRp disk data
       if(this_image() == 1) then
          call chdir(trim(adjustl(num%g2dir)))
          call system('rm gkRp.*')
          call chdir(num%cwd)
       endif
    end if
    sync all
  end subroutine calculate_eph_interaction_ibzk
  
  subroutine calculate_ph_rta_rates(rta_rates_3ph, rta_rates_phe, num, crys, ph, el)
    !! Subroutine for parallel reading of the 3-ph and ph-e transition probabilities
    !! from disk and calculating the relaxation time approximation (RTA)
    !! scattering rates for the 3-ph and ph-e channels.

    real(dp), allocatable, intent(out) :: rta_rates_3ph(:,:)
    real(dp), allocatable, intent(out) :: rta_rates_phe(:,:)
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(phonon), intent(in) :: ph
    type(electron), intent(in), optional :: el

    !Local variables
    integer(k8) :: nstates_irred, procs, istate, nprocs_3ph, nprocs_phe, &
         iproc, chunk, s, iq, im, num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp), allocatable :: rta_rates_3ph_psum(:,:)[:], rta_rates_phe_psum(:,:)[:], &
         W(:), Y(:)
    character(len = 1024) :: filepath_Wm, filepath_Wp, filepath_Y, Wdir, Ydir, tag
    
    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Wdir = trim(adjustl(num%datadumpdir))//'W_T'//trim(adjustl(tag))
    if(present(el)) then
       Ydir = trim(adjustl(num%datadumpdir))//'Y_T'//trim(adjustl(tag))
    end if
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    !Total number of 3-phonon processes for a given phonon state
    nprocs_3ph = ph%nq*ph%numbranches**2

    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Allocate and initialize scattering rates coarrays
    allocate(rta_rates_3ph_psum(ph%nq_irred, ph%numbranches)[*])
    rta_rates_3ph_psum(:, :) = 0.0_dp
    if(present(el)) then
       allocate(rta_rates_phe_psum(ph%nq_irred, ph%numbranches)[*])
       rta_rates_phe_psum(:, :) = 0.0_dp
    end if
    
    !Allocate and initialize scattering rates
    allocate(rta_rates_3ph(ph%nq_irred, ph%numbranches), rta_rates_phe(ph%nq_irred, ph%numbranches))
    rta_rates_3ph(:, :) = 0.0_dp
    rta_rates_phe(:, :) = 0.0_dp
    
    !Allocate transition probabilities variable
    allocate(W(nprocs_3ph))
    
    !Run over first phonon IBZ states
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbranches, s, iq)

       !Set W+ filename
       write(tag, '(I9)') istate
       filepath_Wp = trim(adjustl(Wdir))//'/Wp.istate'//trim(adjustl(tag))
       
       !Read W+ from file
       call read_transition_probs_3ph(trim(adjustl(filepath_Wp)), W)

       do iproc = 1, nprocs_3ph
          rta_rates_3ph_psum(iq, s) = rta_rates_3ph_psum(iq, s) + W(iproc) 
       end do

       !Set W- filename
       filepath_Wm = trim(adjustl(Wdir))//'/Wm.istate'//trim(adjustl(tag))
       
       !Read W- from file
       call read_transition_probs_3ph(trim(adjustl(filepath_Wm)), W)

       do iproc = 1, nprocs_3ph
          rta_rates_3ph_psum(iq, s) = rta_rates_3ph_psum(iq, s) + 0.5_dp*W(iproc)
       end do

       if(present(el)) then
          !Set Y filename
          filepath_Y = trim(adjustl(Ydir))//'/Y.istate'//trim(adjustl(tag))

          !Read Y from file
          if(allocated(Y)) deallocate(Y)
          call read_transition_probs_eph(trim(adjustl(filepath_Y)), nprocs_phe, Y)

          do iproc = 1, nprocs_phe
             rta_rates_phe_psum(iq, s) = rta_rates_phe_psum(iq, s) + el%spindeg*Y(iproc)
          end do
       end if
    end do

    sync all

    !Reduce coarray partial sums
    do im = 1, num_active_images
       rta_rates_3ph(:,:) = rta_rates_3ph(:,:) + rta_rates_3ph_psum(:,:)[im]
       if(present(el)) then
          rta_rates_phe(:,:) = rta_rates_phe(:,:) + rta_rates_phe_psum(:,:)[im]
       end if
    end do
    sync all
  end subroutine calculate_ph_rta_rates

  subroutine calculate_el_rta_rates(rta_rates_eph, num, crys, el)
    !! Subroutine for parallel reading of the e-ph transition probabilities
    !! from disk and calculating the relaxation time approximation (RTA)
    !! scattering rates for the e-ph channel.

    real(dp), allocatable, intent(out) :: rta_rates_eph(:,:)
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    
    !Local variables
    integer(k8) :: nstates_irred, procs, istate, nprocs_eph, &
         iproc, chunk, m, ik, im, num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp), allocatable :: rta_rates_eph_psum(:,:)[:], X(:)
    character(len = 1024) :: filepath_Xp, filepath_Xm, Xdir, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    Xdir = trim(adjustl(num%datadumpdir))//'X_T'//trim(adjustl(tag))

    !Total number of IBZ blocks states
    nstates_irred = el%nk_irred*el%numbands
    
    !Allocate start and end coarrays
    allocate(start[*], end[*])

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Allocate and initialize scattering rates coarrays
    allocate(rta_rates_eph_psum(el%nk_irred, el%numbands)[*])
    rta_rates_eph_psum(:, :) = 0.0_dp

    !Allocate and initialize scattering rates
    allocate(rta_rates_eph(el%nk_irred, el%numbands))
    rta_rates_eph(:, :) = 0.0_dp

    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, el%numbands, m, ik)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle
       
       !Set X+ filename
       write(tag, '(I9)') istate
       filepath_Xp = trim(adjustl(Xdir))//'/Xplus.istate'//trim(adjustl(tag))
       
       !Read X+ from file
       if(allocated(X)) deallocate(X)
       call read_transition_probs_eph(trim(adjustl(filepath_Xp)), nprocs_eph, X)
       
       do iproc = 1, nprocs_eph
          rta_rates_eph_psum(ik, m) = rta_rates_eph_psum(ik, m) + X(iproc) 
       end do

       !Set X- filename
       write(tag, '(I9)') istate
       filepath_Xm = trim(adjustl(Xdir))//'/Xminus.istate'//trim(adjustl(tag))

       !Read X- from file
       if(allocated(X)) deallocate(X)
       call read_transition_probs_eph(trim(adjustl(filepath_Xm)), nprocs_eph, X)
       
       do iproc = 1, nprocs_eph
          rta_rates_eph_psum(ik, m) = rta_rates_eph_psum(ik, m) + X(iproc) 
       end do
    end do

    sync all

    !Reduce coarray partial sums
    do im = 1, num_active_images
       rta_rates_eph(:,:) = rta_rates_eph(:,:) + rta_rates_eph_psum(:,:)[im]
    end do
    sync all
  end subroutine calculate_el_rta_rates
  
  subroutine read_transition_probs_3ph(filepath, T, istate2, istate3)
    !! Subroutine to read the phonon transition probabilities from disk.

    character(len = *), intent(in) :: filepath
    real(dp), intent(out) :: T(:)
    integer(k8), intent(out), optional :: istate2(:), istate3(:)

    !Read data
    open(1, file = trim(adjustl(filepath)), status = 'old', access = 'stream')
    read(1) T
    if(present(istate2) .and. present(istate3)) then
       read(1) istate2
       read(1) istate3
    end if
    close(1)
  end subroutine read_transition_probs_3ph

  subroutine read_transition_probs_eph(filepath, N, TP, istate1, istate2)
    !! Subroutine to read transition probabilities from disk for e-ph/ph-e processes.
    !! This is different from the read_transition_probs_3ph subroutine in that this
    !! returns the number interaction processes and also allocates the return arrays
    !! accordingly.

    character(len = *), intent(in) :: filepath
    integer(k8), intent(out) :: N
    real(dp), allocatable, intent(out) :: TP(:)
    integer(k8), allocatable, intent(out), optional :: istate1(:), istate2(:)

    !Read data
    open(1, file = trim(adjustl(filepath)), status = 'old', access = 'stream')
    read(1) N
    allocate(TP(N))
    read(1) TP
    if(present(istate1) .and. present(istate2)) then
       allocate(istate1(N), istate2(N))
       read(1) istate1
       read(1) istate2
    end if
    close(1)
  end subroutine read_transition_probs_eph
end module interactions

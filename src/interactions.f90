! Copyright 2020 elphbolt contributors.
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

  use params, only: i64, r64, pi, twopi, amu, qe, hbar_eVps, perm0
  use misc, only: exit_with_message, print_message, distribute_points, &
       demux_state, mux_vector, mux_state, expi, Bose, binsearch, Fermi, &
       twonorm, write2file_rank2_real, demux_vector, interpolate
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics
  use delta, only: delta_fn_tetra, delta_fn_triang

  implicit none

  private
  public calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_ph_rta_rates, read_transition_probs_e, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk, &
       calculate_echimp_interaction_ibzk, calculate_el_rta_rates, &
       calculate_bound_scatt_rates, calculate_thinfilm_scatt_rates, &
       calculate_4ph_rta_rates

  !external chdir, system
  
contains

  pure real(r64) function transfac(v1, v2)
    !! Calculate the "transport factor" that suppresses forward scattering
    !! v1, v2: vectors in cartesian coordinates
    
    real(r64), intent(in) :: v1(3),v2(3)
    real(r64) :: v1sc, v2sc, thresh

    thresh = 1.0e-8_r64
    transfac = 0.0_r64
    v1sc = twonorm(v1)
    v2sc = twonorm(v2)
    if(v1sc /= v2sc .and. v1sc > thresh .and. v2sc > thresh) then
       transfac = 1.0_r64 - dot_product(v1,v2)/v1sc/v2sc
    end if
  end function transfac

  pure real(r64) function qdist(q, reclattvecs)
    !! Function to calculate the smallest wave vector distance in the BZ.
    !! q is in crystal coordinates.
    !! qdist will be in nm^-1
    
    real(r64), intent(in) :: q(3), reclattvecs(3, 3)
    real(r64) :: distfromcorners(3**3)
    integer(i64) :: i, j, k, count

    count = 1
    do i = -1, 1
       do j = -1, 1
          do k = -1, 1
             distfromcorners(count) = twonorm(matmul(reclattvecs, q - (/i, j, k/)))
             count = count + 1
          end do
       end do
    end do
    qdist = minval(distfromcorners)
  end function qdist
  
  pure real(r64) function gchimp2(el, crys, q)
    !! Function to calculate the squared electron-charged impurity vertex.
    !!
    !! This is the Fourier transform of the Yukawa potential, c.f. Eq. 33
    !! of RevModPhys.53.745 (1981).

    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    real(r64), intent(in) :: q

    gchimp2 = 1.0e-3_r64/crys%volume/((perm0*crys%epsilon0)*(q**2 + crys%qTF**2))**2*&
         (el%chimp_conc_n*(qe*el%Zn**2)**2 + el%chimp_conc_p*(qe*el%Zp**2)**2) !ev^2
  end function gchimp2

  pure real(r64) function Vm2_3ph(ev1_s1, conjg_ev2_s2, conjg_ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.
    
    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), conjg_ev2_s2(nb), conjg_ev3_s3(nb)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    complex(r64) :: aux1, aux2, aux3, V0
    
    aux1 = (0.0_r64, 0.0_r64)
    do it = 1, ntrip
       aind = 3*(Index_k(it) - 1)
       bind = 3*(Index_j(it) - 1)
       cind = 3*(Index_i(it) - 1)
       V0 = (0.0_r64, 0.0_r64)
       do a = 1, 3
          aux2 = conjg_ev3_s3(a + aind)
          do b = 1, 3
             aux3 = aux2*conjg_ev2_s2(b + bind)
             do c = 1, 3
                if(ifc3(c, b, a, it) /= 0.0_r64) then
                   V0 = V0 + ifc3(c, b, a, it)*ev1_s1(c + cind)*aux3
                end if
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
    integer(i64) :: start, end, chunk, istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), index_minus, index_plus, &
         neg_iq2, neg_q2_indvec(3), num_active_images, plus_count, minus_count
    real(r64) :: en1, en2, en3, massfac, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         occup_fac, const, bose2, bose3, delta_minus, delta_plus
    real(r64), allocatable :: Vm2_1(:), Vm2_2(:), Wm(:), Wp(:)
    integer(i64), allocatable :: istate2_plus(:), istate3_plus(:), istate2_minus(:), istate3_minus(:)
    complex(r64) :: phases_q2q3(ph%numtriplets)
    character(len = 1024) :: filename, filename_Wm, filename_Wp

    if(key /= 'V' .and. key /= 'W') then
       call exit_with_message("Invalid value of key in call to calculate_3ph_interaction. Exiting.")
    end if

    if(key == 'V') then
       call print_message("Calculating 3-ph vertices for all IBZ phonons...")
    else
       call print_message("Calculating 3-ph transition probabilities for all IBZ phonons...")
    end if
   
    !Conversion factor in transition probability expression
    const = pi/4.0_r64*hbar_eVps**5*(qe/amu)**3*1.0d-12

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    !Maximum total number of 3-phonon processes for a given initial phonon state
    nprocs = ph%nwv*ph%numbands**2
    
    !Allocate |V^-|^2
    if(key == 'V') allocate(Vm2_1(nprocs), Vm2_2(nprocs))
    ! Above, we split the |V-|^2 vertices into two parts:
    ! 1. that are non-zero when the minus-type processes are energetically allowed
    ! 2. that are non-zero when the symmetry-related plus-type processes are energetically allowed

    !Allocate W- and W+
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
       !Load |V^-|^2 from disk for scattering rates calculation
       if(key == 'W') then
          !Change to data output directory
          call chdir(trim(adjustl(num%Vdir)))

          !Read data in binary format
          write (filename, '(I9)') istate1
          filename = 'Vm2.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')

          read(1) minus_count
          if(allocated(Vm2_1)) deallocate(Vm2_1)
          allocate(Vm2_1(minus_count))
          if(minus_count > 0) read(1) Vm2_1

          read(1) plus_count
          if(allocated(Vm2_2)) deallocate(Vm2_2)
          allocate(Vm2_2(plus_count))
          if(plus_count > 0) read(1) Vm2_2
          close(1)

          !Change back to working directory
          call chdir(num%cwd)

          !Initialize transition probabilities
          Wp(:) = 0.0_r64
          Wm(:) = 0.0_r64
          istate2_plus(:) = 0_i64
          istate3_plus(:) = 0_i64
          istate2_minus(:) = 0_i64
          istate3_minus(:) = 0_i64
       end if

       !Initialize transition probabilities
       plus_count = 0_i64
       minus_count = 0_i64
       
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate1, ph%numbands, s1, iq1_ibz)

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Energy of phonon 1
       en1 = ph%ens(iq1, s1)
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)
       
       !Run over second (FBZ) phonon wave vectors
       do iq2 = 1, ph%nwv
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%wvmesh)

          !Folded final phonon wave vector
          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.

          !Muxed index of q3_minus
          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)
          
          if(key == 'V') then
             if(en1 /= 0.0_r64) then
                !Calculate the numtriplet number of mass-normalized phases for this (q2,q3) pair
                do it = 1, ph%numtriplets
                   massfac = 1.0_r64/sqrt(&
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
          do s2 = 1, ph%numbands
             !Energy of phonon 2
             en2 = ph%ens(iq2, s2)

             !Get index of -q2
             neg_q2_indvec = modulo(-q2_indvec, ph%wvmesh)
             neg_iq2 = mux_vector(neg_q2_indvec, ph%wvmesh, 0_i64)

             if(key == 'W') then
                !Bose factor for phonon 2
                bose2 = Bose(en2, crys%T)
             end if
             
             !Run over branches of third phonon
             do s3 = 1, ph%numbands                
                !Minus process index
                index_minus = ((iq2 - 1)*ph%numbands + (s2 - 1))*ph%numbands + s3

                !Energy of phonon 3
                en3 = ph%ens(iq3_minus, s3)

                !Evaluate delta functions
                if(num%tetrahedra) then
                   delta_minus = delta_fn_tetra(en1 - en3, iq2, s2, ph%wvmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals) !minus process

                   delta_plus = delta_fn_tetra(en3 - en1, neg_iq2, s2, ph%wvmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals) !plus process
                else
                   delta_minus = delta_fn_triang(en1 - en3, iq2, s2, ph%wvmesh, ph%triangmap, &
                        ph%triangcount, ph%triang_evals) !minus process

                   delta_plus = delta_fn_triang(en3 - en1, neg_iq2, s2, ph%wvmesh, ph%triangmap, &
                        ph%triangcount, ph%triang_evals) !plus process
                end if
                
                if(key == 'V') then
                   if(en1*en2*en3 == 0.0_r64) cycle

                   if(delta_minus > 0.0_r64) then
                      !Increase counter for energetically available minus process
                      minus_count = minus_count + 1

                      !Save the index of this process
                      !V1_indexlist(minus_count) = index_minus

                      !Calculate and save the minus process vertex
                      Vm2_1(minus_count) = Vm2_3ph(ph%evecs(iq1, s1, :), &
                           conjg(ph%evecs(iq2, s2, :)), conjg(ph%evecs(iq3_minus, s3, :)), &
                           ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
                           phases_q2q3, ph%numtriplets, ph%numbands)
                   end if

                   if(delta_plus > 0.0_r64) then
                      !Increase counter for energetically available plus process
                      plus_count = plus_count + 1

                      !Calculate and save the minus process vertex
                      Vm2_2(plus_count) = Vm2_3ph(ph%evecs(iq1, s1, :), &
                           conjg(ph%evecs(iq2, s2, :)), conjg(ph%evecs(iq3_minus, s3, :)), &
                           ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
                           phases_q2q3, ph%numtriplets, ph%numbands)
                   end if
                end if

                if(key == 'W') then
                   if(en1*en2*en3 == 0.0_r64) cycle
                   
                   !Bose factor for phonon 3
                   bose3 = Bose(en3, crys%T)

                   !Calculate W-:
                   
                   !Temperature dependent occupation factor
                   !(bose1 + 1)*bose2*bose3/(bose1*(bose1 + 1))
                   ! = (bose2 + bose3 + 1)
                   occup_fac = (bose2 + bose3 + 1.0_r64)

                   if(delta_minus > 0.0_r64) then
                      !Non-zero process counter
                      minus_count = minus_count + 1

                      !Save W-
                      !Wm(minus_count) = Vm2(index_minus)*occup_fac*delta_minus/en1/en2/en3
                      Wm(minus_count) = Vm2_1(minus_count)*occup_fac*delta_minus/en1/en2/en3
                      istate2_minus(minus_count) = mux_state(ph%numbands, s2, iq2)
                      istate3_minus(minus_count) = mux_state(ph%numbands, s3, iq3_minus)
                   end if
                 
                   !Calculate W+:

                   !Grab index of corresponding plus process using
                   !V-(s1q1|s2q2,s3q3) = V+(s1q1|s2-q2,s3q3)
                   index_plus = ((neg_iq2 - 1)*ph%numbands + (s2 - 1))*ph%numbands + s3

                   !Temperature dependent occupation factor
                   !(bose1 + 1)*(bose2 + 1)*bose3/(bose1*(bose1 + 1))
                   ! = bose2 - bose3.
                   occup_fac = (bose2 - bose3)

                   if(delta_plus > 0.0_r64) then
                      !Non-zero process counter
                      plus_count = plus_count + 1

                      !Save W+
                      Wp(plus_count) = Vm2_2(plus_count)*occup_fac*delta_plus/en1/en2/en3
                      istate2_plus(plus_count) = mux_state(ph%numbands, s2, neg_iq2)
                      istate3_plus(plus_count) = mux_state(ph%numbands, s3, iq3_minus)
                   end if
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
          write(1) minus_count
          write(1) Vm2_1(1:minus_count)
          write(1) plus_count
          write(1) Vm2_2(1:plus_count)
          close(1)
       end if

       if(key == 'W') then
          !Multiply constant factor, unit factor, etc.
          Wm(:) = const*Wm(:) !THz
          Wp(:) = const*Wp(:) !THz

          !Write W+ and W- to disk
          !Change to data output directory
          call chdir(trim(adjustl(num%Wdir)))

          !Write data in binary format
          !Note: this will overwrite existing data!
          write (filename, '(I9)') istate1

          filename_Wm = 'Wm.istate'//trim(adjustl(filename))
          open(1, file = trim(filename_Wm), status = 'replace', access = 'stream')
          write(1) minus_count
          write(1) Wm(1:minus_count)
          write(1) istate2_minus(1:minus_count)
          write(1) istate3_minus(1:minus_count)
          close(1)

          filename_Wp = 'Wp.istate'//trim(adjustl(filename))
          open(1, file = trim(filename_Wp), status = 'replace', access = 'stream')
          write(1) plus_count
          write(1) Wp(1:plus_count)
          write(1) istate2_plus(1:plus_count)
          write(1) istate3_plus(1:plus_count)
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
    integer(i64) :: iq, iqstart, iqend, chunk, num_active_images

    call print_message("Calculating g(Re,Rp) -> g(Re,q) for all IBZ q...")

    call distribute_points(ph%nwv_irred, chunk, iqstart, iqend, num_active_images)

    if(this_image() == 1) then
       print*, "   #q = ", ph%nwv_irred
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
    integer(i64) :: ik, ikstart, ikend, chunk, num_active_images

    call print_message("Calculating g(Re,Rp) -> g(k,Rp) for all IBZ k...")

    call distribute_points(el%nwv_irred, chunk, ikstart, ikend, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #k = ", el%nwv_irred
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
    integer(i64) :: nstates_irred, istate, m, iq, iq_fbz, n, ik, ikp, s, &
         ikp_window, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), nprocs, count, num_active_images
    integer(i64), allocatable :: istate1(:), istate2(:)
    real(r64) :: k(3), q(3), en_ph, en_el, en_elp, const, delta, &
         invboseplus1, fermi1, fermi2, occup_fac
    real(r64), allocatable :: g2_istate(:), Y_istate(:)
    complex(r64), allocatable :: gReq_iq(:,:,:,:)
    character(len = 1024) :: filename
    
    if(key /= 'g' .and. key /= 'Y') then
       call exit_with_message(&
            "Invalid value of key in call to calculate_eph_interaction_ibzq. Exiting.")
    end if

    if(key == 'g') then
       call print_message("Calculating g(Re,q) -> |g(k,q)|^2 for all IBZ phonons...")
    else
       call print_message("Calculating ph-e transition probabilities for all IBZ phonons...")
    end if
    
    !Allocate and initialize gReq_iq and g2_istate
    if(key == 'g') then
       allocate(gReq_iq(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsk))
       
       !Maximum length of g2_istate
       nprocs = el%nstates_inwindow*ph%numbands
       allocate(g2_istate(nprocs))
       g2_istate(:) = 0.0_r64
    end if

    !Conversion factor in transition probability expression
    const = twopi/hbar_eVps
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if
    
    do istate = start, end !over IBZ blocks states
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)

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
          if(en_ph /= 0.0_r64) then
             invboseplus1 = 1.0_r64/(1.0_r64 + Bose(en_ph, crys%T))
          else
             invboseplus1 = 0.0_r64
          end if
       end if
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       q = ph%wavevecs(iq_fbz, :)

       !Convert from crystal to 0-based index vector
       q_indvec = nint(q*ph%wvmesh)

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
          if(nprocs > 0) read(1) g2_istate
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
          
          !Allocate and initialize quantities related to transition probabilities
          allocate(Y_istate(nprocs))
          allocate(istate1(nprocs), istate2(nprocs))
          istate1(:) = -1_i64
          istate2(:) = -1_i64
          Y_istate(:) = 0.0_r64
       end if

       !Initialize process counter
       count = 0
       
       !Run over initial (in-window, FBZ blocks) electron wave vectors
       do ik = 1, el%nwv
          !Initial wave vector (crystal coords.)
          k = el%wavevecs(ik, :)

          !Convert from crystal to 0-based index vector
          k_indvec = nint(k*el%wvmesh)
          
          !Find final electron wave vector
          kp_indvec = modulo(k_indvec + el%mesh_ref_array*q_indvec, el%wvmesh) !0-based index vector

          !Muxed index of kp
          ikp = mux_vector(kp_indvec, el%wvmesh, 0_i64)

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
                   if(num%tetrahedra) then
                      delta = delta_fn_tetra(en_elp - en_ph, ik, m, el%wvmesh, el%tetramap, &
                           el%tetracount, el%tetra_evals)
                   else
                      delta = delta_fn_triang(en_elp - en_ph, ik, m, el%wvmesh, el%triangmap, &
                           el%triangcount, el%triang_evals)
                   end if

                   !Temperature dependent occupation factor
                   occup_fac = fermi1*(1.0_r64 - fermi2)*invboseplus1
                   
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
          call chdir(trim(adjustl(num%Ydir)))

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
    integer(i64) :: nstates_irred, istate, m, ik, n, ikp, s, &
         iq, start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, nprocs, num_active_images
    real(r64) :: k(3), kp(3), q(3), ph_ens_iq(1, ph%numbands), qlist(1, 3), &
         const, bosefac, fermi_minus_fac, fermi_plus_fac, en_ph, en_el, delta, occup_fac
    real(r64), allocatable :: g2_istate(:), Xplus_istate(:), Xminus_istate(:)
    integer(i64), allocatable :: istate_el(:), istate_ph(:)
    complex(r64), allocatable :: gkRp_ik(:, :, :, :)
    complex(r64) :: ph_evecs_iq(1, ph%numbands,ph%numbands)
    character(len = 1024) :: filename
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
    
    !Allocate and initialize gkRp_ik and g2_istate
    if(key == 'g') then
       allocate(gkRp_ik(wann%numwannbands,wann%numwannbands,wann%numbranches,wann%nwsq))
       
       !Length of g2_istate
       nprocs = el%nstates_inwindow*wann%numbranches
       allocate(g2_istate(nprocs))
       g2_istate(:) = 0.0_r64
    end if

    !Conversion factor in transition probability expression
    const = twopi/hbar_eVps

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

       !Electron energy
       en_el = el%ens_irred(ik, m)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(en_el - el%enref) > el%fsthick) cycle

       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el%wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*el%wvmesh)

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
          if(nprocs > 0) read(1) g2_istate
          close(1)

          !Change back to working directory
          call chdir(num%cwd)
          
          !Allocate and initialize quantities related to transition probabilities
          allocate(Xplus_istate(nprocs), Xminus_istate(nprocs))
          allocate(istate_el(nprocs), istate_ph(nprocs))
          istate_el(:) = 0_i64
          istate_ph(:) = 0_i64
          Xplus_istate(:) = 0.0_r64
          Xminus_istate(:) = 0.0_r64
       end if

       !Initialize eligible process counter for this state
       count = 0
       
       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nwv
          !Final wave vector (crystal coords.)
          kp = el%wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*el%wvmesh)

          !Find interacting phonon wave vector
          !Note that q, k, and k' are all on the same mesh
          q_indvec = kp_indvec - k_indvec
          needfinephon = .false.
          if(any(mod(q_indvec(:), el%mesh_ref_array) /= 0_i64)) then
             needfinephon = .true.
             q_indvec = modulo(q_indvec, el%wvmesh) !0-based index vector
             q = q_indvec/dble(el%wvmesh) !crystal coords.
             !Muxed index of q
             iq = mux_vector(q_indvec, el%wvmesh, 0_i64)

             !Calculate the fine mesh phonon.
             qlist(1, :) = q
             call wann%ph_wann_epw(crys, 1_i64, qlist, ph_ens_iq, ph_evecs_iq)
          else !Original (coarser) mesh phonon
             q_indvec = modulo(q_indvec/el%mesh_ref_array, ph%wvmesh) !0-based index vector
             q = q_indvec/dble(ph%wvmesh) !crystal coords.
             !Muxed index of q
             iq = mux_vector(q_indvec, ph%wvmesh, 0_i64)
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
                   if(en_ph /= 0.0_r64) then
                      bosefac = Bose(en_ph, crys%T)
                   else
                      bosefac = 0.0_r64
                   end if
                   fermi_plus_fac = Fermi(en_el + en_ph, el%chempot, crys%T)
                   fermi_minus_fac = Fermi(en_el - en_ph, el%chempot, crys%T)
                   
                   !Calculate X+:

                   !Evaulate delta function
                   if(num%tetrahedra) then
                      delta = delta_fn_tetra(en_el + en_ph, ikp, n, el%wvmesh, el%tetramap, &
                           el%tetracount, el%tetra_evals)
                   else
                      delta = delta_fn_triang(en_el + en_ph, ikp, n, el%wvmesh, el%triangmap, &
                           el%triangcount, el%triang_evals)
                   end if

                   !Temperature dependent occupation factor
                   occup_fac = bosefac + fermi_plus_fac
                   
                   !Save X+
                   if(en_ph >= 0.5e-3) then !Use a small phonon energy cut-off
                      Xplus_istate(count) = g2_istate(count)*occup_fac*delta
                   end if
                   
                   !Calculate X-:

                   !Evaulate delta function
                   if(num%tetrahedra) then
                      delta = delta_fn_tetra(en_el - en_ph, ikp, n, el%wvmesh, el%tetramap, &
                           el%tetracount, el%tetra_evals)
                   else
                      delta = delta_fn_triang(en_el - en_ph, ikp, n, el%wvmesh, el%triangmap, &
                           el%triangcount, el%triang_evals)
                   end if

                   !Temperature dependent occupation factor
                   occup_fac = 1.0_r64 + bosefac - fermi_minus_fac
                   
                   !Save X-
                   if(en_ph >= 0.5e-3) then !Use a small phonon energy cut-off
                      Xminus_istate(count) = g2_istate(count)*occup_fac*delta
                   end if

                   !Save final electron and interacting phonon states (same for + and -)
                   istate_el(count) = mux_state(el%numbands, n, ikp)
                   if(needfinephon) then
                      !Write fine phonon index as negative so that the iterator
                      !knows to interpolate phonon quantities at this wave vector.
                      istate_ph(count) = -mux_state(ph%numbands, s, iq)
                   else
                      istate_ph(count) = mux_state(ph%numbands, s, iq)
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
          call chdir(trim(adjustl(num%Xdir)))

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

  subroutine calculate_echimp_interaction_ibzk(crys, el, num)
    !! Parallel driver of |g_e-chimp(k,k')|^2 over IBZ electron states.
    !!
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside the energy
    !window will be skipped in the calculation as they are irrelevant for transport.

    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(i64) :: nstates_irred, istate, m, ik, n, ikp, &
         start, end, chunk, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, nprocs, num_active_images
    real(r64) :: k(3), kp(3), q_mag, const, en_el, delta, g2
    real(r64), allocatable :: Xchimp_istate(:)
    integer(i64), allocatable :: istate_el(:)
    character(len = 1024) :: filename

    call print_message("Calculating e-ch. imp. transition probabilities for all IBZ electrons...")
    
    !Conversion factor in transition probability expression
    const = twopi/hbar_eVps

    !Number of processes
    nprocs = el%nstates_inwindow
    allocate(Xchimp_istate(nprocs), istate_el(nprocs))
    Xchimp_istate(:) = 0.0_r64
    istate_el(:) = 0_i64
    
    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*el%numbands
    
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    if(this_image() == 1) then
       write(*, "(A, I10)") " #states = ", nstates_irred
       write(*, "(A, I10)") " #states/image = ", chunk
    end if

    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, el%numbands, m, ik)

       !Electron energy
       en_el = el%ens_irred(ik, m)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(en_el - el%enref) > el%fsthick) cycle
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el%wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*el%wvmesh)

       !Initialize eligible process counter for this state
       count = 0
       
       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, el%nwv
          !Final wave vector (crystal coords.)
          kp = el%wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*el%wvmesh)

          !Find interacting phonon wave vector
          !Note that q, k, and k' are all on the same mesh
          q_indvec = kp_indvec - k_indvec !0-based index vector

          !Calculate length of the wave vector
          q_mag = qdist(q_indvec/dble(el%wvmesh), crys%reclattvecs)

          !Calculate matrix element
          g2 = gchimp2(el, crys, q_mag)
         
          !Run over final electron bands
          do n = 1, el%numbands
             !Apply energy window to final electron
             if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle
                      
             !Increment g2 processes counter
             count = count + 1
             
             !Evaulate delta function
             if(num%tetrahedra) then
                delta = delta_fn_tetra(en_el, ikp, n, el%wvmesh, el%tetramap, &
                     el%tetracount, el%tetra_evals)
             else
                delta = delta_fn_triang(en_el, ikp, n, el%wvmesh, el%triangmap, &
                     el%triangcount, el%triang_evals)
             end if

             !Save Xchimp
             Xchimp_istate(count) = g2*transfac(matmul(crys%reclattvecs,k), &
                  matmul(crys%reclattvecs,kp))*delta
             
             !Save final electron state
             istate_el(count) = mux_state(el%numbands, n, ikp)
             end do !n
       end do !ikp

       !Multiply constant factor, unit factor, etc.
       Xchimp_istate(1:count) = const*Xchimp_istate(1:count) !THz
       
       !Change to data output directory
       call chdir(trim(adjustl(num%Xdir)))

       !Write data in binary format
       !Note: this will overwrite existing data!
       write (filename, '(I9)') istate
       filename = 'Xchimp.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'replace', access = 'stream')
       write(1) count
       write(1) Xchimp_istate(1:count)
       write(1) istate_el(1:count)
       close(1)       

       !Change back to working directory
       call chdir(num%cwd)
    end do
    sync all
  end subroutine calculate_echimp_interaction_ibzk
  
  subroutine calculate_ph_rta_rates(rta_rates_3ph, rta_rates_phe, num, crys, ph, el)
    !! Subroutine for parallel reading of the 3-ph and ph-e transition probabilities
    !! from disk and calculating the relaxation time approximation (RTA)
    !! scattering rates for the 3-ph and ph-e channels.

    real(r64), allocatable, intent(out) :: rta_rates_3ph(:,:)
    real(r64), allocatable, intent(out) :: rta_rates_phe(:,:)
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(phonon), intent(in) :: ph
    type(electron), intent(in), optional :: el

    !Local variables
    integer(i64) :: nstates_irred, istate, nprocs_3ph_plus, nprocs_3ph_minus, &
         nprocs_phe, iproc, chunk, s, iq, num_active_images, start, end
    real(r64), allocatable :: W(:), Y(:)
    character(len = 1024) :: filepath_Wm, filepath_Wp, filepath_Y, tag
    
    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands
    
    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    !Allocate and initialize scattering rates
    allocate(rta_rates_3ph(ph%nwv_irred, ph%numbands), rta_rates_phe(ph%nwv_irred, ph%numbands))
    rta_rates_3ph(:, :) = 0.0_r64
    rta_rates_phe(:, :) = 0.0_r64
        
    !Run over first phonon IBZ states
    do istate = start, end
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate, ph%numbands, s, iq)

       !Set W+ filename
       write(tag, '(I9)') istate
       filepath_Wp = trim(adjustl(num%Wdir))//'/Wp.istate'//trim(adjustl(tag))
       
       !Read W+ from file
       if(allocated(W)) deallocate(W)
       call read_transition_probs_e(trim(adjustl(filepath_Wp)), nprocs_3ph_plus, W)

       do iproc = 1, nprocs_3ph_plus
          rta_rates_3ph(iq, s) = rta_rates_3ph(iq, s) + W(iproc) 
       end do

       !Set W- filename
       filepath_Wm = trim(adjustl(num%Wdir))//'/Wm.istate'//trim(adjustl(tag))
       
       !Read W- from file
       if(allocated(W)) deallocate(W)
       call read_transition_probs_e(trim(adjustl(filepath_Wm)), nprocs_3ph_minus, W)

       do iproc = 1, nprocs_3ph_minus
          rta_rates_3ph(iq, s) = rta_rates_3ph(iq, s) + 0.5_r64*W(iproc)
       end do

       if(present(el)) then
          !Set Y filename
          filepath_Y = trim(adjustl(num%Ydir))//'/Y.istate'//trim(adjustl(tag))

          !Read Y from file
          if(allocated(Y)) deallocate(Y)
          call read_transition_probs_e(trim(adjustl(filepath_Y)), nprocs_phe, Y)

          do iproc = 1, nprocs_phe
             rta_rates_phe(iq, s) = rta_rates_phe(iq, s) + el%spindeg*Y(iproc)
          end do
       end if
    end do

    !Reduce partial sums
    sync all
    call co_sum(rta_rates_3ph)
    sync all
    if(present(el)) then
       sync all
       call co_sum(rta_rates_phe)
       sync all
    end if
  end subroutine calculate_ph_rta_rates

  subroutine calculate_4ph_rta_rates(rta_rates, num, crys, ph)
    !! Subroutine for interporlating 4-ph scattering rates from an
    !! external coarser mesh calculation.
    
    real(r64), allocatable, intent(out) :: rta_rates(:,:)
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(phonon), intent(in) :: ph

    !Local variables
    integer(i64) :: chunk, s, iq, coarse_numq_irred, coarse_numq_full, &
         num_active_images, start, end, fineq_indvec(3), mesh_ref_array(3), &
         coarse_qmesh(3), fbz2ibz
    real(r64) :: ignore
    real(r64), allocatable :: coarse_rta_rates_ibz(:, :), coarse_rta_rates_fbz(:, :)
    character(len=1024) :: temp_tag, filename

    allocate(rta_rates(ph%nwv_irred, ph%numbands))
    rta_rates = 0.0_r64

    if(num%fourph) then
       !Set some internal mesh related variables
       mesh_ref_array = num%fourph_mesh_ref
       if(crys%twod) mesh_ref_array(3) = 1
       coarse_qmesh = ph%wvmesh/num%fourph_mesh_ref
       if(crys%twod) coarse_qmesh(3) = 1
       
       write(temp_tag, "(E9.3)") crys%T
       filename = "FourPhonon_BTE.w_4ph_T" // trim(adjustl(temp_tag))

       !Read number of irreducible q-points
       if(this_image() == 1 .and. num%fourph) then
          open(1, file=filename, status="old")
          read(1, *) coarse_numq_irred
       end if

       sync all
       call co_broadcast(coarse_numq_irred, 1)
       sync all
       
       allocate(coarse_rta_rates_ibz(coarse_numq_irred, ph%numbands))
       
       !Read coarse mesh, IBZ 4-ph scattering rates
       if(this_image() == 1 .and. num%fourph) then
          open(1, file=filename, status="old")

          do s = 1, ph%numbands
             do iq = 1, coarse_numq_irred
                read(1, *) ignore, coarse_rta_rates_ibz(iq, s)
             end do
          end do

          close(1)
       end if
       sync all
       call co_broadcast(coarse_rta_rates_ibz, 1)
       sync all

       !Compute FBZ <-> IBZ mapping from BTE.qpoints_full file
       if(this_image() == 1 .and. num%fourph) then
          open(1, file="FourPhonon_BTE.qpoints_full", status="old")

          read(1, *) coarse_numq_full
          if(coarse_numq_full /= product(coarse_qmesh)) then
             call exit_with_message('Wrong q-mesh in external 4-ph calculation. Exiting.')
          end if
       end if

       sync all
       call co_broadcast(coarse_numq_full, 1)
       sync all
          
       allocate(coarse_rta_rates_fbz(coarse_numq_full, ph%numbands))

       if(this_image() == 1 .and. num%fourph) then
          do iq = 1, coarse_numq_full
             read(1, *) ignore, fbz2ibz, ignore, ignore, ignore
             coarse_rta_rates_fbz(iq, :) = coarse_rta_rates_ibz(fbz2ibz, :)
          end do

          close(1)
       end if
       sync all
       call co_broadcast(coarse_rta_rates_fbz, 1)
       sync all
       
       !Divide phonon states among images
       call distribute_points(ph%nwv_irred, chunk, start, end, num_active_images)

       !If needed, parallely interpolate over fine q-mesh
       if(num%fourph_mesh_ref > 1) then
          do iq = start, end
             !Calculate the fine mesh wave vector, 0-based index vector
             call demux_vector(ph%indexlist_irred(iq), fineq_indvec, ph%wvmesh, 0_i64)

             !Interpolate 4-ph scattering rates on this wave vector
             do s = 1, ph%numbands
                call interpolate(coarse_qmesh, mesh_ref_array, coarse_rta_rates_fbz(:, s), &
                     fineq_indvec, rta_rates(iq, s))
             end do
          end do

          sync all
          call co_sum(rta_rates)
          sync all
       else
          rta_rates = coarse_rta_rates_ibz
       end if
    end if
  end subroutine calculate_4ph_rta_rates
  
  subroutine calculate_el_rta_rates(rta_rates_eph, rta_rates_echimp, num, crys, el)
    !! Subroutine for parallel reading of the e-ph transition probabilities
    !! from disk and calculating the relaxation time approximation (RTA)
    !! scattering rates for the e-ph channel.

    real(r64), allocatable, intent(out) :: rta_rates_eph(:,:), rta_rates_echimp(:,:)
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(electron), intent(in) :: el
    
    !Local variables
    integer(i64) :: nstates_irred, istate, nprocs_eph, &
         iproc, chunk, m, ik, num_active_images, start, end
    real(r64), allocatable :: X(:)
    character(len = 1024) :: filepath_Xp, filepath_Xm, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T

    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*el%numbands

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    !Allocate and initialize scattering rates
    allocate(rta_rates_eph(el%nwv_irred, el%numbands))
    rta_rates_eph(:, :) = 0.0_r64
    allocate(rta_rates_echimp(el%nwv_irred, el%numbands))
    rta_rates_echimp(:, :) = 0.0_r64

    do istate = start, end !over IBZ blocks states
       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate, el%numbands, m, ik)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle
       
       !Set X+ filename
       write(tag, '(I9)') istate
       filepath_Xp = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))
       
       !Read X+ from file
       if(allocated(X)) deallocate(X)
       call read_transition_probs_e(trim(adjustl(filepath_Xp)), nprocs_eph, X)
       
       do iproc = 1, nprocs_eph
          rta_rates_eph(ik, m) = rta_rates_eph(ik, m) + X(iproc) 
       end do

       !Set X- filename
       write(tag, '(I9)') istate
       filepath_Xm = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

       !Read X- from file
       if(allocated(X)) deallocate(X)
       call read_transition_probs_e(trim(adjustl(filepath_Xm)), nprocs_eph, X)

       do iproc = 1, nprocs_eph
          rta_rates_eph(ik, m) = rta_rates_eph(ik, m) + X(iproc) 
       end do

       if(num%elchimp) then
          !Set Xchimp filename
          write(tag, '(I9)') istate
          filepath_Xm = trim(adjustl(num%Xdir))//'/Xchimp.istate'//trim(adjustl(tag))

          !Read Xchimp from file
          if(allocated(X)) deallocate(X)
          call read_transition_probs_e(trim(adjustl(filepath_Xm)), nprocs_eph, X)

          do iproc = 1, nprocs_eph
             rta_rates_echimp(ik, m) = rta_rates_echimp(ik, m) + X(iproc) 
          end do
       end if
    end do

    !Reduce partial sums
    sync all
    call co_sum(rta_rates_eph)
    sync all
    if(num%elchimp) then
       sync all
       call co_sum(rta_rates_echimp)
       sync all
    end if    
  end subroutine calculate_el_rta_rates
  
  subroutine read_transition_probs_e(filepath, N, TP, istate1, istate2)
    !! Subroutine to read transition probabilities from disk for interaction processes.

    character(len = *), intent(in) :: filepath
    integer(i64), intent(out) :: N
    real(r64), allocatable, intent(out) :: TP(:)
    integer(i64), allocatable, intent(out), optional :: istate1(:), istate2(:)

    !Read data
    open(1, file = trim(adjustl(filepath)), status = 'old', access = 'stream')
    read(1) N
    allocate(TP(N))
    if(N > 0) read(1) TP
    if(present(istate1)) then
       allocate(istate1(N))
       if(N > 0) then
          read(1) istate1
       end if
       if(present(istate2)) then
          allocate(istate2(N))
          if(N > 0) then
             read(1) istate2
          end if
       end if
    end if
    close(1)
  end subroutine read_transition_probs_e

  subroutine calculate_bound_scatt_rates(prefix, finite_crys, length, vels_fbz, &
       indexlist_irred, scatt_rates)
    !! Subroutine to calculate the phonon/electron-boundary scattering rates.
    !!
    !! prefix Type of particle
    !! finite_crys Is the crystal finite?
    !! length Characteristic boundary scattering length scale in mm
    !! vels Velocities on the FBZ
    !! indexlist_irred List of muxed indices of the IBZ wedge.
    !! scatt_rates Boundary scattering rates on the IBZ

    character(len = 2), intent(in) :: prefix
    logical, intent(in) :: finite_crys
    real(r64), intent(in) :: length
    real(r64), intent(in) :: vels_fbz(:,:,:)
    integer(i64), intent(in) :: indexlist_irred(:)
    real(r64), allocatable, intent(out) :: scatt_rates(:,:)

    !Local variables
    integer(i64) :: ik, ib, nk_irred, nb

    !Number of IBZ wave vectors and bands
    nk_irred = size(indexlist_irred(:))
    nb = size(vels_fbz(1,:,1))
    
    !Allocate boundary scattering rates and initialize to infinite crystal values
    allocate(scatt_rates(nk_irred, nb))
    scatt_rates = 0.0_r64
    
    !Check finiteness of crystal
    if(finite_crys) then
       do ik = 1, nk_irred
          do ib = 1, nb
             scatt_rates(ik, ib) = twonorm(vels_fbz(indexlist_irred(ik), ib, :))&
                  /length*1.e-6_r64 !THz
          end do
       end do
    end if

    !Write to file
    call write2file_rank2_real(prefix // '.W_rta_'//prefix//'bound', scatt_rates)
  end subroutine calculate_bound_scatt_rates

  subroutine calculate_thinfilm_scatt_rates(prefix, finite_crys, height, normal, vels_fbz, &
       indexlist_irred, scatt_rates)
    !! Subroutine to calculate the phonon/electron-thin-film scattering rates.
    !!
    !! prefix Type of particle
    !! finite_crys Is the crystal finite?
    !! height Height of thin-film in mm
    !! normal Normal direction to thin-film
    !! vels Velocities on the FBZ
    !! indexlist_irred List of muxed indices of the IBZ wedge.
    !! scatt_rates Thin-film scattering rates on the IBZ

    character(len = 2), intent(in) :: prefix
    logical, intent(in) :: finite_crys
    real(r64), intent(in) :: height
    character(1), intent(in) :: normal
    real(r64), intent(in) :: vels_fbz(:,:,:)
    integer(i64), intent(in) :: indexlist_irred(:)
    real(r64), allocatable, intent(out) :: scatt_rates(:,:)

    !Local variables
    integer(i64) :: ik, ib, nk_irred, nb, dir

    !Number of IBZ wave vectors and bands
    nk_irred = size(indexlist_irred(:))
    nb = size(vels_fbz(1,:,1))

    !Allocate boundary scattering rates and initialize to infinite crystal values
    allocate(scatt_rates(nk_irred, nb))
    scatt_rates = 0.0_r64

    if(normal == 'x') then
       dir = 1_i64
    else if(normal == 'y') then
       dir = 2_i64
    else if(normal == 'z') then
       dir = 3_i64
    else
       call exit_with_message("Bad thin-film normal direction in calculate_thinfilm_scattrates. Exiting.")
    end if
    
    !Check finiteness of crystal
    if(finite_crys) then
       do ik = 1, nk_irred
          do ib = 1, nb
             scatt_rates(ik, ib) = abs(vels_fbz(indexlist_irred(ik), ib, dir)) &
                  /height*1.e-6_r64 !THz
          end do
       end do
    end if
    scatt_rates = 2.0_r64*scatt_rates

    !Write to file
    call write2file_rank2_real(prefix // '.W_rta_'//prefix//'thinfilm', scatt_rates)
  end subroutine calculate_thinfilm_scatt_rates

!!$  subroutine calculate_defect_scatt_rates(prefix, def_frac, indexlist_ibz, ens_fbz, diagT)!, scatt_rates)
!!$    !! Subroutine to calculate the phonon-defect scattering rate given
!!$    !! the diagonal of the scattering T-matrix.
!!$    !!
!!$    !! prefix Particle type label
!!$    !! def_frac Elemental fraction of defects
!!$    !! ens IBZ energies
!!$    !! diagT Diagonal of the IBZ T-matrix
!!$    !! scatt_rates IBZ Scattering rates
!!$
!!$    character(len = 2), intent(in) :: prefix
!!$    real(r64), intent(in) :: def_frac
!!$    real(r64), intent(in) :: ens_fbz(:, :)
!!$    integer(i64), intent(in) :: indexlist_ibz(:)
!!$    complex(r64), intent(in) :: diagT(:, :)
!!$    !real(r64), allocatable, intent(out) :: scatt_rates(:, :)
!!$    
!!$    !Local variables
!!$    integer(i64) :: nk_ibz, nbands, ik
!!$    real(r64), allocatable :: scatt_rates(:, :)
!!$
!!$    nk_ibz = size(diagT, 1)
!!$    nbands = size(diagT, 2)
!!$
!!$    print*, 'def_frac = ', def_frac
!!$    
!!$    allocate(scatt_rates(nk_ibz, nbands))
!!$
!!$    do ik = 1, nk_ibz
!!$       scatt_rates(ik, :) = imag(diagT(ik, :))/ens_fbz(indexlist_ibz(ik), :)
!!$    end do
!!$
!!$    scatt_rates = -def_frac*scatt_rates/hbar_eVps
!!$
!!$    !Deal with Gamma point acoustic phonons! and zero-velocity optic phonons
!!$    scatt_rates(1, 1:3) = 0.0_r64
!!$    
!!$    !Write to file
!!$    call write2file_rank2_real(prefix // '.W_rta_'//prefix//'defect', scatt_rates)
!!$  end subroutine calculate_defect_scatt_rates
end module interactions

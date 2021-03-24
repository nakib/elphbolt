module interactions
  !! Module containing the procedures related to the computation of interactions.

  use params, only: k4, dp
  use misc, only: print_message, distribute_points, demux_state, mux_vector, expi
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal
  use electron_module, only: electron
  use phonon_module, only: phonon
  use numerics_module, only: numerics

  implicit none

  private
  public calculate_g_mixed, calculate_g2_bloch, calculate_V_3ph

contains
  
  subroutine calculate_V_phases(ph, crys, phases)
    !! Parallel calculator of the mass-normalized phase factor
    !! needed for the 3-phonon vertex.
    !!
    !! Warning: precalculation of phases requires a massive amount of ram.
    !! TODO: Think about whether the phases should be saved to disk instead.
    
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    complex(dp), intent(out) :: phases(ph%numtriplets, ph%nq, ph%nq)

    !Local variables
    integer(k4) :: chunk, it, counter, im, iq2, iq3
    integer(k4), allocatable :: start[:], end[:]
    real(dp) :: q2_cart(3), q3_cart(3), massfac
    complex(dp), allocatable :: phases_chunk(:,:,:)[:]

    !Allocate the start and end coarrays
    allocate(start[*], end[*])
    
    !Divide number of triplets among images
    call distribute_points(ph%numtriplets, chunk, start, end)
    
    !Allocate small work variable chunk for each image
    allocate(phases_chunk(end-start+1, ph%nq, ph%nq)[*])

    !Calculate phases
    counter = 0
    do it = start, end
       counter = counter + 1
       massfac = 1.0_dp/sqrt(&
            crys%masses(crys%atomtypes(ph%Index_i(it)))*&
            crys%masses(crys%atomtypes(ph%Index_j(it)))*&
            crys%masses(crys%atomtypes(ph%Index_k(it))))
       do iq2 = 1, ph%nq
          q2_cart = matmul(crys%reclattvecs, ph%wavevecs(iq2,:))
          do iq3 = 1, ph%nq
             q3_cart = matmul(crys%reclattvecs, ph%wavevecs(iq3,:))
             phases_chunk(counter, iq2, iq3) = massfac*&
                  expi(-dot_product(q2_cart, ph%R_j(:,it)))*&
                  expi(-dot_product(q3_cart, ph%R_k(:,it)))
          end do
       end do
    end do
    sync all
    
    !Collect chunks into globally accessible phases
    do im = 1, num_images()
       phases(start[im]:end[im],:,:) = phases_chunk(:,:,:)[im]
    end do
    sync all        
  end subroutine calculate_V_phases

!!$  subroutine test_phase_calculator(ph, crys)
!!$    type(phonon), intent(in) :: ph
!!$    type(crystal), intent(in) :: crys
!!$    complex(dp) :: phases(ph%numtriplets, ph%nq, ph%nq)
!!$
!!$    phases(:,:,:) = (-1.0_dp,1.0_dp)
!!$    
!!$    call calculate_V_phases(ph, crys, phases)
!!$
!!$    if(this_image() == 1) then
!!$       print*, '....'
!!$       print*, phases(1:5,ph%nq,ph%nq)
!!$       print*, phases(1:5,1,ph%nq)
!!$       print*, phases(1:5,ph%nq,1)
!!$       print*, phases(ph%numtriplets,ph%nq,ph%nq)
!!$    end if
!!$  end subroutine test_phase_calculator
  
!!$  function Vp_minus(i,j,k,q,qprime,qdprime,realqprime,realqdprime,eigenvect,&
!!$       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
!!$    implicit none
!!$
!!$    integer(kind=4),intent(in) :: i
!!$    integer(kind=4),intent(in) :: j
!!$    integer(kind=4),intent(in) :: k
!!$    integer(kind=4),intent(in) :: q
!!$    integer(kind=4),intent(in) :: qprime
!!$    integer(kind=4),intent(in) :: qdprime
!!$    real(kind=8),intent(in) :: realqprime(3)
!!$    real(kind=8),intent(in) :: realqdprime(3)    
!!$    complex(kind=8),intent(in) :: eigenvect(nptq,Nbranches,Nbranches)
!!$    integer(kind=4),intent(in) :: Ntri
!!$    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
!!$    real(kind=8),intent(in) :: R_j(3,Ntri)
!!$    real(kind=8),intent(in) :: R_k(3,Ntri)
!!$    integer(kind=4),intent(in) :: Index_i(Ntri)
!!$    integer(kind=4),intent(in) :: Index_j(Ntri)
!!$    integer(kind=4),intent(in) :: Index_k(Ntri)
!!$
!!$    !complex(kind=8) :: Vp_minus
!!$    complex(kind=8) :: aux
!!$    real(kind=8) :: Vp_minus
!!$
!!$    integer(kind=4) :: ll
!!$    integer(kind=4) :: rr
!!$    integer(kind=4) :: ss
!!$    integer(kind=4) :: tt
!!$    complex(kind=8) :: prefactor
!!$    complex(kind=8) :: Vp0
!!$
!!$    Vp_minus=0.d0
!!$
!!$    do ll=1,Ntri
!!$       prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
!!$            masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
!!$            phexp(-dot_product(realqprime,R_j(:,ll)))*&
!!$            phexp(-dot_product(realqdprime,R_k(:,ll)))
!!$       Vp0=0.
!!$       do rr=1,3
!!$          do ss=1,3
!!$             do tt=1,3
!!$                Vp0=Vp0+Phi(tt,ss,rr,ll)*&
!!$                     eigenvect(q,i,tt+3*(Index_i(ll)-1))*&
!!$                     conjg(eigenvect(qprime,j,ss+3*(Index_j(ll)-1)))*&
!!$                     conjg(eigenvect(qdprime,k,rr+3*(Index_k(ll)-1)))
!!$             end do
!!$          end do
!!$       end do
!!$       !Vp_minus=Vp_minus+prefactor*Vp0
!!$       aux=aux+prefactor*Vp0
!!$    end do
!!$    Vp_minus=abs(aux)**2
!!$  end function Vp_minus
!!$  
  function Vm2_3ph(ph, s1, s2, s3, iq1, iq2, iq3, phases_q2q3)
    !! Function to calculate the 3-ph interaction vertex |V-|^2.
    !! TODO This is a bottleneck. Think about how to speed this function up.
    
    type(phonon), intent(in) :: ph
    integer(k4), intent(in) :: s1, s2, s3, iq1, iq2, iq3
    complex(dp), intent(in) :: phases_q2q3(ph%numtriplets)
    real(dp) :: Vm2_3ph

    !Local variables
    integer(k4) :: it, a, b, c
    complex(dp) :: aux, V0

    aux = (0.0_dp, 0.0_dp)
    do it = 1, ph%numtriplets
       V0 = (0.0_dp, 0.0_dp)
       do a = 1, 3
          do b = 1, 3
             do c = 1, 3
                V0 = V0 + ph%ifc3(c, b, a, it)*&
                     ph%evecs(iq1, s1, c + 3*(ph%Index_i(it) - 1))*&
                     conjg(ph%evecs(iq2, s2, b + 3*(ph%Index_j(it) - 1)))*&
                     conjg(ph%evecs(iq3, s3, a + 3*(ph%Index_k(it) - 1)))
             end do
          end do
       end do
       aux = aux + V0*phases_q2q3(it)
    end do
    Vm2_3ph = abs(aux)**2
  end function Vm2_3ph
  
  subroutine calculate_V_3ph(ph, crys, num)
    !! Parallel driver of the 3-ph vertex calculator for all IBZ phonon wave vectors.
    !! This subroutine calculates |V-(<q1>|q2,q3)|^2 for each irreducible phonon
    !! and saves the result to disk.

    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(k4) :: start, end, chunk, istate1, nstates_irred, nstates, &
         s1, s2, s3, iq1_ibz, iq1, iq2, iq3, count, it, &
         q1_indvec(3), q2_indvec(3), q3_indvec(3)
    real(dp) :: en1, en2, en3, massfac, q1(3), q2(3), q3(3), q2_cart(3), q3_cart(3)
    real(dp), allocatable :: Vm2(:)
    complex(dp) :: phases_q2q3(ph%numtriplets)
    character(len = 1024) :: filename
    
    call print_message("Calculating 3-ph vertices for all IBZ q...")

!!$    !Calculate mass-normalized phases needed for the V calculation
!!$    call calculate_V_phases(ph)
    
    !Total number of IBZ blocks states
    nstates_irred = ph%nq_irred*ph%numbranches

    !Total number of processes for a given initial phonon state
    nstates = ph%nq*ph%numbranches**2

    allocate(Vm2(nstates))

    !Divide phonon states among images
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
       call demux_state(istate1, ph%numbranches, s1, iq1_ibz)

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)
       
       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%qmesh)

       !Energy of phonon 1
       en1 = ph%ens(iq1, s1)
       
       !Run over second (FBZ) phonon wave vectors
       do iq2 = 1, ph%nq
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%qmesh)

          !Folded final phonon wave vector
          q3_indvec = modulo(q1_indvec - q2_indvec, ph%qmesh) !0-based index vector
          q3 = q3_indvec/dble(ph%qmesh) !crystal coords.

          !Muxed index of q3
          iq3 = mux_vector(q3_indvec, ph%qmesh, 0_dp)
          
          !Calculate the numtriplet number of mass-normalized phases for this (q2,q3) pair
          do it = 1, ph%numtriplets
             massfac = 1.0_dp/sqrt(&
                  crys%masses(crys%atomtypes(ph%Index_i(it)))*&
                  crys%masses(crys%atomtypes(ph%Index_j(it)))*&
                  crys%masses(crys%atomtypes(ph%Index_k(it))))
             q2_cart = matmul(crys%reclattvecs, ph%wavevecs(iq2,:))
             q3_cart = matmul(crys%reclattvecs, ph%wavevecs(iq3,:))
             phases_q2q3(it) = massfac*&
                  expi(-dot_product(q2_cart, ph%R_j(:,it)))*&
                  expi(-dot_product(q3_cart, ph%R_k(:,it)))
          end do
          
          !Run over branches of second phonon
          do s2 = 1, ph%numbranches
             !Energy of phonon 2
             en2 = ph%ens(iq2, s2)

             !Run over branches of third phonon
             do s3 = 1, ph%numbranches
                !Energy of phonon 3
                en3 = ph%ens(iq3, s3)
                
                !Count process
                count = count + 1

                !Minus process index
                !index_minus = (iq2 - 1)*nstates + (s2 - 1)*ph%numbranches + s3
                
                !Calculate the minus process vertex
                Vm2(count) = Vm2_3ph(ph, s1, s2, s3, iq1, iq2, iq3, phases_q2q3)

!!$                !Calculate corresponding plus process using
!!$                !V-(s1q1|s2q2,s3q3) = V+(s1q1|s2-q2,s3q3)
!!$                neg_q2_indvec = modulo(-q2_indvec, ph%qmesh)
!!$                neg_iq2 = mux_vector(neg_q2_indvec, ph%qmesh, 0_dp)
!!$                index_plus = (neg_iq2 - 1)*nstates + (s2 - 1)*ph%numbranches + s3
!!$                Vp2(index_plus) = Vm2(index_minus)                
             end do !s3
          end do !s2
       end do !iq2

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
    sync all
  end subroutine calculate_g2_bloch
end module interactions

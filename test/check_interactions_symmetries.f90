program check_interactions_symmetries

  use precision, only: i64, r64
  use electron_module, onlY: electron
  use phonon_module, only: phonon
  use wannier_module, only: Wannier
  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use symmetry_module, only: symmetry
  use interactions, only: Vm2_3ph
  use misc, only: mux_vector, expi, binsearch
  use testify_m, only: testify

  implicit none

  ! TODO set up minimal scaffolding to test interactions_on_star
  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(phonon) :: ph
  type(electron) :: el
  type(wannier) :: wann

  integer(i64) :: iq1_ibz, iq2, ib1, ib2, ib3, &
       ik1_ibz, ik2, ieq, m, n, s
  real(r64), allocatable :: V2(:, :, :, :), g2(:, :, :, :)
  real(r64) :: aux

  !Set up crystal
  call crys%initialize
  
  !Set up numerics data
  call num%initialize(crys)
  
  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)
  
  !Calculate phonons
  call ph%initialize(crys, sym, num)
  
  !A trivial case.
  !iq1_ibz = 2
  !iq2 = 1

!!$  !A non-trivial case.
!!$  iq1_ibz = 9
!!$  iq2 = 5
!!$  
!!$  allocate(V2(ph%numbands, ph%numbands, ph%numbands, &
!!$         ph%nequiv(iq1_ibz)))
!!$  V2 = V2_on_star(iq1_ibz, iq2, ph, crys)
!!$
!!$  do ieq = 1, ph%nequiv(iq1_ibz)
!!$     do ib3 = 1, ph%numbands
!!$        do ib2 = 1, ph%numbands
!!$           do ib1 = 1, ph%numbands
!!$              !These won't be equal as this quantity is not gauge invariant
!!$              print*, ib1, ib2, ib3, V2(ib1, ib2, ib3, ieq)
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  print*, '...'
!!$
!!$  do ib1 = 1, ph%numbands
!!$     print*, 'branch number:', ib1
!!$     print*, '              image#   \sum_{s2 s3}|V^2(s1q1,s2q2|s3q3)|'
!!$     do ieq = 1, ph%nequiv(iq1_ibz)
!!$        aux = 0.0_r64
!!$        do ib3 = 1, ph%numbands
!!$           do ib2 = 1, ph%numbands
!!$              aux = aux + V2(ib1, ib2, ib3, ieq)
!!$           end do
!!$        end do
!!$        print*, ieq, aux
!!$     end do
!!$  end do

  !Test g2
  !Read EPW Wannier data
  call wann%read(num)

  !Calculate electrons
  call el%initialize(wann, crys, sym, num)

  ik1_ibz = 4
  ik2 = 10
  
  allocate(g2(wann%numwannbands, wann%numwannbands, ph%numbands, &
         el%nequiv(ik1_ibz)))
  g2 = g2_on_star(ik1_ibz, ik2, el, ph, wann, crys, num)

!!$  do ieq = 1, el%nequiv(ik1_ibz)
!!$     do m = 1, 1!wann%numwannbands
!!$        do n = 1, 1!wann%numwannbands
!!$           do s = 1, 1!ph%numbands
!!$              !These won't be equal as this quantity is not gauge invariant
!!$              print*, m, n, s, g2(m, n, s, ieq)
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  print*, '...'
  
  do m = 1, wann%numwannbands
     print*, 'band number:', m
     print*, '              image#  \sum_{s n}|g(mk,nk+q|sq)|^2'
     do ieq = 1, el%nequiv(ik1_ibz)
        aux = 0.0_r64
        do n = 1, wann%numwannbands
           do s = 1, ph%numbands
              aux = aux + g2(m, n, s, ieq)
           end do
        end do
        print*, ieq, aux
     end do
  end do
  
contains

  function V2_on_star(iwv1_ibz, iwv2, ph, crys) result(intstar)
    !! Directly computes the ph-ph vertex on the star
    !! of the given wave vector.
    !!
    !! iwv1_ibz Index of 1st, IBZ wave vector
    !! iwv2_ibz Index of 2nd wave vector
    !! ph Phonon object
    !! crys Crystal object

    integer(i64), intent(in) :: iwv1_ibz, iwv2
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys

    integer(i64) :: ieq, s1, s2, s3, it, isym, iwv1_image, iwv2_image, iwv3, &
         q1_indvec(3), q2_indvec(3), q3_indvec(3)
    real(r64) :: intstar(ph%numbands, ph%numbands, ph%numbands, &
         ph%nequiv(iwv1_ibz)), &
         q1(3), q2(3), q3(3), q2_cart(3), q3_cart(3)
    complex(r64) :: phases(ph%numtriplets)

    print*, 'Initial IBZ wave vector', ph%wavevecs_irred(iwv1_ibz, :)
    print*, 'Number of equivalent points', ph%nequiv(iwv1_ibz)
    do ieq = 1, ph%nequiv(iwv1_ibz)
       !Image of the IBZ (on the FBZ) point under this symmetry operation
       isym = ph%ibz2fbz_map(ieq, iwv1_ibz, 1) !symmetry
       iwv1_image = ph%ibz2fbz_map(ieq, iwv1_ibz, 2) !image due to symmetry

       !Initial wave vector (crystal coords.)
       q1 = ph%wavevecs(iwv1_image, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)
       
       !2nd wave vector (crystal coords.) under this symmetry
       iwv2_image = ph%equiv_map(isym, iwv2)
       q2 = ph%wavevecs(iwv2_image, :)

       !Convert from crystal to 0-based index vector
       q2_indvec = nint(q2*ph%wvmesh)

       !Folded final phonon wave vector
       q3_indvec = modulo(q1_indvec - q2_indvec, &
            ph%wvmesh) !0-based index vector
       q3 = q3_indvec/dble(ph%wvmesh) !crystal coords.
       
       !Muxed index of q3
       iwv3 = mux_vector(q3_indvec, ph%wvmesh, 0_i64)

       q2_cart = matmul(crys%reclattvecs, q2)
       q3_cart = matmul(crys%reclattvecs, q3)
       do it = 1, ph%numtriplets
          phases(it) = &
               expi(-dot_product(q2_cart, (ph%R_j(:, it))) &
                    -dot_product(q3_cart, (ph%R_k(:, it))))
       end do
       
       do s3 = 1, ph%numbands
          do s2 = 1, ph%numbands       
             do s1 = 1, ph%numbands
                intstar(s1, s2, s3, ieq) = Vm2_3ph(ph%evecs(iwv1_image, s1, :), &
                     ph%evecs(iwv2_image, s2, :), ph%evecs(iwv3, s3, :), &
                     ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), &
                     ph%ifc3(:,:,:,:), phases(:), &
                     ph%numtriplets, ph%numbands)
             end do
          end do
       end do
    end do
  end function V2_on_star

  function g2_on_star(iwv1_ibz, iwv2, &
       el, ph, wann, crys, num) result(intstar)
    !! Directly computes the e-ph vertex on the star
    !! of the given wave vector.
    !!
    !! iwv1_ibz Index of 1st, IBZ wave vector
    !! iwv2_ibz Index of 2nd wave vector
    !! el Electron object
    !! ph Phonon object
    !! wann Wannier object
    !! crys Crystal object
    !! num Numerics object

    integer(i64), intent(in) :: iwv1_ibz, iwv2
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num

    integer(i64) :: ieq, s, m, n, isym, iwv1_image, iwv2_image, iq, &
         k1_indvec(3), k2_indvec(3), q_indvec(3)
    real(r64) :: intstar(wann%numwannbands, wann%numwannbands, ph%numbands, &
         el%nequiv(iwv1_ibz)), ph_ens_iq(1, ph%numbands), &
         k1(3), k2(3), q(3), qlist(1, 3), discardible(1, ph%numbands), &
         el_ens_ik1(wann%numwannbands), el_ens_ik2(wann%numwannbands)
    complex(r64) :: gmixed_k(wann%numwannbands, &
         wann%numwannbands, wann%numbranches, wann%nwsq), &
         gkRp_iwv1_image(wann%numwannbands, wann%numwannbands, &
         wann%numbranches, wann%nwsq), &
         ph_evecs_iq(1, ph%numbands,ph%numbands)
    character(1024) :: filetag, filename

    !Conform gwann to the best shape for the contraction in gkRp.
    call wann%reshape_gwann_for_gkRp
    sync all
    
    print*, 'Initial IBZ wave vector', el%wavevecs_irred(iwv1_ibz, :)
    print*, 'Number of equivalent points', el%nequiv(iwv1_ibz)
    do ieq = 1, el%nequiv(iwv1_ibz)
       !Image of the IBZ (on the FBZ) point under this symmetry operation
       isym = el%ibz2fbz_map(ieq, iwv1_ibz, 1) !symmetry
       call binsearch(el%indexlist, el%ibz2fbz_map(ieq, iwv1_ibz, 2), &
            iwv1_image) !image due to symmetry

       !Calculate g(k, Rp) from g(Re, Rp)
       call wann%gkRp(num, iwv1_image, el%wavevecs(iwv1_image, :))

       !Load gkRp from file
       !Change to data output directory
       call chdir(trim(adjustl(num%g2dir)))
       write (filetag, '(I9)') iwv1_image
       filename = 'gkRp.ik'//adjustl(trim(filetag))
       open(1, file = filename, status = "old", access = 'stream')
       read(1) gkRp_iwv1_image
       close(1)
       !Change back to working directory
       call chdir(num%cwd)

       !Initial wave vector (crystal coords.)
       k1 = el%wavevecs(iwv1_image, :)

       !Convert from crystal to 0-based index vector
       k1_indvec = nint(k1*el%wvmesh)
       
       !2nd wave vector (crystal coords.) under this symmetry
       !print*, isym, iwv2, el%equiv_map(isym, iwv2)
       call binsearch(el%indexlist, el%equiv_map(isym, iwv2), &
            iwv2_image)
       k2 = el%wavevecs(iwv2_image, :)

       !Convert from crystal to 0-based index vector
       k2_indvec = nint(k2*el%wvmesh)

       !This is the crudest way to symmetrize the electron energies.
       !Note that the eigenvectors will not be the same.
       if(ieq == 1) then
          el_ens_ik1 = el%ens(iwv1_image, :)
          el_ens_ik2 = el%ens(iwv2_image, :)
       end if
       
       !Folded final phonon wave vector
       q_indvec = modulo(k2_indvec - k1_indvec, el%wvmesh) !0-based index vector
       q = q_indvec/dble(el%wvmesh) !crystal coords.
       
       !Muxed index of q
       iq = mux_vector(q_indvec, el%wvmesh, 0_i64)

       qlist(1, :) = q
       
       !This is the crudest way to symmetrize the phonon energies.
       !Note that the eigenvectors will not be the same.
       if(ieq == 1) then
          call wann%ph_wann(crys, 1_i64, qlist, ph_ens_iq, ph_evecs_iq)
       else
          call wann%ph_wann(crys, 1_i64, qlist, discardible, ph_evecs_iq)
       end if
              
       do s = 1, ph%numbands !band index of phonon
          do n = 1, wann%numwannbands !band index of final electron
             do m = 1, wann%numwannbands !band index of initial electron
                intstar(m, n, s, ieq) = &
                     wann%g2(crys, k1, q, el%evecs(iwv1_image, m, :), &
                     el%evecs(iwv2_image, n, :), ph_evecs_iq(1, s, :), &
                     ph_ens_iq(1, s), gkRp_iwv1_image, 'ph')
             end do
          end do
       end do

       call average_over_degenerate_subspace(intstar(:, :, :, ieq), &
            ph_ens_iq(1, :), el_ens_ik1, el_ens_ik2)
    end do
    
    !Put gwann back to original shape
    call wann%reshape_gwann_for_gkRp(revert = .true.)
    sync all

!!$    !Stress test
!!$    !Multiply with a phonon dependent random function
!!$    call srand(1236)
!!$    do s = 1, ph%numbands
!!$       intstar(:, :, s, :) = rand()*intstar(:, :, s, :)
!!$    end do
  end function g2_on_star

  subroutine average_over_degenerate_subspace(g2, ph_ens_q, el_ens_k1, el_ens_k2)

    real(r64), intent(inout) :: g2(:, :, :)
    real(r64), intent(in) :: ph_ens_q(:), el_ens_k1(:), el_ens_k2(:)

    real(r64) :: numbands, numbranches, aux, ph_en, el_en, thres
    integer :: m, mp, n, np, s, sp, deg_count

    thres = 1.0e-6_r64 !closeness
    
    numbands = size(el_ens_k1)
    numbranches = size(ph_ens_q)

    !Average over degenerate phonon branches
    do m = 1, numbands
       do n = 1, numbands
          do s = 1, numbranches
             deg_count = 0
             aux = 0.0_r64
             ph_en = ph_ens_q(s)
             do sp = 1, numbranches
                if(abs(ph_en - ph_ens_q(sp)) < thres) then
                   deg_count = deg_count + 1
                   aux = aux + g2(m, n, sp)
                end if
             end do
             g2(m, n, s) = aux/dble(deg_count)
          end do
       end do
    end do
    
    !Average over initial electron bands
    do s = 1, numbranches
       do n = 1, numbands
          do m = 1, numbands
             deg_count = 0
             aux = 0.0_r64
             el_en = el_ens_k1(m)
             do mp = 1, numbands
                if(abs(el_en - el_ens_k1(mp)) < thres) then
                   deg_count = deg_count + 1
                   aux = aux + g2(mp, n, s)
                end if
             end do
             g2(m, n, s) = aux/dble(deg_count)
          end do
       end do
    end do

    !Average over final electron bands
    do s = 1, numbranches
       do m = 1, numbands
          do n = 1, numbands
             deg_count = 0
             aux = 0.0_r64
             el_en = el_ens_k2(n)
             do np = 1, numbands
                if(abs(el_en - el_ens_k2(np)) < thres) then
                   deg_count = deg_count + 1
                   aux = aux + g2(m, np, s)
                end if
             end do
             g2(m, n, s) = aux/dble(deg_count)
          end do
       end do
    end do
end subroutine average_over_degenerate_subspace
  
end program check_interactions_symmetries

module screening_module

  use precision, only: r64, i64
  use params, only: kB, qe, pi, perm0, oneI, hbar, hbar_evps, me
  use electron_module, only: electron
  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use misc, only: linspace, mux_vector, binsearch, Fermi, print_message, &
       compsimps, twonorm, write2file_rank2_real, write2file_rank1_real, &
       distribute_points, sort, qdist, operator(.umklapp.), Bose, &
       Hilbert_transform, demux_vector
  use wannier_module, only: wannier
  use delta, only: delta_fn, get_delta_fn_pointer
  use vector_allreps_module, only: vec=>vector_allreps, &
       vec_add=>vector_allreps_add

  implicit none

  private
  public calculate_qTF, &
       spectral_head_polarizability_3d_q, calculate_RPA_dielectric_2d_model!, &
       !calculate_RPA_dielectric_3d_G0_scratch, &
  
contains

  subroutine calculate_qTF(crys, el)
    !! Calculate Thomas-Fermi screening wave vector from the static
    !! limit of the Lindhard function.
    !
    !Captain's log May 7, 2024. I would like to turn this into a pure function.
    !Why do we even need crystal to have qTF as a member?
    !It is only ever accessed by gchimp2.
    !This is a super cheap calculation anyway...

    type(crystal), intent(inout) :: crys
    type(electron), intent(in) :: el

    !Local variables
    real(r64) :: beta, fFD
    integer(i64) :: ib, ik

    beta = 1.0_r64/kB/crys%T/qe !1/J
    crys%qTF = 0.0_r64

    call print_message("Calculating Thomas-Fermi screening wave vector...")

    do ib = 1, el%numbands
       do ik = 1, el%nwv
          fFD = Fermi(el%ens(ik, ib), el%chempot, crys%T)
          crys%qTF = crys%qTF + fFD*(1.0_r64 - fFD)
       end do
    end do

    !Free-electron gas Thomas-Fermi model
    ! qTF**2 = spindeg*e^2*beta/nptq/vol_pcell/perm0*Sum_{BZ}f0_{k}(1-f0_{k})
    crys%qTF = sqrt(1.0e9_r64*crys%qTF*el%spindeg*beta*qe**2/product(el%wvmesh)&
         /crys%volume/perm0) !nm^-1
    if(crys%twod) crys%qTF = crys%thickness*crys%qTF**2/2  ! nm^-1

    if(this_image() == 1) then
       write(*, "(A, 1E16.8, A)") ' Thomas-Fermi screening wave vector = ', crys%qTF, ' 1/nm'
    end if
  end subroutine calculate_qTF

  subroutine spectral_head_polarizability_3d_q(spec_eps, Omegas, qvec, &
       el, crys, tetrahedra)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! qvec Transfer wave vector
    !! el Electron data type
    !! crys Crystal data type
    !! tetrahedra Delta evaulator selector

    real(r64), intent(in) :: Omegas(:)
    type(vec), intent(in) :: qvec
    type(electron), intent(in) :: el
    type(crystal), intent(in) :: crys
    logical, intent(in) :: tetrahedra
    real(r64), allocatable, intent(out) :: spec_eps(:)

    !Locals
    integer(i64) :: m, n, ik, iOmega, nOmegas, k_indvec(3), kp_indvec(3), where_in_indexlist
    real(r64) :: overlap, ek, ekp, el_ens_kp(1, el%numbands), kppathvecs(1, 3)
    complex(r64) :: el_evecs_kp(1, el%numbands, el%numbands)
    procedure(delta_fn), pointer :: delta_fn_ptr => null()
    type(vec) :: kvec, kpvec
    real(r64) :: dOmega

    nOmegas = size(Omegas)

    dOmega = Omegas(2) - Omegas(1)

    allocate(spec_eps(nOmegas))

    !Associate delta function procedure pointer
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra)

    spec_eps = 0.0
    !Below, we will sum over k, m, and n
    do ik = 1, el%nwv
       !This k vector
       kvec = vec(el%indexlist(ik), el%wvmesh, crys%reclattvecs)

       !This k' = k + q vector
       kpvec = vec_add(kvec, qvec, el%wvmesh, crys%reclattvecs)

       !Binary search for k' in in k-vectors (full BZ) index list
       call binsearch(el%indexlist, kpvec%muxed_index, where_in_indexlist)

       !Ignore terms for which k' is outside the Fermi shell
       if(where_in_indexlist < 0) cycle

       el_ens_kp(1, :) = el%ens(where_in_indexlist, :)
       el_evecs_kp(1, :, :) = el%evecs(where_in_indexlist, :, :)

       do m = 1, el%numbands
          ek = el%ens(ik, m)

          !Apply energy window to initial electron
          if(abs(ek - el%enref) > el%fsthick) cycle

          do iOmega = nOmegas/2 + 2, nOmegas !positive energy sector
             do n = 1, el%numbands

                ekp = el_ens_kp(1, n)

                !Apply energy window to final electron
                if(abs(ekp - el%enref) > el%fsthick) cycle

                !This is |U(k')U^\dagger(k)|_nm squared
                !(Recall that U^\dagger(k) is the diagonalizer of the electronic hamiltonian.)
                overlap = (abs(dot_product(el_evecs_kp(1, n, :), el%evecs(ik, m, :))))**2

                spec_eps(iOmega) = spec_eps(iOmega) + &
                     (Fermi(ek, el%chempot, crys%T) - &
                     Fermi(ekp, el%chempot, crys%T))*overlap* &
                     delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
                     el%wvmesh, el%simplex_map, &
                     el%simplex_count, el%simplex_evals)
             end do
          end do
       end do
    end do

!   do iOmega = 1, nOmegas
!      !Recall that the resolvent is already normalized in the full wave vector mesh.
!      !As such, the 1/product(el%wvmesh) is not needed in the expression below.
!      spec_eps(iOmega) = spec_eps(iOmega)*el%spindeg/crys%volume
!   end do
    spec_eps = spec_eps*el%spindeg/crys%volume
    if(crys%twod) spec_eps = spec_eps*crys%thickness
    !At this point [spec_eps] = nm^-3.eV^-1 or nm^-2.eV^-1 for 3D or 2D case

    !The negative energy sector
    do iOmega = 1, nOmegas/2
       spec_eps(iOmega) = -spec_eps(nOmegas + 1 - iOmega)
    end do

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end subroutine spectral_head_polarizability_3d_q

  subroutine spectral_head_polarizability_3d_qpath(spec_eps, Omegas, qvec, &
       el, wann, crys, tetrahedra)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! qvec Transfer wave vector
    !! el Electron data type
    !! wann Wannier data type
    !! crys Crystal data type
    !! tetrahedra Delta evaulator selector

    real(r64), intent(in) :: Omegas(:)
    type(vec), intent(in) :: qvec
    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    logical, intent(in) :: tetrahedra
    real(r64), allocatable, intent(out) :: spec_eps(:)

    !Locals
    integer(i64) :: m, n, ik, iOmega, nOmegas, k_indvec(3), kp_indvec(3), where_in_indexlist
    real(r64) :: overlap, ek, ekp, el_ens_kp(1, el%numbands), kppathvecs(1, 3)
    complex(r64) :: el_evecs_kp(1, el%numbands, el%numbands)
    procedure(delta_fn), pointer :: delta_fn_ptr => null()
    type(vec) :: kvec, kpvec

    real(r64) :: dOmega

    nOmegas = size(Omegas)

    dOmega = Omegas(2) - Omegas(1)

    allocate(spec_eps(nOmegas))

    !Associate delta function procedure pointer
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra)

    spec_eps = 0.0
    !Below, we will sum over k, m, and n
    do ik = 1, el%nwv
       !This k vector
       kvec = vec(el%indexlist(ik), el%wvmesh, crys%reclattvecs)

       !This k' = k + q vector
       kpvec = vec_add(kvec, qvec, el%wvmesh, crys%reclattvecs)

       !Binary search for k' in in k-vectors (full BZ) index list
       !call binsearch(el%indexlist, kpvec%muxed_index, where_in_indexlist)

       if(where_in_indexlist < 0) then !k' does not exist in index list
          !Put this k' vector (in fractional coordinates) in a list
          kppathvecs(1, :) = kpvec%frac

          !Compute electronic properties
          call wann%el_wann(crys = crys, &
               nk = 1_i64, &
               kvecs = kppathvecs, &
               energies = el_ens_kp, &
               evecs = el_evecs_kp, &
               scissor = el%scissor)
       else !k' already exists in index list
          el_ens_kp(1, :) = el%ens(where_in_indexlist, :)
          el_evecs_kp(1, :, :) = el%evecs(where_in_indexlist, :, :)
       end if

       do m = 1, wann%numwannbands
          ek = el%ens(ik, m)

          !Apply energy window to initial electron
          if(abs(ek - el%enref) > el%fsthick) cycle

          do iOmega = nOmegas/2 + 2, nOmegas !positive energy sector
             do n = 1, wann%numwannbands

                ekp = el_ens_kp(1, n)

                !Apply energy window to final electron
                if(abs(ekp - el%enref) > el%fsthick) cycle

                !This is |U(k')U^\dagger(k)|_nm squared
                !(Recall that U^\dagger(k) is the diagonalizer of the electronic hamiltonian.)
                overlap = (abs(dot_product(el_evecs_kp(1, n, :), el%evecs(ik, m, :))))**2

                spec_eps(iOmega) = spec_eps(iOmega) + &
                     (Fermi(ek, el%chempot, crys%T) - &
                     Fermi(ekp, el%chempot, crys%T))*overlap* &
                     delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
                     el%wvmesh, el%simplex_map, &
                     el%simplex_count, el%simplex_evals)
             end do
          end do
       end do
    end do

    do iOmega = 1, nOmegas
       !Recall that the resolvent is already normalized in the full wave vector mesh.
       !As such, the 1/product(el%wvmesh) is not needed in the expression below.
       spec_eps(iOmega) = spec_eps(iOmega)*el%spindeg/crys%volume
    end do
    !At this point [spec_eps] = nm^-3.eV^-1

    !The negative energy sector
    do iOmega = 1, nOmegas/2
       spec_eps(iOmega) = -spec_eps(nOmegas + 1 - iOmega)
    end do

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end subroutine spectral_head_polarizability_3d_qpath

  subroutine spectral_head_polarizability_2d_qpath_old(spec_eps, Omegas, qcrys, &
       el, wann, crys, tetrahedra)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! qvec Wave vector (fractional) of excitation in the electron gas
    !! el Electron data type

    real(r64), intent(in) :: Omegas(:), qcrys(3)
    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys

    real(r64), allocatable, intent(out) :: spec_eps(:)
    logical, intent(in) :: tetrahedra

    !Locals
    integer(i64) :: m, n, ik, iOmega, nOmegas, k_indvec(3), kp_indvec(3), no, npo
    real(r64) :: overlap, ek, ekp, delta, Omega_l, Omega_r, &
         el_ens_kp(1, el%numbands), kppathvecs(1, 3), kvec(3), kpqvec(3),&
         kpvel(1, el%numbands, 3)
    complex(r64) :: el_evecs_kp(1, el%numbands, el%numbands)
    procedure(delta_fn), pointer :: delta_fn_ptr => null()
    real(r64) :: onebyroot2pi, sigma 

    onebyroot2pi = 1.0_r64/sqrt(2.0*pi)

    nOmegas = size(Omegas)

    allocate(spec_eps(nOmegas))

    !Associate delta function procedure pointer ! 2D
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra = .false.)
    !delta_fn_ptr => get_delta_fn_pointer(tetrahedra)

    spec_eps = 0.0
    do ik = 1, el%nwv
       
       kppathvecs(1, :) = el%wavevecs(ik, :).umklapp.qcrys
       
       call wann%el_wann(crys = crys, &
        nk = 1_i64, &
        kvecs = kppathvecs, &
        energies = el_ens_kp, velocities = kpvel, &
        evecs = el_evecs_kp, &
        scissor = el%scissor)
!$!        call wann%el_wann(crys = crys, &
!$!         nk = 1_i64, &
!$!         kvecs = kppathvecs, &
!$!         energies = el_ens_kp, &
!$!         evecs = el_evecs_kp, &
!$!         scissor = el%scissor)

       !Below, we will sum out m, n, and k
       do m = 1, wann%numwannbands

          ek = el%ens(ik, m)

          !Apply energy window to initial electron
          if(abs(ek - el%enref) > el%fsthick) cycle

          do iOmega = nOmegas/2 + 2, nOmegas
             do n = 1, wann%numwannbands

                !kvec = matmul(crys%reclattvecs, el%wavevecs(ik, :))
                !kpqvec = matmul(crys%reclattvecs, kppathvecs(1, :))
!$!                 overlap = (1.0_r64 + no*npo*dot_product(el%vels(ik, m, :), kpvel(1, n, :))&
!$!                                /twonorm(el%vels(ik, m, :))/twonorm(kpvel(1, n, :)))/2
                
                ekp = el_ens_kp(1, n)

                !Apply energy window to final electron
                if(abs(ekp - el%enref) > el%fsthick) cycle

                !This is |U(k')U^\dagger(k)|_nm squared
                !(Recall that U^\dagger(k) is the diagonalizer of the electronic hamiltonian.)
                overlap = (abs(dot_product(el_evecs_kp(1, n, :), el%evecs(ik, m, :))))**2
                ! overlap = 1.0_r64

!$!                 spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                      (Fermi(ek, el%chempot, crys%T) - &
!$!                      Fermi(ekp, el%chempot, crys%T))*overlap* &
!$!                      delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                      el%wvmesh, el%simplex_map, &
!$!                      el%simplex_count, el%simplex_evals)

!!$                spec_eps(iOmega) = spec_eps(iOmega) + &
!!$                     (Fermi(ek, el%chempot, crys%T) - &
!!$                     Fermi(ek + Omegas(iOmega), el%chempot, crys%T))*overlap* &
!!$                     delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!!$                     el%wvmesh, el%simplex_map, &
!!$                     el%simplex_count, el%simplex_evals)
!!$
!$!                 spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                      (Fermi(ekp - Omegas(iOmega), el%chempot, crys%T) - &
!$!                      Fermi(ekp, el%chempot, crys%T))*overlap* &
!$!                      delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                      el%wvmesh, el%simplex_map, &
!$!                      el%simplex_count, el%simplex_evals)
                ! Delta replaced by Gaussian
                spec_eps(iOmega) = spec_eps(iOmega) + &
                     (Fermi(ekp - Omegas(iOmega), el%chempot, crys%T) - &
                     Fermi(ekp, el%chempot, crys%T))*overlap* &
                     deltafunc_pol(ik, m, kpvel, ekp - ek - Omegas(iOmega), crys, el) 

!!$                spec_eps(iOmega) = spec_eps(iOmega) + &
!!$                     Fermi(ekp, el%chempot, crys%T)* &
!!$                     (1.0_r64 - Fermi(ekp - Omegas(iOmega), el%chempot, crys%T))/ &
!!$                     Bose(Omegas(iOmega), crys%T)* &
!!$                     overlap* &
!!$                     delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!!$                     el%wvmesh, el%simplex_map, &
!!$                     el%simplex_count, el%simplex_evals)
             end do
          end do
       end do
    end do

    spec_eps = spec_eps*el%spindeg/crys%volume*crys%thickness
    
    do iOmega = 1, nOmegas/2 ! negative sector
       !Recall that the resolvent is already normalized in the full wave vector mesh.
       !As such, the 1/product(el%wvmesh) is not needed in the expression below.
       spec_eps(iOmega) = -spec_eps(nOmegas + 1 -iOmega)
    end do
    !At this point [spec_eps] = nm^-2.eV^-1

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end subroutine spectral_head_polarizability_2d_qpath_old

  subroutine spectral_head_polarizability_2d_lowqpath_old(spec_eps, Omegas, qcrys, &
       el, wann, crys, tetrahedra)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! qvec Wave vector (fractional) of excitation in the electron gas
    !! el Electron data type

    real(r64), intent(in) :: Omegas(:), qcrys(3)
    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys

    real(r64), allocatable, intent(out) :: spec_eps(:)
    logical, intent(in) :: tetrahedra

    !Locals
    integer(i64) :: m, n, ik, iOmega, nOmegas, k_indvec(3), kp_indvec(3), no, npo
    real(r64) :: overlap, ek, ekp, delta, Omega_l, Omega_r, &
         el_ens_kp(1, el%numbands), kppathvecs(1, 3), kvec(3), kpqvec(3),&
         kpvel(1, el%numbands, 3), beta, fnk, qcart(3), kunit(3), delta_func
    complex(r64) :: el_evecs_kp(1, el%numbands, el%numbands)
    procedure(delta_fn), pointer :: delta_fn_ptr => null()
    real(r64) :: delta_sum, fermidiff_sum, fermidiff, qmag

    beta = 1.0_r64/kB/crys%T
    qcart = matmul(crys%reclattvecs, qcrys)
    qmag = twonorm(qcart)

    nOmegas = size(Omegas)

    allocate(spec_eps(nOmegas))

    !Associate delta function procedure pointer ! 2D
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra = .false.)
    !delta_fn_ptr => get_delta_fn_pointer(tetrahedra)

    spec_eps = 0.0
    fermidiff_sum = 0.0
    delta_sum = 0.0
    do ik = 1, el%nwv
       
       kppathvecs(1, :) = el%wavevecs(ik, :).umklapp.qcrys
       ! kunit = el%wavevecs(ik, :)/twonorm(el%wavevecs(ik, :))
       
       call wann%el_wann(crys = crys, &
        nk = 1_i64, &
        kvecs = kppathvecs, &
        energies = el_ens_kp, velocities = kpvel, &
        evecs = el_evecs_kp, &
        scissor = el%scissor)
!$!        call wann%el_wann(crys = crys, &
!$!         nk = 1_i64, &
!$!         kvecs = kppathvecs, &
!$!         energies = el_ens_kp, &
!$!         evecs = el_evecs_kp, &
!$!         scissor = el%scissor)

       !Below, we will sum out m, n, and k
       do m = 1, wann%numwannbands

          ek = el%ens(ik, m)

          !Apply energy window to initial electron
          if(abs(ek - el%enref) > el%fsthick) cycle

          do iOmega = nOmegas/2 + 2, nOmegas
             do n = 1, wann%numwannbands

                ekp = el_ens_kp(1, n)

                fnk = Fermi(el%ens(ik, n), el%chempot, crys%T)

                !Apply energy window to final electron
                if(abs(ekp - el%enref) > el%fsthick) cycle

                !This is |U(k')U^\dagger(k)|_nm squared
                !(Recall that U^\dagger(k) is the diagonalizer of the electronic hamiltonian.)
                overlap = (abs(dot_product(el_evecs_kp(1, n, :), el%evecs(ik, m, :))))**2
                ! overlap = 1.0_r64

!$!                 spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                      (Fermi(ek, el%chempot, crys%T) - &
!$!                      (fnk - beta*fnk*(1.0_r64 - fnk)*hbar_eVps*&
!$!                      dot_product(qcart, el%vels(ik, n, :))))*overlap* &
!$!                      delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                      el%wvmesh, el%simplex_map, &
!$!                      el%simplex_count, el%simplex_evals)
                
                ! Delta replaced by Gaussian
                delta_func = deltafunc_pol(ik, m, kpvel, ekp - ek - Omegas(iOmega), crys, el)
                fermidiff = Fermi(ek, el%chempot, crys%T) - &
                            (fnk - beta*fnk*(1.0_r64 - fnk)*hbar_eVps*&
                            dot_product(qcart, el%vels(ik, n, :)))
                if(Omegas(iOmega)>0.0606 .and. Omegas(iOmega)<0.1818) then
                   delta_sum = delta_sum + delta_func
                   fermidiff_sum = fermidiff_sum + fermidiff
                end if
                spec_eps(iOmega) = spec_eps(iOmega) + &
                                   fermidiff*overlap*delta_func

                !if (delta_func > 1e-5) print *,"Gauss overlap==>", delta_func
                !write(1000, *) delta_func
!$!                 if(Omegas(iOmega)>0.02) then 
!$!                   spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                        (Fermi(ek, el%chempot, crys%T) - &
!$!                        (fnk - beta*fnk*(1.0_r64 - fnk)*hbar_eVps*&
!$!                        dot_product(qcart, el%vels(ik, n, :))))*overlap* &
!$!                        delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                        el%wvmesh, el%simplex_map, &
!$!                        el%simplex_count, el%simplex_evals)
!$!                 else
!$!                   ! Delta replaced by 1, linear in omega
!$!                   spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                        (Fermi(ek, el%chempot, crys%T) - &
!$!                        (fnk - beta*fnk*(1.0_r64 - fnk)*Omegas(iOmega)))*&
!$!                        overlap/product(el%wvmesh)
!$!                 end if

!$!                 spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                      beta*fnk*(1.0_r64 - fnk)*hbar_eVps* &
!$!                      dot_product(qcart, el%vels(ik, n, :))*overlap* &
!$!                      delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                      el%wvmesh, el%simplex_map, &
!$!                      el%simplex_count, el%simplex_evals)

!$!                 spec_eps(iOmega) = spec_eps(iOmega) + &
!$!                      beta*fnk*(1.0_r64 - fnk)*hbar_eVps*0.3e3* &
!$!                      dot_product(qcart, kunit)*overlap* &
!$!                      delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!$!                      el%wvmesh, el%simplex_map, &
!$!                      el%simplex_count, el%simplex_evals)
             end do
          end do
       end do
    end do

    write(170, '(F8.3, 1X, E12.3, 1X, E12.3)') qmag, fermidiff_sum, delta_sum
    spec_eps = spec_eps*el%spindeg/crys%volume*crys%thickness
    
    do iOmega = 1, nOmegas/2 ! negative sector
       !Recall that the resolvent is already normalized in the full wave vector mesh.
       !As such, the 1/product(el%wvmesh) is not needed in the expression below.
       spec_eps(iOmega) = -spec_eps(nOmegas + 1 -iOmega)
    end do
    !At this point [spec_eps] = nm^-2.eV^-1

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end subroutine spectral_head_polarizability_2d_lowqpath_old


  pure function deltafunc_pol(ik, m, kpvel, en, crys, el) 
    integer(i64), intent(in) :: ik, m
    real(r64), intent(in) :: en, kpvel(3)
    type(electron), intent(in) :: el
    type(crystal), intent(in) :: crys
    real(r64) :: sigma, deltafunc_pol
   
    integer(i64) :: dim
    real (i64) :: onebyroot2pi, onebyroot12, Qs(2, 3), aux
    
    ! Compute Qs for the sigma Gaussian
    onebyroot2pi = 1.0_r64/sqrt(2.0*pi)
    onebyroot12 = 1.0_r64/sqrt(12.0_r64)
    do dim = 1, 2
       !Qs(dim, :) = crys%reclattvecs(dim, :)/el%wvmesh(dim)
       Qs(dim, :) = crys%reclattvecs(:, dim)/el%wvmesh(dim)
    end do
    ! Calculate adaptive smearing
    aux = 0.0_r64
    do dim = 1, 2
       aux = aux + &
             dot_product(el%vels(ik, m, :) - kpvel, Qs(dim, :))**2
    end do
    sigma = hbar_eVps*onebyroot12*sqrt(aux)
    
    !sigma = 1e-2_r64
    deltafunc_pol = max(onebyroot2pi/sigma/product(el%wvmesh)*&
                    exp(-0.5_r64*(en/sigma)**2), 1.0e-9_r64)
!$!     deltafunc_pol = onebyroot2pi/sigma/product(el%wvmesh)*&
!$!                     exp(-0.5_r64*(en/sigma)**2)
  end function deltafunc_pol

!!$  !DEBUG/TEST
!!$  subroutine calculate_RPA_dielectric_3d_G0_scratch(el, crys, num, wann)
!!$    !! ??
!!$    !!
!!$    !! el Electron data type
!!$    !! crys Crystal data type
!!$    !! num Numerics data type
!!$
!!$    type(electron), intent(in) :: el
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$    type(wannier), intent(in) :: wann
!!$
!!$    !Locals
!!$    real(r64), allocatable :: energylist(:), qlist(:, :), qmaglist(:)
!!$    real(r64) :: qcrys(3)
!!$    integer(i64) :: iq, iOmega, numomega, numq, &
!!$         start, end, chunk, num_active_images, qxmesh
!!$    real(r64), allocatable :: spec_X0(:), ImX0(:), ReX0(:)
!!$    complex(r64), allocatable :: diel(:, :)
!!$    character(len = 1024) :: filename
!!$    real(r64) :: omega_plasma
!!$
!!$    !Silicon
!!$    !omega_plasma = 1.0e-9_r64*hbar*sqrt(el%conc_el/perm0/crys%epsiloninf/(0.267*me)) !eV
!!$
!!$    !wGaN
!!$    !omega_plasma = 1.0e-9_r64*hbar*sqrt(el%conc_el/perm0/crys%epsiloninf/(0.259_r64*me)) !eV
!!$        
!!$    !if(this_image() == 1) then
!!$    !   print*, "plasmon energy = ", omega_plasma
!!$    !   print*, "epsilon infinity = ", crys%epsiloninf
!!$    !end if
!!$
!!$    !TEST
!!$    numq = el%wvmesh(1)
!!$    qxmesh = numq
!!$    !Create qlist in crystal coordinates
!!$    allocate(qlist(numq, 3), qmaglist(numq))
!!$    do iq = 1, numq
!!$       qlist(iq, :) = [(iq - 1.0_r64)/qxmesh, (iq - 1.0_r64)/qxmesh, 0.0_r64]
!!$       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
!!$    end do
!!$    call sort(qmaglist)
!!$
!!$    !Create energy grid
!!$    numomega = 601 !6 !5 !1001 !1001
!!$    allocate(energylist(numomega))
!!$    call linspace(energylist, -0.5_r64, 0.5_r64, numomega)
!!$    
!!$    !Allocate diel_ik to hold maximum possible Omega
!!$    allocate(diel(numq, numomega))
!!$    diel = 0.0_r64
!!$
!!$    !Allocate spectral and imaginary X0 (defined on uniform frequency mesh)
!!$    allocate(spec_X0(numomega), ImX0(numomega))
!!$    
!!$    !Distribute points among images
!!$    call distribute_points(numq, chunk, start, end, num_active_images)
!!$
!!$    if(this_image() == 1) then
!!$       write(*, "(A, I10)") " #q-vecs = ", numq
!!$       write(*, "(A, I10)") " #q-vecs/image <= ", chunk
!!$    end if
!!$
!!$    do iq = start, end !Over IBZ k points
!!$       qcrys = qlist(iq, :) !crystal coordinates
!!$
!!$       call spectral_head_polarizability_3d_q(&
!!$            spec_X0, energylist, qcrys, el, wann, crys, num%tetrahedra)
!!$       
!!$       ImX0 = -pi*spec_X0
!!$       call hilbert_transform(-ImX0, ReX0)
!!$          
!!$       !Calculate RPA dielectric (diagonal in G-G' space)
!!$       diel(iq, :) = crys%epsiloninf - &
!!$            1.0_r64/qmaglist(iq)**2* &
!!$            (ReX0 + oneI*ImX0)/perm0*qe*1.0e9_r64
!!$    end do
!!$
!!$    call co_sum(diel)
!!$    
!!$    !Print to file
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_qpath", qlist)
!!$    call write2file_rank1_real("RPA_dielectric_3D_G0_qmagpath", qmaglist)
!!$    call write2file_rank1_real("RPA_dielectric_3D_G0_Omega", energylist)
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_real", real(diel))
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_imag", imag(diel))
!!$  end subroutine calculate_RPA_dielectric_3d_G0_scratch

! debugging
  subroutine calculate_RPA_dielectric_2d_model(el, crys, num, wann)
    !! ??
    !!
    !! el Electron data type
    !! crys Crystal data type
    !! num Numerics data type

    type(electron), intent(in) :: el
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    type(wannier), intent(in) :: wann

    !Locals
    real(r64), allocatable :: energylist(:), qlist(:, :), qmaglist(:)
    real(r64) :: qcrys(3), qcart(3)
    integer(i64) :: iq, io, iOmega, numomega, numq, &
         start, end, chunk, num_active_images, qxmesh
    real(r64), allocatable :: spec_X0(:), ImX0(:), ReX0(:), Ls(:, :)
    complex(r64), allocatable :: diel_rpa(:, :), diel_tf(:, :), X0_qw(:)
    complex(r64), allocatable :: pol(:, :)
    character(len = 1024) :: filename
    real(r64) :: omega_plasma, prefac, q2norm, dim_norm
    real(r64), parameter :: plus0 = 1
    integer :: ik1, ik2, ik3
    real(r64) :: W_qw_msq, Gplusq(3), Gplusq_2norm
    real(r64), allocatable :: diel_qw(:) 

    !Silicon
    !omega_plasma = 1.0e-9_r64*hbar*sqrt(el%conc_el/perm0/crys%epsiloninf/(0.267*me)) !eV

    !wGaN
    !omega_plasma = 1.0e-9_r64*hbar*sqrt(el%conc_el/perm0/crys%epsiloninf/(0.259_r64*me)) !eV
        
    !if(this_image() == 1) then
    !   print*, "plasmon energy = ", omega_plasma
    !   print*, "epsilon infinity = ", crys%epsiloninf
    !end if

    !TEST
    numq = el%wvmesh(1)
    qxmesh = numq
    !Create qlist in crystal coordinates
    allocate(qlist(numq, 3), qmaglist(numq))
    do iq = 1, numq
       ! Gamma -> 1, 1, 0
       !qlist(iq, :) = [(iq - 1.0_r64)/qxmesh, (iq - 1.0_r64)/qxmesh, 0.0_r64]
       ! Gamma -> K (0.333, 0.333, 0)
       !qlist(iq, :) = [(iq - 1.0_r64)/(qxmesh - 1)/3, (iq - 1.0_r64)/(qxmesh - 1)/3, &
       !                 0.0_r64] 
       ! Gamma -> K (0.1, 0.1, 0)
       qlist(iq, :) = [(iq - 1.0_r64)/(qxmesh - 1)/10, (iq - 1.0_r64)/(qxmesh - 1)/10, &
                        0.0_r64] 
       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
    end do
    call sort(qmaglist)

    !Create energy grid
    numomega = num%ncont_mesh !601 !6 !5 !1001
    allocate(energylist(numomega))
    allocate(diel_qw(numomega))
    call linspace(energylist, -3.5_r64, 3.5_r64, numomega)
    
    !Allocate diel_ik to hold maximum possible Omega
    allocate(diel_rpa(numq, numomega), diel_tf(numq, numomega), Ls(numq, numomega))
    allocate(pol(numq, numomega))
    pol = 0.0_r64
    diel_rpa = 0.0_r64

    !Allocate spectral and imaginary X0 (defined on uniform frequency mesh)
    allocate(spec_X0(numomega), ImX0(numomega), ReX0(numomega))
    
    !Distribute points among images
    call distribute_points(numq, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #q-vecs = ", numq
       write(*, "(A, I10)") " #q-vecs/image <= ", chunk
    end if

    prefac = 1.0e9_r64*qe/perm0 ! ev.nm 
    if(crys%twod) prefac = prefac/2

    do iq = start, end !Over IBZ k points
       qcrys = qlist(iq, :) !crystal coordinates
       qcart = matmul(crys%reclattvecs, qcrys) !cartesian coordinates
       
       ! print *, "q->", iq
       ! |G + q|^2 or |G + q|, for 3D or 2D case
       q2norm = qmaglist(iq)**(crys%dim-1)
       !q2norm = qmaglist(iq)**2 

!$!        call spectral_head_polarizability_2d_qpath_old(&
!$!             spec_X0, energylist, qcrys, el, wann, crys, num%tetrahedra)
      
       !if(q2norm<1e-1) then
       if(qmaglist(iq)<crys%qTF/29.0_r64) then
          call spectral_head_polarizability_2d_lowqpath_old(&
            spec_X0, energylist, qcrys, el, wann, crys, num%tetrahedra)
       else
          call spectral_head_polarizability_2d_qpath_old(&
            spec_X0, energylist, qcrys, el, wann, crys, num%tetrahedra)
       end if

       ImX0 = -pi*spec_X0
       call hilbert_transform(-ImX0, ReX0)
       X0_qw = ReX0 + oneI*ImX0
       pol(iq, :) = X0_qw
          
!$!        W_qw_msq = 0.0_r64
!$!        diel_qw = 0.0_r64
!$!        do concurrent(ik1 = -1:1, ik2 = -1:1, ik3 = -1:1)
!$!           Gplusq = (ik1*crys%reclattvecs(:, 1) &
!$!                + ik2*crys%reclattvecs(:, 2) &
!$!                + ik3*crys%reclattvecs(:, 3)) + qcart
!$! 
!$!           !|G + q|^2 or |G + q|, for 3D or 2D case
!$!           Gplusq_2norm = twonorm(Gplusq)**(crys%dim-1)
!$! 
!$!           !Dielectric matrix elements 
!$!           diel_qw = diel_qw + prefac*X0_qw/Gplusq_2norm
!$! 
!$!           !Squared Coulomb matrix elements without the prefactor
!$!           !W_qw_msq = W_qw_msq + abs(1.0_r64/diel_qw/Gplusq_2norm)**2
!$!        end do
!$! 
!$!        !gCoul2_RPA = W_qw_msq*prefac**2*overlap/dim_norm ! eV^2
!$!        diel_rpa(iq, :) = 1.0_r64 - diel_qw

       !Calculate RPA dielectric (diagonal in G-G' space)
!$!        diel_rpa(iq, :) = crys%epsiloninf - &
!$!             1.0_r64/qmaglist(iq)**2* &
!$!             (ReX0 + oneI*ImX0)/perm0*qe*1.0e9_r64

       !k_star = 4.436
       !diel_rpa(iq, :) = 1.0_r64 - prefac*X0_qw/q2norm
       diel_rpa(iq, :) = crys%epsiloninf - prefac*X0_qw/q2norm
       !diel_tf(iq, :) = crys%epsiloninf + crys%qTF/q2norm + oneI*energylist*plus0
       !Ls(iq, :) = -prefac*q2norm*ImX0/((q2norm*crys%epsiloninf - &
       !  prefac*ReX0)**2 + (prefac*ImX0)**2)
       !if(iq==1) Ls(iq, :) = 1e-20_r64
       !do io=1, numomega
       !  if(Ls(iq, io)<1e-50_r64) Ls(iq, io) = 1e-50_r64
       !end do

    end do

    call co_sum(pol)
    call co_sum(diel_rpa)
    !call co_sum(diel_tf)
    !call co_sum(Ls)
    
    !Print to file
    call write2file_rank2_real("RPA_dielectric_2D_G0_qpath", qlist)
    call write2file_rank1_real("RPA_dielectric_2D_G0_qmagpath", qmaglist)
    call write2file_rank1_real("RPA_dielectric_2D_G0_Omega", energylist)
    call write2file_rank2_real("RPA_polarisability_2D_G0_real", real(pol))
    call write2file_rank2_real("RPA_polarisability_2D_G0_imag", imag(pol))
    call write2file_rank2_real("RPA_dielectric_2D_G0_real", real(diel_rpa))
    call write2file_rank2_real("RPA_dielectric_2D_G0_imag", imag(diel_rpa))
    !call write2file_rank2_real("TF_dielectric_2D_G0_real", real(diel_tf))
    !call write2file_rank2_real("TF_dielectric_2D_G0_imag", imag(diel_tf))
    !call write2file_rank2_real("RPA_dielectric_2D_G0_loss", Ls)
  end subroutine calculate_RPA_dielectric_2d_model
end module screening_module

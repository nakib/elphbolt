module screening_module

  use precision, only: r64, i64
  use params, only: kB, qe, pi, perm0, oneI, hbar, me
  use electron_module, only: electron
  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use misc, only: linspace, mux_vector, binsearch, Fermi, print_message, &
       compsimps, twonorm, write2file_rank2_real, write2file_rank1_real, &
       distribute_points, sort, qdist
  use delta, only: delta_fn, get_delta_fn_pointer

  implicit none

  private
  public calculate_qTF, calculate_RPA_dielectric_3d_G0_qpath
  
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

    if(crys%epsilon0 /= 0) then
       call print_message("Calculating Thomas-Fermi screening wave vector...")

       do ib = 1, el%numbands
          do ik = 1, el%nwv
             fFD = Fermi(el%ens(ik, ib), el%chempot, crys%T)
             crys%qTF = crys%qTF + fFD*(1.0_r64 - fFD)
          end do
       end do

       !Free-electron gas Thomas-Fermi model
       ! qTF**2 = spindeg*e^2*beta/nptq/vol_pcell/perm0/epsilon0*Sum_{BZ}f0_{k}(1-f0_{k})
       crys%qTF = sqrt(1.0e9_r64*crys%qTF*el%spindeg*beta*qe**2/product(el%wvmesh)&
            /crys%volume/perm0/crys%epsilon0) !nm^-1

       if(this_image() == 1) then
          write(*, "(A, 1E16.8, A)") ' Thomas-Fermi screening wave vector = ', crys%qTF, ' 1/nm'
       end if
    end if
  end subroutine calculate_qTF

  pure elemental real(r64) function finite_element(w, wl, w0, wr)
    !! Triangular finite elements.
    !!
    !! w Continuous sampling variable
    !! wl, w0, wr 3-point stencil left, center, right, respectively

    real(r64), intent(in) :: w
    real(r64), intent(in) :: wl, w0, wr
    
    if(wl <= w .and. w <= w0) then
       finite_element = (w - wl)/(w0 - wl)
    else if(w0 <= w .and. w <= wr) then
       finite_element = (wr - w)/(wr - w0)
    else
       finite_element = 0.0_r64
    end if
  end function finite_element

  subroutine calculate_Hilbert_weights(w_disc, w_cont, zeroplus, Hilbert_weights)
    !! Calculator of the weights for the Hilbert-Kramers-Kronig transform.
    !!
    !! w_disc Discrete sampling points
    !! w_cont Continuous sampling variable
    !! zeroplus A small positive number
    !! Hilbert_weights The Hilbert weights

    real(r64), intent(in) :: w_disc(:), w_cont(:), zeroplus
    complex(r64), allocatable, intent(out) :: Hilbert_weights(:, :)

    !Locals
    integer :: n_disc, n_cont, i, j
    real(r64) :: phi_i(size(w_cont)), dw, realpart, imagpart
    complex(r64) :: integrand(size(w_cont))

    n_disc = size(w_disc)
    n_cont = size(w_cont)

    allocate(Hilbert_weights(n_disc, n_disc))

    !Grid spacing of "continuous" variable
    dw = w_cont(2) - w_cont(1)    

    do i = 1, n_disc - 1
       if(i == 1) then
          Phi_i = finite_element(w_cont, w_disc(i) - dw, w_disc(i), w_disc(i + 1))
       else if(i == n_disc) then
          Phi_i = finite_element(w_cont, w_disc(i - 1), w_disc(i), w_disc(i) + dw)
       else
          Phi_i = finite_element(w_cont, w_disc(i - 1), w_disc(i), w_disc(i + 1))
       end if
       do j = 1, n_disc
!!$          integrand = Phi_i*(&
!!$               (w_disc(j) - w_cont)/((w_disc(j) - w_cont)**2 + zeroplus) - &
!!$               (w_disc(j) + w_cont)/((w_disc(j) + w_cont)**2 + zeroplus))

!!$          integrand = Phi_i*2.0_r64*w_cont/&
!!$               (w_disc(j)**2 - w_cont**2 - 2.0_r64*w_cont*oneI*zeroplus)

!!$          integrand = Phi_i/(w_disc(j)**2 - w_cont**2)

!!$          integrand = Phi_i*(&
!!$               1.0_r64/(w_disc(j) - w_cont - oneI*zeroplus) - &
!!$               1.0_r64/(w_disc(j) + w_cont + oneI*zeroplus))

          !DBG
          integrand = Phi_i*(&
               1.0_r64/(-w_disc(j) + w_cont - oneI*zeroplus) - &
               1.0_r64/(-w_disc(j) - w_cont + oneI*zeroplus))
          
          !TODO Would be nice to have a compsimps function instead of
          !a subroutine.
          !call compsimps(integrand, dw, Hilbert_weights(j, i))
          call compsimps(real(integrand), dw, realpart)
          call compsimps(imag(integrand), dw, imagpart)
          Hilbert_weights(j, i) = realpart + oneI*imagpart
       end do
    end do
  end subroutine calculate_Hilbert_weights

  subroutine head_polarizability_3d_T(eps_T, Omegas, spec_eps_T, Hilbert_weights_T)
    !! Head of the bare polarizability of the 3d Kohn-Sham system using
    !! Hilbert transform for a given set of temperature-dependent quantities.
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! eps_T Real part of bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! spec_eps_T Spectral head of bare polarizability
    !! Hilbert_weights_T Hilbert transform weights
    
    complex(r64), allocatable, intent(out) :: eps_T(:)
    real(r64), intent(in) :: Omegas(:)
    real(r64), intent(in) :: spec_eps_T(:)
    complex(r64), intent(in) :: Hilbert_weights_T(:, :)

    allocate(eps_T(size(Omegas)))
    
    !TODO Can optimize this sum with blas
    eps_T = matmul(Hilbert_weights_T, spec_eps_T)
  end subroutine head_polarizability_3d_T
  
  !subroutine Im_head_polarizability_3d(Imeps, Omegas, q_indvec, el, pcell_vol, T)
  subroutine spectral_head_polarizability_3d(spec_eps, Omegas, q_indvec, el, pcell_vol, T)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! q_indvec Wave vector of excitation in the electron gas (0-based integer triplet)
    !! el Electron data type
    !! pcell_vol Primitive unit cell volume
    !! T Temperature (K)

    real(r64), intent(in) :: Omegas(:), pcell_vol, T
    integer(i64), intent(in) :: q_indvec(3)
    type(electron), intent(in) :: el
    real(r64), allocatable, intent(out) :: spec_eps(:)
    
    !Locals
    integer(i64) :: m, n, ik, ikp, ikp_window, iOmega, nOmegas, k_indvec(3), kp_indvec(3)
    real(r64) :: overlap, ek, ekp, dOmega, delta, Omega_l, Omega_r
    procedure(delta_fn), pointer :: delta_fn_ptr => null()

    nOmegas = size(Omegas)

    dOmega = Omegas(2) - Omegas(1)
    
    allocate(spec_eps(nOmegas))

    !TODO don't have to use tetrahedron since I only need the imag part.
    !May use triangles also => easily extensible to the 2D case
    !
    !Associate delta function procedure pointer
    !delta_fn_ptr => get_delta_fn_pointer(tetrahedra = .true.)
    !delta_fn_ptr => get_delta_fn_pointer(tetrahedra = .false.)
    !Comment: The delta-fn methods above are giving very different numbers.
    
    spec_eps = 0.0
    do iOmega = 1, nOmegas
       !Below, we will sum out m, n, and k
       do m = 1, el%numbands
          do ik = 1, el%nwv
             ek = el%ens(ik, m)

             !Apply energy window to initial electron
             if(abs(ek - el%enref) > el%fsthick) cycle

             !Calculate index vector of final electron: k' = k + q
             k_indvec = nint(el%wavevecs(ik, :)*el%wvmesh)    
             kp_indvec = modulo(k_indvec + q_indvec, el%wvmesh) !0-based index vector

             !DEBUG: Disallow Umklapp
             !if(any(k_indvec/el%wvmesh >= 1.0)) cycle
             
             !Multiplex kp_indvec
             ikp = mux_vector(kp_indvec, el%wvmesh, 0_i64)

             !Check if final electron wave vector is within FBZ blocks
             call binsearch(el%indexlist, ikp, ikp_window)
             if(ikp_window < 0) cycle

             do n = 1, el%numbands
                ekp = el%ens(ikp_window, n)

                !Apply energy window to final electron
                if(abs(ekp - el%enref) > el%fsthick) cycle

                !This is |U(k')U^\dagger(k)|_nm squared
                !(Recall that U^\dagger(k) is the diagonalizer of the electronic hamiltonian.)
                overlap = (abs(dot_product(el%evecs(ikp_window, n, :), el%evecs(ik, m, :))))**2
                !overlap = 1.0

!!$                spec_eps(iOmega) = spec_eps(iOmega) + &
!!$                     (Fermi(ek, el%chempot, T) - &
!!$                     Fermi(ekp, el%chempot, T))*overlap* &
!!$                     delta_fn_ptr(ekp - Omegas(iOmega), ik, m, &
!!$                     el%wvmesh, el%simplex_map, &
!!$                     el%simplex_count, el%simplex_evals)
                
                
                !DBG
                Omega_l = Omegas(iOmega) - dOmega
                Omega_r = Omegas(iOmega) + dOmega
                
                if(Omega_l < ekp - ek .and. ekp - ek < Omega_r) continue
                
                delta = finite_element(ekp - ek, &
                     Omega_l, Omegas(iOmega), Omega_r)
!!$                
!!$                if(iOmega == 1) then
!!$                   delta = finite_element(ekp - ek, &
!!$                        Omegas_l, Omegas(iOmega), Omegas(iOmega + 1))
!!$                else if(iOmega == nOmegas) then
!!$                   delta = finite_element(ekp - ek, &
!!$                        Omegas(iOmega - 1), Omegas(iOmega), Omegas(iOmega) + dOmega)
!!$                else
!!$                   delta = finite_element(ekp - ek, &
!!$                        Omegas(iOmega - 1), Omegas(iOmega), Omegas(iOmega + 1))
!!$                end if
                spec_eps(iOmega) = spec_eps(iOmega) + &
                     (Fermi(ek, el%chempot, T) - &
                     Fermi(ekp, el%chempot, T))*overlap*delta/product(el%wvmesh)
             end do
          end do
       end do
    end do

    do iOmega = 1, nOmegas
       !Recall that the resolvent is already normalized in the full wave vector mesh.
       !As such, the 1/product(el%wvmesh) is not needed in the expression below.
!!$       spec_eps(iOmega) = &
!!$            -spec_eps(iOmega)*pi*sign(1.0_r64, Omegas(iOmega))*el%spindeg/pcell_vol

       spec_eps(iOmega) = &
            -spec_eps(iOmega)*sign(1.0_r64, Omegas(iOmega))*el%spindeg/pcell_vol

       !Zero out extremely small numbers !TODO What is small? In what units?
       !if(abs(spec_eps(iOmega)) < 1.0e-30_r64) spec_eps(iOmega) = 0.0_r64
    end do
    !At this point [spec_eps] = nm^-3.eV^-1

    if(associated(delta_fn_ptr)) nullify(delta_fn_ptr)
  end subroutine spectral_head_polarizability_3d

  !DEBUG/TEST
  subroutine calculate_RPA_dielectric_3d_G0_qpath(el, crys, num)
    !! ??
    !!
    !! el Electron data type
    !! crys Crystal data type
    !! num Numerics data type

    type(electron), intent(in) :: el
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num

    !Locals
    real(r64), allocatable :: energylist(:), qlist(:, :), qmaglist(:)
    real(r64) :: qcrys(3), zeroplus
    integer(i64) :: iq, iOmega, numomega, numq, &
         start, end, chunk, num_active_images, qxmesh
    real(r64), allocatable :: spec_eps(:)
    complex(r64), allocatable :: diel(:, :), eps(:), Hilbert_weights(:, :)
    character(len = 1024) :: filename
    real(r64) :: eps0_q0_prefac, omega_plasma

    !Silicon
    !omega_plasma = 1.0e-9*hbar*sqrt(el%conc_el/perm0/crys%epsilon0/(0.267*me)) !eV

    !wGaN
    omega_plasma = 1.0e-9*hbar*sqrt(el%conc_el/perm0/crys%epsiloninf/(0.22*me)) !eV
    
    if(this_image() == 1) then
       print*, "plasmon energy = ", omega_plasma
    end if

    !TEST
    !Material: Si
!!$    numq = el%wvmesh(1)
!!$    qxmesh = el%wvmesh(1) !numq
!!$    !Create qlist in crystal coordinates
!!$    allocate(qlist(numq, 3), qmaglist(numq))
!!$    do iq = 1, numq
!!$       qlist(iq, :) = [(iq - 1.0_r64)/qxmesh, (iq - 1.0_r64)/qxmesh, 0.0_r64]
!!$       !qlist(iq, :) = [(iq - 1.0_r64)/qxmesh - 0.5, (iq - 1.0_r64)/qxmesh - 0.5, 0.0_r64]
!!$       !qmaglist(iq) = qdist(qlist(iq, :), crys%reclattvecs)
!!$       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
!!$    end do
!!$    !call sort(qmaglist)

!!$    numq = el%nwv_irred
!!$    !Create qlist in crystal coordinates
!!$    allocate(qlist(numq, 3), qmaglist(numq))
!!$    do iq = 1, numq
!!$       qlist = el%wavevecs_irred
!!$       !qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
!!$       qmaglist(iq) = qdist(qlist(iq, :), crys%reclattvecs)
!!$    end do
!!$    call sort(qmaglist)

    numq = 2
    !Create qlist in crystal coordinates
    allocate(qlist(numq, 3), qmaglist(numq))
    do iq = 1, numq
       qlist(iq, :) = el%wavevecs_irred(iq, :)
       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
       !qmaglist(iq) = qdist(qlist(iq, :), crys%reclattvecs)
    end do

    !Create energy grid
    numomega = 200
    allocate(energylist(numomega))
    call linspace(energylist, 0.0_r64, 0.2_r64, numomega)

    !Small number
    zeroplus = (energylist(2) - energylist(1))*0.1
    
    !Allocate diel_ik to hold maximum possible Omega
    allocate(diel(numq, numomega))
    diel = 0.0_r64

    !Allocate polarizability
    allocate(spec_eps(numomega), eps(numomega))

    !Allocate the Hilbert weights
    allocate(Hilbert_weights(numomega, numomega))
    
    !Distribute points among images
    call distribute_points(numq, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #q-vecs = ", numq
       write(*, "(A, I10)") " #q-vecs/image <= ", chunk
    end if

    do iq = start, end !Over IBZ k points
       qcrys = qlist(iq, :) !crystal coordinates

       !TODO Check both limits:
       !1. q -> 0
       !2. Omega -> 0
          call spectral_head_polarizability_3d(&
               spec_eps, energylist, nint(qcrys*el%wvmesh)+0_i64, el, crys%volume, crys%T)

          !Calculate re_eps with Hilbert-Kramers-Kronig transform
          call calculate_Hilbert_weights(&
               w_disc = energylist, &
               w_cont = energylist, &
               zeroplus = zeroplus, & !Can this magic "small" number be removed?
               Hilbert_weights = Hilbert_weights)

          !call Re_head_polarizability_3d_T(Reeps, energylist, Imeps, Hilbert_weights)
          call head_polarizability_3d_T(eps, energylist, spec_eps, Hilbert_weights)
          
          !Calculate RPA dielectric (diagonal in G-G' space)
!!$          diel(iq, :) = 1.0_r64 - &
!!$               (1.0_r64/twonorm(matmul(crys%reclattvecs, qcrys)))**2* &
!!$               eps/perm0*qe*1.0e9_r64

!!$          diel(iq, :) = 1.0_r64 - &
!!$               (1.0_r64/twonorm(matmul(crys%reclattvecs, qcrys)))**2* &
!!$               (real(eps) + oneI*pi*spec_eps)/perm0*qe*1.0e9_r64

!!$          diel(iq, :) = 1.0_r64 - &
!!$               1.0_r64/qmaglist(iq)**2* &
!!$               (real(eps) + oneI*pi*spec_eps)/perm0*qe*1.0e9_r64

       if(all(qcrys == 0.0_r64)) then
          !DEBUG: This is not a generally computable limit since the plasmon
          !energy expression requires a single effective mass.
          !Just leaving it here for now to get the plasmon peak in the loss function.
          diel(iq, 1:numomega) = 1.0_r64 - &
               (omega_plasma/energylist(1:numomega))**2 !+ oneI*1.0e-9
       else
          !DBG
          diel(iq, :) = crys%epsiloninf - &
               1.0_r64/qmaglist(iq)**2* &
               (real(eps) + oneI*pi*(spec_eps))/perm0*qe*1.0e9_r64
!!$          
!!$          diel(iq, :) = crys%epsiloninf - &
!!$               1.0_r64/qmaglist(iq)**2*eps/perm0*qe*1.0e9_r64
       end if
    end do

    call co_sum(diel)

    !Handle Omega = 0 case
    !diel(:, 1) = 1.0_r64 + (crys%qTF/qmaglist(:))**2
    
    !Print to file
    call write2file_rank2_real("RPA_dielectric_3D_G0_qpath", qlist)
    call write2file_rank1_real("RPA_dielectric_3D_G0_qmagpath", qmaglist)
    call write2file_rank1_real("RPA_dielectric_3D_G0_Omega", energylist)
    call write2file_rank2_real("RPA_dielectric_3D_G0_real", real(diel))
    call write2file_rank2_real("RPA_dielectric_3D_G0_imag", imag(diel))
  end subroutine calculate_RPA_dielectric_3d_G0_qpath
  
!!$  subroutine calculate_RPA_dielectric_3d(el, crys, num)
!!$    !! Dielectric function of the 3d Kohn-Sham system in the
!!$    !! random-phase approximation (RPA) using Eq. B1 and B2
!!$    !! of Knapen, Kozaczuk and Lin Phys. Rev. D 104, 015031 (2021).
!!$    !!
!!$    !! Here we calculate the diagonal in G-G' space. Moreover,
!!$    !! we use the approximation G.r -> 0.
!!$    !!
!!$    !! el Electron data type
!!$    !! crys Crystal data type
!!$    !! num Numerics data type
!!$
!!$    type(electron), intent(in) :: el
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$
!!$    !Locals
!!$    real(r64) :: Omega, k(3), kp(3), q(3), ek, ekp, Gplusq_crys(3)
!!$    integer(i64) :: ik, ikp, m, n, &
!!$         k_indvec(3), kp_indvec(3), q_indvec(3), count, &
!!$         start, end, chunk, num_active_images
!!$    integer :: ig1, ig2, ig3, igmax, num_gmax
!!$    complex(r64) :: polarizability
!!$    complex(r64), allocatable :: diel_ik(:)
!!$    character(len = 1024) :: filename
!!$
!!$    print*, el%numbands
!!$
!!$    !This sets how large the G vectors can be. Choosing 3
!!$    !means including -3G to +3G in the calculations. This
!!$    !should be a safe range.
!!$    igmax = 3
!!$    num_gmax = (2*igmax + 1)**3
!!$
!!$    !Allocate diel_ik to hold maximum possible Omega x G points
!!$    allocate(diel_ik(num_gmax*el%nstates_inwindow*el%numbands))
!!$
!!$    !Distribute points among images
!!$    call distribute_points(el%nwv_irred, chunk, start, end, num_active_images)
!!$
!!$    if(this_image() == 1) then
!!$       write(*, "(A, I10)") " #k-vecs = ", el%nwv_irred
!!$       write(*, "(A, I10)") " #k-vecs/image <= ", chunk
!!$    end if
!!$
!!$    do ik = start, end !Over IBZ k points
!!$       !Initiate counter for Omega x G points
!!$       count = 0
!!$
!!$       !Initial (IBZ blocks) wave vector (crystal coords.)
!!$       k = el%wavevecs_irred(ik, :)
!!$
!!$       !Convert from crystal to 0-based index vector
!!$       k_indvec = nint(k*el%wvmesh)
!!$
!!$       do m = 1, el%numbands
!!$          !IBZ electron energy
!!$          ek = el%ens_irred(ik, m)
!!$
!!$          !Check energy window
!!$          if(abs(ek - el%enref) > el%fsthick) cycle
!!$
!!$          do ikp = 1, el%nwv
!!$             !Final wave vector (crystal coords.)
!!$             kp = el%wavevecs(ikp, :)
!!$
!!$             !Convert from crystal to 0-based index vector
!!$             kp_indvec = nint(kp*el%wvmesh)
!!$
!!$             !Calculate q_indvec (folded back to the 1BZ)
!!$             q_indvec = modulo(kp_indvec - k_indvec, el%wvmesh) !0-based index vector
!!$
!!$             !Calculate q_indvec (without folding back to the 1BZ)
!!$             !q_indvec = kp_indvec - k_indvec !0-based index vector
!!$
!!$             do n = 1, el%numbands
!!$                ekp = el%ens(ikp, n)
!!$
!!$                !Check energy window
!!$                if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle
!!$
!!$                !Electron-hole pair/plasmon energy
!!$                Omega = ekp - ek
!!$
!!$                !Calculate RPA polarizability
!!$                polarizability = &
!!$                     RPA_polarizability_3d(Omega, q_indvec, el, crys%volume, crys%T)
!!$
!!$                do ig1 = -igmax, igmax
!!$                   do ig2 = -igmax, igmax
!!$                      do ig3 = -igmax, igmax
!!$                         !Counter for Omega x G points
!!$                         count = count + 1
!!$
!!$                         !Calculate G+q
!!$                         Gplusq_crys = ([ig1, ig2, ig3] + q_indvec)/dble(el%wvmesh)
!!$
!!$                         !Calculate RPA dielectric (diagonal in G-G' space)
!!$                         diel_ik(count) = 1.0_r64 - &
!!$                              (1.0_r64/twonorm(matmul(crys%reclattvecs, Gplusq_crys)))**2* &
!!$                              polarizability/perm0*qe*1.0e9_r64
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       !Change to data output directory
!!$       call chdir(trim(adjustl(num%epsilondir)))
!!$
!!$       !Write data in binary format
!!$       !Note: this will overwrite existing data!
!!$       write (filename, '(I9)') ik
!!$       filename = 'epsilon_RPA.ik'//trim(adjustl(filename))
!!$       open(1, file = trim(filename), status = 'replace', access = 'stream')
!!$       write(1) count
!!$       write(1) diel_ik(1:count)
!!$       close(1)
!!$    end do
!!$
!!$    sync all
!!$  end subroutine calculate_RPA_dielectric_3d
!!$
!!$  !DEBUG/TEST
!!$  subroutine calculate_RPA_dielectric_3d_G0_qpath(el, crys, num)
!!$    !! Dielectric function of the 3d Kohn-Sham system in the
!!$    !! random-phase approximation (RPA) using Eq. B1 and B2
!!$    !! of Knapen, Kozaczuk and Lin Phys. Rev. D 104, 015031 (2021).
!!$    !!
!!$    !! Here we calculate the diagonal in G-G' space. Moreover,
!!$    !! we use the approximation G.r -> 0.
!!$    !!
!!$    !! el Electron data type
!!$    !! crys Crystal data type
!!$    !! num Numerics data type
!!$
!!$    type(electron), intent(in) :: el
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$
!!$    !Locals
!!$    real(r64), allocatable :: energylist(:), qlist(:, :), qmaglist(:)
!!$    real(r64) :: qcrys(3)
!!$    integer(i64) :: iq, iOmega, numomega, numq, &
!!$         start, end, chunk, num_active_images, qxmesh
!!$    complex(r64) :: polarizability
!!$    complex(r64), allocatable :: diel(:, :)
!!$    character(len = 1024) :: filename
!!$    real(r64) :: a0, eps0_q0_prefac, omega_plasma
!!$
!!$    !a0 = 1.0e-9_r64*bohr2nm !Bohr radius in m
!!$    !eps0_q0_prefac = (4.0_r64/9.0_r64/pi)/a0* &
!!$    !     (3.0_r64*abs(sum(el%conc))*1.0e6_r64)**(-1.0_r64/3.0_r64) !
!!$    omega_plasma = 0.5025125628E-01 !eV
!!$
!!$    !TEST
!!$    !Material: Si
!!$    !Uniform energy mesh from 0-1.0 eV with uniform q-vecs in Gamma-Gamma along x
!!$    numomega = 200
!!$    numq = 50
!!$    qxmesh = 200
!!$    !Create qlist in crystal coordinates
!!$    allocate(qlist(numq, 3), qmaglist(numq))
!!$    do iq = 1, numq
!!$       qlist(iq, :) = [(iq - 1.0_r64)/qxmesh, 0.0_r64, 0.0_r64]
!!$       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
!!$    end do
!!$
!!$    !Create energy grid
!!$    allocate(energylist(numomega))
!!$    call linspace(energylist, 0.0_r64, 0.5_r64, numomega)
!!$
!!$    !Allocate diel_ik to hold maximum possible Omega
!!$    allocate(diel(numq, numomega))
!!$    diel = 0.0_r64
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
!!$       if(all(qcrys == 0.0_r64)) then
!!$          diel(iq, :) = 1.0_r64 - (omega_plasma/energylist)**2 + 1.0e-3*oneI
!!$       else
!!$          do iOmega = 1, numomega
!!$             !Calculate RPA polarizability
!!$             polarizability = &
!!$                  RPA_polarizability_3d(energylist(iOmega), &
!!$                  nint(qcrys*numq)+0_i64, el, crys%volume, crys%T)
!!$
!!$             !Calculate RPA dielectric (diagonal in G-G' space)
!!$             diel(iq, iOmega) = 1.0_r64 - &
!!$                  (1.0_r64/twonorm(matmul(crys%reclattvecs, qcrys)))**2* &
!!$                  polarizability/perm0*qe*1.0e9_r64
!!$          end do
!!$       end if
!!$    end do
!!$
!!$    sync all
!!$    call co_sum(diel)
!!$    sync all
!!$
!!$    !Print to file
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_qpath", qlist)
!!$    call write2file_rank1_real("RPA_dielectric_3D_G0_qmagpath", qmaglist)
!!$    call write2file_rank1_real("RPA_dielectric_3D_G0_Omega", energylist)
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_real", real(diel)*1.0_r64)
!!$    call write2file_rank2_real("RPA_dielectric_3D_G0_imag", imag(diel)*1.0_r64)
!!$  end subroutine calculate_RPA_dielectric_3d_G0_qpath
  !!
end module screening_module

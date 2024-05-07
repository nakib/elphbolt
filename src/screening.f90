module screening_module

  use precision, only: r64, i64
  use params, only: kB, qe, perm0
  use electron_module, only: electron
  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use misc, only: linspace, mux_vector, binsearch, Fermi, print_message
  use Green_function, only: resolvent

  implicit none

  private
  public calculate_qTF
  
contains
  
  subroutine calculate_qTF(crys, el)
    !! Calculate Thomas-Fermi screening wave vector from the static
    !! limit of the Lindhard function.
    !
    !Captain's log May 7, 2024. I would like to turn this into a pure function.
    !Why we even need crystal to have qTF as a member? It is only ever accessed by gchimp2.
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
  
  complex(r64) function Im_polarizability_3d(Omega, q_indvec, el, pcell_vol, T)
    !! Imaginary part of bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 18 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! Omega Energy of excitation in the electron gas
    !! q_indvec Wave vector of excitation in the electron gas (0-based integer triplet)
    !! el Electron data type
    !! pcell_vol Primitive unit cell volume
    !! T Temperature (K)

    real(r64), intent(in) :: Omega, pcell_vol, T
    integer(i64), intent(in) :: q_indvec(3)
    type(electron), intent(in) :: el

    !Locals
    integer(i64) :: m, n, ik, ikp, ikp_window, k_indvec(3), kp_indvec(3)
    real(r64) :: overlap, ek, ekp
    complex(r64) :: aux

    !Below, we will sum out m, n, and k, accumulating the result into aux.
    aux = 0.0
    do m = 1, el%numbands
       do ik = 1, el%nwv
          ek = el%ens(ik, m)

          !Apply energy window to initial electron
          if(abs(ek - el%enref) > el%fsthick) cycle

          !Calculate index vector of final electron: k' = k + q
          k_indvec = nint(el%wavevecs(ik, :)*el%wvmesh)
          kp_indvec = modulo(k_indvec + q_indvec, el%wvmesh) !0-based index vector

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

             aux = aux + &
                  (Fermi(ek, el%chempot, T) - &
                  Fermi(ekp, el%chempot, T))* &
                  resolvent(el, m, ik, ekp - Omega)*overlap
          end do
       end do
    end do
    !Recall that the resolvent is already normalized in the full wave vector mesh.
    !As such, the 1/product(el%wvmesh) is not needed in the expression below.
    aux = aux*sign(1.0_r64, omega)*el%spindeg/pcell_vol

    !Zero out extremely small numbers
    Im_polarizability_3d = imag(aux)
    if(abs(Im_polarizability_3d) < 1.0e-30) Im_polarizability_3d = 0.0_r64
  end function Im_polarizability_3d

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

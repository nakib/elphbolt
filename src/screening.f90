module screening_module

  use precision, only: r64, i64
  use params, only: kB, qe, pi, perm0, oneI, hbar, hbar_evps, me
  use electron_module, only: electron
  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use misc, only: linspace, mux_vector, binsearch, Fermi, print_message, &
       compsimps, twonorm, write2file_rank2_real, write2file_rank1_real, &
       distribute_points, sort, qdist, operator(.umklapp.), Bose, &
       Hilbert_transform
  use wannier_module, only: wannier
  use delta, only: delta_fn, get_delta_fn_pointer

  implicit none

  private
  public calculate_qTF, calculate_RPA_dielectric_3d_G0_scratch
  
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
       ! qTF**2 = spindeg*e^2*beta/nptq/vol_pcell/perm0*Sum_{BZ}f0_{k}(1-f0_{k})
       crys%qTF = sqrt(1.0e9_r64*crys%qTF*el%spindeg*beta*qe**2/product(el%wvmesh)&
            /crys%volume/perm0) !nm^-1

       if(this_image() == 1) then
          write(*, "(A, 1E16.8, A)") ' Thomas-Fermi screening wave vector = ', crys%qTF, ' 1/nm'
       end if
    end if
  end subroutine calculate_qTF

!!$  subroutine head_polarizability_real_3d_T(Reeps_T, Omegas, spec_eps_T, Hilbert_weights_T)
!!$    !! Head of the bare real polarizability of the 3d Kohn-Sham system using
!!$    !! Hilbert transform for a given set of temperature-dependent quantities.
!!$    !!
!!$    !! Here we calculate the diagonal in G-G' space. Moreover,
!!$    !! we use the approximation G.r -> 0.
!!$    !!
!!$    !! Reeps_T Real part of bare polarizability
!!$    !! Omega Energy of excitation in the electron gas
!!$    !! spec_eps_T Spectral head of bare polarizability
!!$    !! Hilbert_weights_T Hilbert transform weights
!!$    
!!$    real(r64), allocatable, intent(out) :: Reeps_T(:)
!!$    real(r64), intent(in) :: Omegas(:)
!!$    real(r64), intent(in) :: spec_eps_T(:)
!!$    real(r64), intent(in) :: Hilbert_weights_T(:, :)
!!$    
!!$    allocate(Reeps_T(size(Omegas)))
!!$    
!!$    !TODO Can optimize this sum with blas
!!$    Reeps_T = matmul(Hilbert_weights_T, spec_eps_T)
!!$  end subroutine head_polarizability_real_3d_T

  subroutine head_polarizability_imag_3d_T(Imeps_T, Omegas_samp, Omegas_cont, spec_eps_T)
    !! Head of the bare Imagniary polarizability of the 3d Kohn-Sham system using
    !! linear interpolation from the same evaluated on a fine, continuous mesh.
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! Imeps_T Imaginary part of bare polarizability
    !! Omegas_samp Sampling energies of excitation in the electron gas
    !! Omegas_cont Continous energies of excitation in the electron gas
    !! spec_eps_T Spectral head of bare polarizability
    
    real(r64), intent(in) :: Omegas_samp(:), Omegas_cont(:)
    real(r64), intent(in) :: spec_eps_T(:)
    real(r64), allocatable, intent(out) :: Imeps_T(:)

    integer :: isamp, ncont, nsamp, ileft, iright
    real(r64) :: dOmega, w

    ncont = size(Omegas_cont)
    nsamp = size(Omegas_samp)
    
    allocate(Imeps_T(nsamp))

    dOmega = Omegas_cont(2) - Omegas_cont(1)
    
    do isamp = 1, nsamp
       w = Omegas_samp(isamp)
       ileft = minloc(abs(Omegas_cont - w), dim = 1)
       if(ileft == ncont) then
          Imeps_T(isamp) = spec_eps_T(ileft)
       else
          iright = ileft + 1
          Imeps_T(isamp) = ((Omegas_cont(iright) - w)*spec_eps_T(ileft) + &
               (w - Omegas_cont(ileft))*spec_eps_T(iright))/dOmega
       end if
    end do

    Imeps_T = -pi*Imeps_T
  end subroutine head_polarizability_imag_3d_T
  
  subroutine spectral_head_polarizability_3d_qpath(spec_eps, Omegas, qcrys, &
       el, wann, crys, tetrahedra)
    !! Spectral head of the bare polarizability of the 3d Kohn-Sham system using
    !! Eq. 16 of Shishkin and Kresse Phys. Rev. B 74, 035101 (2006).
    !!
    !! Here we calculate the diagonal in G-G' space. Moreover,
    !! we use the approximation G.r -> 0.
    !!
    !! spec_eps Spectral head of the bare polarizability
    !! Omega Energy of excitation in the electron gas
    !! qcrys Transfer wave vector in fractionl coordinates
    !! el Electron data type
    !! wann Wannier data type
    !! crys Crystal data type
    !! tetrahedra Delta evaulators selector

    real(r64), intent(in) :: Omegas(:), qcrys(3)
    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    logical, intent(in) :: tetrahedra
    real(r64), allocatable, intent(out) :: spec_eps(:)
    
    !Locals
    integer(i64) :: m, n, ik, iOmega, nOmegas, k_indvec(3), kp_indvec(3)
    real(r64) :: overlap, ek, ekp, el_ens_kp(1, el%numbands), kppathvecs(1, 3)
    complex(r64) :: el_evecs_kp(1, el%numbands, el%numbands)
    procedure(delta_fn), pointer :: delta_fn_ptr => null()

    real(r64) :: dOmega

    nOmegas = size(Omegas)

    dOmega = Omegas(2) - Omegas(1)
    
    allocate(spec_eps(nOmegas))

    !Associate delta function procedure pointer
    delta_fn_ptr => get_delta_fn_pointer(tetrahedra)

    spec_eps = 0.0
    !Below, we will sum over k, m, and n
    do ik = 1, el%nwv
       kppathvecs(1, :) = el%wavevecs(ik, :) + qcrys
       call wann%el_wann(crys = crys, &
            nk = 1_i64, &
            kvecs = kppathvecs, &
            energies = el_ens_kp, &
            evecs = el_evecs_kp, &
            scissor = el%scissor)

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

  !DEBUG/TEST
  subroutine calculate_RPA_dielectric_3d_G0_scratch(el, crys, num, wann)
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
    real(r64) :: qcrys(3)
    integer(i64) :: iq, iOmega, numomega, numq, &
         start, end, chunk, num_active_images, qxmesh
    real(r64), allocatable :: spec_X0(:), ImX0(:), ReX0(:)
    complex(r64), allocatable :: diel(:, :)
    character(len = 1024) :: filename
    real(r64) :: omega_plasma

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
       qlist(iq, :) = [(iq - 1.0_r64)/qxmesh, (iq - 1.0_r64)/qxmesh, 0.0_r64]
       qmaglist(iq) = twonorm(matmul(crys%reclattvecs, qlist(iq, :)))
    end do
    call sort(qmaglist)

    !Create energy grid
    numomega = 601 !6 !5 !1001 !1001
    allocate(energylist(numomega))
    call linspace(energylist, -0.5_r64, 0.5_r64, numomega)
    
    !Allocate diel_ik to hold maximum possible Omega
    allocate(diel(numq, numomega))
    diel = 0.0_r64

    !Allocate spectral and imaginary X0 (defined on uniform frequency mesh)
    allocate(spec_X0(numomega), ImX0(numomega))
    
    !Distribute points among images
    call distribute_points(numq, chunk, start, end, num_active_images)

    if(this_image() == 1) then
       write(*, "(A, I10)") " #q-vecs = ", numq
       write(*, "(A, I10)") " #q-vecs/image <= ", chunk
    end if

    do iq = start, end !Over IBZ k points
       qcrys = qlist(iq, :) !crystal coordinates

       call spectral_head_polarizability_3d_qpath(&
            spec_X0, energylist, qcrys, el, wann, crys, num%tetrahedra)
       
       ImX0 = -pi*spec_X0
       call hilbert_transform(-ImX0, ReX0)
          
       !Calculate RPA dielectric (diagonal in G-G' space)
!!$       if(all(qcrys == 0.0_r64)) then
!!$          diel(iq, :) = 1.0_r64 - &
!!$               (omega_plasma/energylist(:))**2
!!$       else
       diel(iq, :) = crys%epsiloninf - &
            1.0_r64/qmaglist(iq)**2* &
            (ReX0 + oneI*ImX0)/perm0*qe*1.0e9_r64
!!$       end if
    end do

    call co_sum(diel)
    
    !Print to file
    call write2file_rank2_real("RPA_dielectric_3D_G0_qpath", qlist)
    call write2file_rank1_real("RPA_dielectric_3D_G0_qmagpath", qmaglist)
    call write2file_rank1_real("RPA_dielectric_3D_G0_Omega", energylist)
    call write2file_rank2_real("RPA_dielectric_3D_G0_real", real(diel))
    call write2file_rank2_real("RPA_dielectric_3D_G0_imag", imag(diel))
  end subroutine calculate_RPA_dielectric_3d_G0_scratch
end module screening_module

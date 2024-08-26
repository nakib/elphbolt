program screening_comparison
  !! Test for the screening stuff

  use precision, only: r64, i64
  use params, only: hbar, hbar_eVps, me, twopi, pi, kB, qe, bohr2nm, perm0
  use misc, only: qdist, linspace, compsimps, outer, sort, &
       write2file_rank2_real, write2file_rank1_real, twonorm, exit_with_message

  implicit none

  !integer :: itest
  !integer, parameter :: num_tests = 1
  !type(testify) :: test_array(num_tests), tests_all

  integer(i64) :: ik, numomega, numq
  real(r64) :: mu, eF, kF, kTF, beta, en_plasmon
  real(r64), allocatable :: el_ens_parabolic(:), qmags(:)

  !concentration and temperature
  real(r64), parameter :: conc = 1.0e18 !cm^-3
  real(r64), parameter :: T = 1.0

  !real(r64), parameter :: m_eff = 0.267*me !Si
  !real(r64), parameter :: m_eff = 0.2*me !wGaN

  !GaAs
  real(r64), parameter :: m_eff = 0.07*me
  real(r64), parameter :: epsiloninf = 11.1
  !real(r64), parameter :: epsilon0 = 12.9
  
  real(r64), allocatable :: imeps(:, :), reeps(:, :), Omegas(:)
  
  if(this_image() == 1) then
     !write(*, '(A)')  'Screening test for wGaN'
     write(*, '(A)')  'Screening test for GaAs'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if 
  
  !Test counter
  !itest = 0

  !Set inverse temperature energy
  beta = 1.0_r64/kB/T

  !Calculate Fermi wave vector for model band (degenerate limit)
  kF = (3.0_r64*pi**2*conc)**(1.0_r64/3.0_r64)*1.0e-7_r64 !nm^-1
  print*, 'Fermi wave vector = ', kF, ' nm^-1'

  !Calculate Fermi energy for model band (degenerate limit)
  eF =  energy_parabolic(kF, m_eff)
  print*, 'Fermi energy = ', eF, ' eV'
  
  !Calculate chemical potential for model band to match carrier conc.
  mu = chempot(conc, m_eff, beta, eF)
  
  !Calculate Plasmon energy
  !en_plasmon = 1.0e-9_r64*hbar_evps*qe*sqrt(conc/perm0/epsilon0/m_eff) !eV

  en_plasmon = 1.0e-9_r64*hbar_evps*qe*sqrt(conc/perm0/epsiloninf/m_eff) !eV
  print*, 'Plasmon energy = ', en_plasmon, ' eV'

  !Calculate Thomas-Fermi screening wave vector (T << T_F limit)
  !TODO Replace it with the more general expression involving
  !energy integral over f0(1 - f0).
  kTF = 1.0/kF/bohr2nm/pi**2
  print*, 'Thomas-Fermi screening wave vector = ', kTF, ' nm^-1'

  !Create wave vector mesh
  numq = 400
  call linspace(qmags, 0.0_r64, 3.0_r64*kF, numq)
  call write2file_rank1_real("RPA_test_qmags", qmags)
  
  !Create bosonic energy mesh
  numomega = 400
  call linspace(Omegas, 0.0_r64, 3.0_r64*eF, numomega)
  call write2file_rank1_real("RPA_test_Omegas", Omegas)
  
  !Calculate analytic Im RPA dielectric function
  call calculate_Imeps(qmags, Omegas, mu, m_eff, eF, kF, beta, Imeps)
  call write2file_rank2_real("model_RPA_dielectric_3D_imag", Imeps)
  
  !Calculate analytic Re RPA dielectric function
  call calculate_Reeps(qmags, Omegas, mu, m_eff, eF, en_plasmon, &
       kF, kTF, epsiloninf, beta, Reeps)
  call write2file_rank2_real("model_RPA_dielectric_3D_real", Reeps)
  
contains
  
  pure real(r64) elemental function energy_parabolic(k, m_eff)
    !! Parabolic band energy for a given wave vector magnitude
    !!
    !! k Wave vector magnitude in nm^-1
    !! m_eff Effective mass in Kg

    real(r64), intent(in) :: k, m_eff
    
    energy_parabolic = 0.5_r64*(hbar*k)**2/m_eff*1.0e-6_r64/qe !eV
  end function energy_parabolic
  
  real(r64) function chempot(conc, m_eff, beta, eF)
    !!Use bisection method to find chemical potential
    !!for a given carrier concentration
    !!
    !! conc Carrier concentration in cm^-3
    !! m_eff Effective mass in Kg
    !! beta Inverse temperature energy in eV^-1

    real(r64), intent(in) :: conc
    real(r64), intent(in) :: m_eff
    real(r64), intent(in) :: beta, eF

    integer :: i
    integer, parameter :: maxiter = 100
    real(r64) :: a, b, tmp, aux, thresh, upper

    upper = eF*beta + 15.0 !to be used as inifinity of Fermi integral
    
    a = -5.0 !eV, "Safe" lower bound
    b = 5.0 !eV, "Safe" upper bound

    thresh = 1.0e-12_r64
    
    tmp = 2.0_r64*(0.5_r64*m_eff/hbar_eVps**2/beta/pi/qe)**1.5_r64*1.0e30_r64 !cm^-3

    do i = 1, maxiter
       chempot = 0.5*(a + b)

       aux = tmp*fdi(0.5_r64, chempot*beta, upper)
       
       if(abs(aux - conc)/conc < thresh) then
          exit
       else if(aux < conc) then
          a = chempot
       else
          b = chempot
       end if
    end do
    write(*, "(A, 1E16.8, A)") 'requested ne = ', conc, ' cm^-3'
    write(*, "(A, 1E16.8, A)") 'calculated ne in parabolic model = ', aux, ' cm^-3'
    write(*, "(A, 1E16.8, A)") 'calculated chem. pot. in parabolic model = ', chempot, ' eV'
  end function chempot

  real(r64) function fdi(j, eta, upper)
    !! Fermi-Dirac integral for positive j
    
    real(r64), intent(in) :: j, eta, upper

    integer(i64), parameter :: ngrid = 100000_i64
    real(r64), allocatable :: x(:)
    real(r64) :: dx
    
    if(j < 0.0_r64) call exit_with_message("Negative j is not allowed. Exiting.")
    
    call linspace(x, 0.0_r64, upper, ngrid)

    dx = x(2) - x(1)
    call compsimps(x**j/(exp(x - eta) + 1.0_r64), dx, fdi)
    fdi = fdi/gamma(j + 1.0_r64)
  end function fdi

  real(r64) function fdi_minus1half(eta)
    !! Fermi-Dirac integral j = -1/2
    !! Source: Frank G. Lether
    !! Journal of Scientific Computing, Vol. 15, No. 4, 2000
    
    real(r64), intent(in) :: eta

    integer :: k
    real(r64) :: last_term, aux

    last_term = 0.0_r64
    do k = 1, 20
       aux = sqrt(eta**2 + (2.0_r64*k - 1.0_r64)**2*pi**2)
       last_term = last_term + sqrt(aux - eta)/aux
    end do
    
    fdi_minus1half = 8220.0_r64/919 + 3923.0_r64/110242*eta + &
         27.0_r64/381503*eta**2 - eta**3/3553038.0_r64 - &
         eta**4/714900621.0_r64 - eta**5/128458636383.0_r64 - &
         sqrt(twopi)*last_term
  end function fdi_minus1half

  subroutine calculate_Imeps(qmags, ens, chempot, m_eff, eF, kF, beta, Imeps)
    !! Imaginary part of RPA dielectric for the isotropic band model.
    !!
    !! qmags Magnitude of probe wave vectors in nm^-1
    !! ens Probe energies in eV
    !! chempot Chemical potential in eV
    !! m_eff Effective mass of model band in Kg
    !! eF Fermi energy (degenerate limit => 0K) in eV
    !! kF Fermi wave vector (degenerate limit => 0K) in nm^-1
    !! beta Inverse temperature energy in eV^-1
    !! Imeps Imaginary part of RPA dielectric for the isotropic band model

    real(r64), intent(in) :: qmags(:), ens(:), chempot, m_eff, eF, kF, beta
    real(r64), allocatable :: Imeps(:, :)

    !Locals
    integer :: iOmega
    real(r64) :: u(size(qmags), size(ens))
    
    allocate(Imeps(size(qmags), size(ens)))
    
    call outer(0.5_r64*kF/qmags, ens/eF, u)

    do iOmega = 1, size(ens)
       Imeps(:, iOmega) = &
            real(log(&
            (1.0_r64 + exp(beta*(chempot - eF*(u(:, iOmega) - qmags(:)/kF/2.0_r64)**2)))/ &
            (1.0_r64 + exp(beta*(chempot - eF*(u(:, iOmega) + qmags(:)/kF/2.0_r64)**2)))))/ &
            (qmags(:)/kF)**3
    end do
    Imeps(1, :) = 0.0_r64
    Imeps = (m_eff/me/bohr2nm/kF/eF/beta)*Imeps

!!$    !This gives the same result as the expression above:
!!$    !Locals
!!$    integer :: iOmega, iq
!!$    real(r64) :: E1(size(qmags), size(ens)), E2(size(qmags), size(ens)), Eq(size(qmags))
!!$
!!$    allocate(Imeps(size(qmags), size(ens)))
!!$
!!$    Eq = energy_parabolic(qmags, m_eff)
!!$
!!$    do iq = 1, size(qmags)
!!$       E1(iq, :) = (ens(:) - Eq(iq))**2/4.0/Eq(iq)
!!$       E2(iq, :) = (ens(:) + Eq(iq))**2/4.0/Eq(iq)
!!$    end do    
!!$
!!$    do iOmega = 1, size(ens)
!!$       Imeps(:, iOmega) = &
!!$            log(&
!!$            (1.0_r64 + exp(beta*(chempot - E1(:, iOmega))))/ &
!!$            (1.0_r64 + exp(beta*(chempot - E2(:, iOmega)))))/ &
!!$            (qmags(:))**3
!!$    end do
!!$    !Imeps(1, :) = 0.0_r64
!!$    Imeps = 2.0*m_eff**2/(me*bohr2nm*beta*hbar**2)*1.0e6*qe*Imeps
  end subroutine calculate_Imeps

  subroutine calculate_Reeps(qmags, ens, chempot, m_eff, eF, eplasmon, &
       kF, kTF, epsinf, beta, Reeps)
    !! Real part of RPA dielectric for the isotropic band model.
    !!
    !! qmags Magnitude of probe wave vectors in nm^-1
    !! ens Probe energies in eV
    !! chempot Chemical potential in eV
    !! m_eff Effective mass of model band in Kg
    !! eF Fermi energy in eV
    !! eplasmon Plasmon energy in eV
    !! kF Fermi wave vector (degenerate limit => 0K) in nm^-1
    !! kTF Thomas-Fermi wave vector in nm^-1
    !! beta Inverse temperature energy in eV^-1
    !! Reeps Real part of RPA dielectric for the isotropic band model.

    real(r64), intent(in) :: qmags(:), ens(:), chempot, m_eff, eF, kF, &
         kTF, epsinf, beta, eplasmon
    real(r64), allocatable :: Reeps(:, :)
    
    !Locals
    integer(i64) :: iOmega, iq, ngrid
    real(r64) :: ymax, dy, aux0, aux1, aux2, D, x, eta, ks_squared
    real(r64) :: u(size(qmags), size(ens)), z(size(qmags))
    real(r64), allocatable :: y(:), I0(:)
    
    !Here I use pra 29 1471.
    !Note that in the paper above the Bohr radius is renormalized.
    !Here we need an extra factor of ms/me.

    !Magic numbers?
    ngrid = 1000
    ymax = 10.0!20.0_r64

    allocate(y(ngrid), I0(ngrid), Reeps(size(qmags), size(ens)))
    
    call linspace(y, 0.0_r64, ymax, ngrid)

    dy = y(2) - y(1)
    D = ef*beta
    eta = chempot*beta
    I0 = y/(exp(D*y**2 - eta) + 1.0_r64)

    z = 0.5_r64*qmags/kF

    call outer(0.5_r64*kF/qmags, ens/eF, u)
        
    !Calculate screening wave vector (squared)
    !ks_squared = 0.5_r64*kTF**2*sqrt(1.0_r64/D)*fdi_minus1half(eta)
    !print*, 'ks = ', sqrt(ks_squared), ' nm^-1'

    !Non-zero energy and momentum
    do iOmega = 1, size(ens)
       do iq = 2, size(qmags)
          aux0 = 1.0_r64/z(iq)**3

          x = u(iq, iOmega) + z(iq)
          call compsimps(I0*real(log(abs((x + y)/(x - y)))), dy, aux1)

          x = u(iq, iOmega) - z(iq)
          call compsimps(I0*real(log(abs((x + y)/(x - y)))), dy, aux2)

          Reeps(iq, iOmega) = aux0*(aux1 - aux2)
       end do
    end do

    !       valence band contribution + electronic contribution
    Reeps = epsinf + (0.25_r64/pi/kF/bohr2nm*m_eff/me)*Reeps

    !Omega -> 0 limit
    !Reeps(2:size(qmags), 1) = epsinf*&
    !     (1.0_r64 + ks_squared/qmags(2:size(qmags))**2)

    !q -> limit
    Reeps(1, :) = 1.0_r64 - (eplasmon/ens)**2

!!$    !Check where the plasmon mode is at the Gamma point
!!$    do iOmega = 1, size(ens)
!!$       if(abs(Reeps(1, iOmega)) <= 1.0e-1) print*, ens(iOmega), ' eV'
!!$    end do
  end subroutine calculate_Reeps

end program screening_comparison

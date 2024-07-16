program screening_comparison
  !! Test for the screening stuff

  use precision, only: r64, i64
  use params, only: hbar, hbar_eVps, me, twopi, pi, kB, qe, bohr2nm
  use misc, only: qdist, linspace, compsimps, outer, sort, &
       write2file_rank2_real, write2file_rank1_real, twonorm
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use wannier_module, only: wannier

  implicit none

  !integer :: itest
  !integer, parameter :: num_tests = 1
  !type(testify) :: test_array(num_tests), tests_all

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(wannier) :: wann
  type(electron) :: el
  !type(timer) :: t_all, t_event

  integer(i64) :: ik, numomega, numq
  real(r64) :: mu, eF, kF, beta
  real(r64), allocatable :: el_ens_parabolic(:), qmags(:)
  real(r64), parameter :: m_eff = 0.267*me
  real(r64), allocatable :: imeps(:, :), reeps(:, :), Omegas(:)
  
  if(this_image() == 1) then
     write(*, '(A)')  'Screening test'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if 
  
  !Test counter
  !itest = 0

  !call t_all%start_timer('elphbolt: BTE')

  !call t_event%start_timer('Initialization')

  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Read EPW Wannier data
  call wann%read(num)

  !Calculate electrons
  call el%initialize(wann, crys, sym, num)

  !Set inverse temperature energy
  beta = 1.0_r64/kB/crys%T
  
  !Create grid of probe |q|
!!$  allocate(qmags(el%nwv_irred))
!!$  do ik = 1, el%nwv_irred
!!$     qmags(ik) = qdist(el%wavevecs_irred(ik, :), crys%reclattvecs)
!!$  end do
!!$  call sort(qmags)
  numq = 200
  call linspace(qmags, 0.0_r64, twopi/twonorm(crys%lattvecs(:, 1)), numq)
  call write2file_rank1_real("RPA_test_qmags", qmags)

  !Calculate parabolic dispersion 
  !allocate(el_ens_parabolic(el%nwv_irred))
  !allocate(el_ens_parabolic(numq))  
  !el_ens_parabolic = energy_parabolic(qmags, m_eff)
  !call write2file_rank1_real("model_el_ens_parabolic", el_ens_parabolic)
  
  !print*, qmags(1:10)
  !print*, el_ens_parabolic(1:10)

  !Calculate chemical potential for model band to match carrier conc.
  mu = chempot(el%conc_el, m_eff, beta)
  
  !Calculate Fermi wave vector for model band (degenerate limit)
  kF = (3.0_r64*pi**2*el%conc_el)**(1.0_r64/3.0_r64)*1.0e-7_r64 !nm^-1
  print*, 'Fermi wave vector = ', kF, ' nm^-1'

  !Calculate Fermi energy for model band (degenerate limit)
  eF =  energy_parabolic(kF, m_eff)
  print*, 'Fermi energy = ', eF, ' eV'

  !Create bosonic energy mesh
  numomega = 100
  call linspace(Omegas, 0.0_r64, 10.0_r64*eF, numomega)
  call write2file_rank1_real("RPA_test_Omegas", Omegas)
  
  !Calculate analytic Im RPA dielectric function
  call calculate_Imeps(qmags, Omegas, mu, m_eff, eF, kF, beta, Imeps)
  call write2file_rank2_real("model_RPA_dielectric_3D_imag", Imeps)
  
  !Calculate analytic Re RPA dielectric function
  call calculate_Reeps(qmags, Omegas, mu, m_eff, eF, kF, beta, Reeps)
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
  
  real(r64) function chempot(conc, m_eff, beta)
    !!Use bisection method to find chemical potential
    !!for a given carrier concentration
    !!
    !! conc Carrier concentration in cm^-3
    !! m_eff Effective mass in Kg
    !! beta Inverse temperature energy in eV^-1

    real(r64), intent(in) :: conc
    real(r64), intent(in) :: m_eff
    real(r64), intent(in) :: beta

    integer :: i
    integer, parameter :: maxiter = 100
    real(r64) :: a, b, tmp, aux, thresh

    a = -5.0 !eV, "Safe" lower bound
    b = 5.0 !eV, "Safe" upper bound

    thresh = 1.0e-12_r64
    
    tmp = 2.0_r64*(0.5_r64*m_eff/hbar_eVps**2/beta/pi/qe)**1.5_r64*1.0e30_r64 !cm^-3

    do i = 1, maxiter
       chempot = 0.5*(a + b)

       aux = tmp*fdi(0.5_r64, chempot*beta)
       
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

  real(r64) function fdi(j, x)
    !! Fermi-Dirac integral
    
    real(r64), intent(in) :: j, x

    integer(i64), parameter :: ngrid = 10000_i64
    real(r64), allocatable :: t(:)
    real(r64) :: dt

    !Here 10 is infinity...
    call linspace(t, 0.0_r64, 10.0_r64, ngrid)

    dt = t(2) - t(1)
    call compsimps(t**j/(exp(t - x) + 1.0_r64), dt, fdi)
    fdi = fdi/gamma(j + 1.0_r64)    
  end function fdi

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
    
    call outer(0.5_r64/qmags/kF, ens/eF, u)

    do iOmega = 1, size(ens)
       Imeps(:, iOmega) = &
            real(log(&
            (1.0_r64 + exp(beta*(chempot - eF*(u(:, iOmega) - qmags(:)/kF/2.0_r64)**2)))/ &
            (1.0_r64 + exp(beta*(chempot - eF*(u(:, iOmega) + qmags(:)/kF/2.0_r64)**2)))))/ &
            (qmags(:)/kF)**3
    end do
    Imeps(1, :) = 0.0_r64
    Imeps = (m_eff/me/bohr2nm/kF/eF/beta)*Imeps
  end subroutine calculate_Imeps

  subroutine calculate_Reeps(qmags, ens, chempot, m_eff, eF, kF, beta, Reeps)
    !! Real part of RPA dielectric for the isotropic band model.
    !!
    !! qmags Magnitude of probe wave vectors in nm^-1
    !! ens Probe energies in eV
    !! chempot Chemical potential in eV
    !! m_eff Effective mass of model band in Kg
    !! kF Fermi wave vector (degenerate limit => 0K) in nm^-1
    !! beta Inverse temperature energy in eV^-1
    !! Reeps Real part of RPA dielectric for the isotropic band model.

    real(r64), intent(in) :: qmags(:), ens(:), chempot, m_eff, eF, kF, beta
    real(r64), allocatable :: Reeps(:, :)
    
    !Locals
    integer(i64) :: iOmega, iq, ngrid
    real(r64) :: ymax, dy, aux0, aux1, aux2, D, x, eta
    real(r64) :: u(size(qmags), size(ens)), z(size(qmags))
    real(r64), allocatable :: y(:), I0(:)
    
    !Here I use pra 29 1471.
    !Note that in the paper above the Bohr radius is renormalized.
    !Here we need an extra factor of ms/me.

    !Magic numbers?
    ngrid = 200
    ymax = 10.0_r64

    allocate(y(ngrid), I0(ngrid), Reeps(size(qmags), size(ens)))
    
    call linspace(y, 0.0_r64, ymax, ngrid)

    dy = y(2) - y(1)
    D = ef*beta
    eta = chempot*beta
    I0 = y/(exp(D*y**2 - eta) + 1.0_r64)

    z = 0.5_r64*qmags/kF

    call outer(0.5_r64/qmags/kF, ens/eF, u)
    
    Reeps = 0.0_r64
    do iOmega = 2, size(ens)
       do iq = 2, size(qmags)
          aux0 = 1.0_r64/z(iq)**3

          x = u(iq, iOmega) + z(iq)
          call compsimps(I0*real(log(abs((x + y)/(x - y)))), dy, aux1)

          x = u(iq, iOmega) - z(iq)
          call compsimps(I0*real(log(abs((x + y)/(x - y)))), dy, aux2)

          Reeps(iq, iOmega) = aux0*(aux1 - aux2)
       end do
    end do
    
    Reeps = 1.0_r64 + (0.25_r64/pi/kF/bohr2nm*m_eff/me)*Reeps
  end subroutine calculate_Reeps

end program screening_comparison

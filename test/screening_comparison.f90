program screening_comparison
  !! Test for the screening stuff

  use precision, only: r64, i64
  use params, only: hbar, hbar_eVps, me, pi, kB, qe
  use misc, only: qdist, linspace, compsimps
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

  integer :: ik
  real(r64) :: mu
  real(r64), allocatable :: el_ens_parabolic(:), kmags(:)
  real(r64), parameter :: m_eff = 0.267*me
  
  if(this_image() == 1) then
     write(*, '(A)')  'Screening test'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if 

  !Set effective mass of model band
  !m_eff = 0.267*me
  
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

  !Create grid of |k|
  allocate(kmags(el%nwv_irred))
  do ik = 1, el%nwv_irred
     kmags(ik) = qdist(el%wavevecs_irred(ik, :), crys%reclattvecs)
  end do

  !Calculate parabolic dispersion 
  allocate(el_ens_parabolic(el%nwv_irred))  
  el_ens_parabolic = energy_parabolic(kmags, m_eff)
  print*, kmags(1:10)
  print*, el_ens_parabolic(1:10)

  !Calculate chemical potential for model band to match carrier conc.
  mu = chempot(el%conc_el, m_eff, 1.0_r64/kB/crys%T)
  print*, 'Chemical potential = ', mu, ' eV'

  call exit
  !Calculate Fermi wave vector for model band
  !TODO
  
  !TODO Calculate analytic RPA dielectric function
  
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
    integer, parameter :: maxiter = 1000
    real(r64) :: a, b, tmp, aux

    a = -5.0 !eV, "Safe" lower bound
    b = 5.0 !eV, "Safe" upper bound

    !TODO Check units
    tmp = 2.0_r64*(0.5_r64*m_eff/hbar_eVps**2/beta/pi/qe)**1.5_r64*1.0e30_r64 !cm^-3
    print*, 'tmp = ', tmp

    do i = 1, maxiter
       chempot = 0.5*(a + b)
       
       aux = tmp*fdi(0.5_r64, chempot*beta)
       if(abs(aux - conc)/conc < 1e-12_r64) then
          exit
       else if(aux < conc) then
          a = chempot
       else
          b = chempot
       end if
    end do
    !print'(A,E16.8,A)', 'calculated ne =', aux*1d-6, ' 1/cm^3'
    !print'(A,E16.8,A)', 'calculated chempot =', chempot, ' J'
    print*,'calculated ne = ', aux, ' cm^-3'
    print*,'calculated chem. pot. = ', chempot, ' eV'
  end function chempot

  !TODO Remove all "magic" parameters in the function below
  real(r64) function fdi(j, x)
    !! Fermi-Dirac integral
    
    real(r64), intent(in) :: j, x

    integer(i64), parameter :: ngrid = 1000_i64
    real(r64), allocatable :: t(:)
    real(r64) :: dt

    call linspace(t, 0.0_r64, 5.0_r64, ngrid)

    dt = t(2) - t(1)
    call compsimps(t**j/(exp(t - x) + 1.0_r64), dt, fdi)
    fdi = fdi/gamma(j + 1.0_r64)    
  end function fdi
!!$  subroutine calculate_Reeps
!!$  end subroutine calculate_Reeps
!!$
!!$  subroutine calculate_Imeps
!!$  end subroutine calculate_Imeps
end program screening_comparison

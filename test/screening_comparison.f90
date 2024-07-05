program screening_comparison
  !! Test for the screening stuff

  use precision, only: r64, i64
  use params, only: hbar_eVps, me
  use misc, only: qdist
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

  real(r64), allocatable :: el_ens_parabolic(:), kmags(:)
  real(r64) :: m_eff
  
  if(this_image() == 1) then
     write(*, '(A)')  'Screening test'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if 

  !Set effective mass of model band
  m_eff = 0.267*me
  
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
  allocate(kmags(el%nwv_irred)
  do ik = 1, el%nwv_irred
     kmags(ik) = qdist(el%wavevecs_irred(ik))
  end do

  !Calculate parabolic dispersion 
  allocate(el_ens_parabolic(el%nwv_irred))  
  el_ens_parabolic = energy_parabolic(kmags, m_eff)

  !Calculate chemical potential for model band to match carrier conc.
  !TODO

  !Calculate Fermi wave vector for model band
  !TODO
  
  !TODO Calculate analytic RPA dielectric function
  
contains
  
  pure elemental function energy_parabolic(k, m_eff)
    !! Parabolic band energy for a given wave vector magnitude
    !!
    !! k Wave vector magnitude in nm^-1
    !! m_eff Effective mass in Kg

    real(r64), intent(in) :: k, m_eff
    
    energy_parabolic = 0.5_r64*(hbar_eVps*k)**2/m_eff
  end function energy_parabolic

!!$  subroutine calculate_Reeps
!!$  end subroutine calculate_Reeps
!!$
!!$  subroutine calculate_Imeps
!!$  end subroutine calculate_Imeps
end program screening_comparison

program test_autodiff

  use params, only : r64, i64
  use misc, only : linspace
  use testify_m, only : testify
  use autodiff_m, only : autodiff, &
       operator(+), operator(-), operator(*), operator(/)

  implicit none
  
  integer :: itest
  integer, parameter :: num_tests = 2
  type(testify) :: test_array(num_tests), tests_all

  integer :: i
  integer(i64), parameter :: N = 3
  real(r64), allocatable :: x(:)
  type(autodiff) :: ad_x(N), ad_fx(N)
    
  call linspace(x, 0.0_r64, 1.0_r64, N)

  itest = 1
  test_array(itest) = testify("f(x) = x^2 + 2x")
  do i = 1, N
     ad_x(i)%f = x(i)
     ad_x(i)%df = 1.0_r64 !take derivative w.r.t x
     ad_fx(i) = ad_x(i)*ad_x(i) + 2.0_r64*ad_x(i)
  end do
  call test_array(itest)%assert(&
       [(ad_fx(i)%df, i = 1, N)], &
       [2.0_r64, 3.0_r64, 4.0_r64], &
       tol = 1.0e-12_r64)

  itest = 2
  test_array(itest) = testify("f(x) = 1/(1 + x^2) - x/2")
  do i = 1, N
     ad_x(i)%f = x(i)
     ad_x(i)%df = 1.0_r64 !take derivative w.r.t x
     ad_fx(i) = 1.0_r64/(1.0_r64 + ad_x(i)*ad_x(i)) &
          - ad_x(i)/2.0_r64

     !print*, ad_x(i)%f, ad_fx(i)%f, ad_fx(i)%df
  end do
  call test_array(itest)%assert(&
       [(ad_fx(i)%df, i = 1, N)], &
       [-0.5_r64, -16.0_r64/25.0_r64 - 0.5_r64, -1.0_r64], &
       tol = 1.0e-12_r64)
  
  tests_all = testify(test_array)
  call tests_all%report

  if(tests_all%get_status() .eqv. .false.) error stop -1
end program test_autodiff

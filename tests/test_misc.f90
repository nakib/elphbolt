program test_misc

  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use misc, only: int_div
  
  implicit none

  integer :: itest
  integer, parameter :: num_tests = 3
  type(testify) :: test_array(num_tests), tests_all
  integer(i64) :: quotient, remainder
  
  test_array(1) = testify("5/2")
  call int_div(5_i64, 2_i64, quotient, remainder)
  call test_array(1)%assert([quotient, remainder], [2_i64, 1_i64])

  test_array(2) = testify("9/3")
  call int_div(9_i64, 3_i64, quotient, remainder)
  call test_array(2)%assert([quotient, remainder], [3_i64, 0_i64])

  test_array(3) = testify("3/10")
  call int_div(3_i64, 10_i64, quotient, remainder)
  call test_array(3)%assert([quotient, remainder], [0_i64, 3_i64])

  tests_all = testify(test_array)
  call tests_all%report
  
  if (tests_all%get_status() .eqv. .false.) error stop
end program test_misc

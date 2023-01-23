program test_misc

  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use params, only: pi
  use misc, only: int_div, expi, trace, kronecker, sort, cross_product, &
       twonorm, binsearch, mux_vector, demux_vector, interpolate
  
  implicit none

  integer :: itest
  integer, parameter :: num_tests = 13
  type(testify) :: test_array(num_tests), tests_all
  integer(i64) :: index, quotient, remainder, int_array(5), v1(3), v2(3), v1_muxed, v2_muxed
  real(r64) :: pauli1(2, 2), ipauli2(2, 2), pauli3(2, 2), &
       real_array(5)

  !Some data to be used in the tests below
  pauli1 = reshape([0.0_r64, 1.0_r64, 1.0_r64, 0.0_r64], [2, 2])
  ipauli2 = reshape([0.0_r64, -1.0_r64, 1.0_r64, 0.0_r64], [2, 2])
  pauli3 = reshape([1.0_r64, 0.0_r64, 0.0_r64, -1.0_r64], [2, 2])
 
  !int_div
  itest = 1
  test_array(itest) = testify("int_div 5/2")
  call int_div(5_i64, 2_i64, quotient, remainder)
  call test_array(itest)%assert([quotient, remainder], [2_i64, 1_i64])

  itest = itest + 1
  test_array(itest) = testify("int_div 9/3")
  call int_div(9_i64, 3_i64, quotient, remainder)
  call test_array(itest)%assert([quotient, remainder], [3_i64, 0_i64])

  itest = itest + 1
  test_array(itest) = testify("int_div 3/10")
  call int_div(3_i64, 10_i64, quotient, remainder)
  call test_array(itest)%assert([quotient, remainder], [0_i64, 3_i64])

  !distribute_points
  !TODO This is a coarray dependent test. Will revisit.

  !cross_product
  itest = itest + 1
  test_array(itest) = testify("cross_product i x j, j x k, i x k")
  call test_array(itest)%assert(&
       [cross_product([1.0_r64, 0.0_r64, 0.0_r64], [0.0_r64, 1.0_r64, 0.0_r64]), &
        cross_product([0.0_r64, 1.0_r64, 0.0_r64], [0.0_r64, 0.0_r64, 1.0_r64]), &
        cross_product([1.0_r64, 0.0_r64, 0.0_r64], [0.0_r64, 0.0_r64, 1.0_r64])], &
       [[0.0_r64, 0.0_r64, 1.0_r64], &
        [1.0_r64, 0.0_r64, 0.0_r64], &
        [0.0_r64, -1.0_r64, 0.0_r64]])

  !kronecker
  itest = itest + 1
  test_array(itest) = testify("kronecker")
  call test_array(itest)%assert(&
       [kronecker(0_i64, 1_i64), kronecker(-1_i64, -1_i64)], &
       [0_i64, 1_i64]) 
  
  !expi
  itest = itest + 1
  test_array(itest) = testify("expi 0, pi")
  call test_array(itest)%assert(&
       [expi(0.0_r64), expi(pi)], &
       [(1.0_r64, 0.0_r64), cmplx(cos(pi), sin(pi), r64)])
  
  !twonorm_real_rank1, *_rank2
  itest = itest + 1
  test_array(itest) = testify("twonorm rank-1, rank-2")
  call test_array(itest)%assert(&
       [twonorm([sqrt(1.0_r64/3.0_r64), -sqrt(1.0_r64/3.0_r64), sqrt(1.0_r64/3.0_r64)]), twonorm(pauli1)], &
       [1.0_r64, sqrt(2.0_r64)])
  
  !trace
  itest = itest + 1
  test_array(itest) = testify("trace pauli1, ipauli2, pauli3, -ipauli1..3")
  call test_array(itest)%assert(&
       [trace(pauli1), trace(ipauli2), trace(pauli3), trace(-matmul(pauli1, matmul(ipauli2, pauli3)))], &
       [0.0_r64, 0.0_r64, 0.0_r64, 2.0_r64])
  
  !sort_int
  itest = itest + 1
  test_array(itest) = testify("sort_int")
  int_array = [105_i64, 976_i64, 276_i64, -7865_i64, 0_i64]
  call sort(int_array)
  call test_array(itest)%assert(&
       int_array, &
       [-7865_i64, 0_i64, 105_i64, 276_i64, 976_i64])

  !sort_real
  itest = itest + 1
  test_array(itest) = testify("sort_real")
  real_array = [105.0_r64, 976.0_r64, 276.0_r64, -7865.0_r64, 0.0_r64]
  call sort(real_array)
  call test_array(itest)%assert(&
       real_array, &
       [-7865.0_r64, 0.0_r64, 105.0_r64, 276.0_r64, 976.0_r64])
  
  !binsearch
  itest = itest + 1
  test_array(itest) = testify("sort_real")
  int_array = [-105_i64, 105_i64, 105_i64, 105_i64, 0_i64]
  call sort(int_array)
  call binsearch(int_array, 105_i64, index)
  call test_array(itest)%assert(index, 3_i64)

  !TODO
  
  !compsimps

  !mux_vector
  itest = itest + 1
  test_array(itest) = testify("mux_vector base 0, base 1")
  call test_array(itest)%assert(&
       [mux_vector(1_i64*[1, 2, 3], 1_i64*[4, 4, 4], 0_i64), &
        mux_vector(1_i64*[1, 2, 3], 1_i64*[4, 4, 4], 1_i64)], &
       [58_i64, 37_i64])

  !demux_vector
  itest = itest + 1
  test_array(itest) = testify("demux_vector base 0, base 1")
  v1_muxed = 58_i64 !a test muxed vector
  v2_muxed = 37_i64 !another test muxed vector
  call demux_vector(v1_muxed, v1, 1_i64*[4, 4, 4], 0_i64)
  call demux_vector(v2_muxed, v2, 1_i64*[4, 4, 4], 1_i64)
  call test_array(itest)%assert(&
       [v1, v2], 1_i64*[1, 2, 3, 1, 2, 3])
  
  !demux_mesh

  !mux_state

  !demux_state

  !Bose

  !Fermi

  !Interpolate

  !Pade_coeffs

  !Pade_continued
  
  tests_all = testify(test_array)
  call tests_all%report
  
  if (tests_all%get_status() .eqv. .false.) error stop
end program test_misc

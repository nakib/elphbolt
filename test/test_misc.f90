program test_misc
  
  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use params, only: pi, kB, oneI
  use misc, only: int_div, expi, trace, kronecker, sort, cross_product, &
       twonorm, binsearch, mux_vector, demux_vector, interpolate, coarse_grained, &
       unique, linspace, compsimps, mux_state, demux_state, demux_mesh, expm1, &
       Fermi, Bose, Pade_continued, precompute_interpolation_corners_and_weights, &
       interpolate_using_precomputed, operator(.umklapp.), shrink
  
  implicit none

  integer :: itest
  integer, parameter :: num_tests = 28
  type(testify) :: test_array(num_tests), tests_all
  integer(i64) :: index, quotient, remainder, int_array(5), v1(3), v2(3), &
       v1_muxed, v2_muxed, ik, ik1, ik2, ik3, ib1, ib2, ib3, wvmesh(3), &
       mesh_ref_array(3), nk_coarse, ninterp
  integer(i64), allocatable :: index_mesh_0(:, :), index_mesh_1(:, :), &
       ksint(:, :), idc(:, :), ik_interp(:), array_of_ints(:)
  real(r64) :: pauli1(2, 2), ipauli2(2, 2), pauli3(2, 2), &
       real_array(5), result, q1(3, 4), q2(3, 4), q3(3, 4)
  real(r64), allocatable :: integrand(:), domain(:), im_axis(:), real_func(:), &
       widc(:, :), f_coarse(:), f_interp(:), array_of_reals(:)

  print*, '<<module misc unit tests>>'
  
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
  test_array(itest) = testify("binsearch")
  int_array = [-105_i64, 105_i64, 105_i64, 105_i64, 0_i64]
  call sort(int_array)
  call binsearch(int_array, 105_i64, index)
  call test_array(itest)%assert(index, 3_i64)
  
  !compsimps
  itest = itest + 1
  test_array(itest) = testify("compsims gaussian")
  allocate(domain(300), integrand(300))
  call linspace(domain, 0.0_r64, 1.0_r64, 300_i64)
  integrand = 2.0_r64/sqrt(pi)*exp(-domain**2)
  call compsimps(integrand, domain(2) - domain(1), result)
  call test_array(itest)%assert(&
       [result], &
       [erf(1.0_r64)], tol = 1.0e-8_r64)

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
  itest = itest + 1
  test_array(itest) = testify("demux_mesh base 0, base 1")
  allocate(index_mesh_0(3, 8), index_mesh_1(3, 8))
  call demux_mesh(index_mesh_0, 1_i64*[2, 2, 2], 0_i64)
  call demux_mesh(index_mesh_1, 1_i64*[2, 2, 2], 1_i64)
  call test_array(itest)%assert(&
       [index_mesh_0, index_mesh_1], &
       1_i64*[ 0, 0, 0, &
               1, 0, 0, &
               0, 1, 0, &
               1, 1, 0, &
               0, 0, 1, &
               1, 0, 1, &
               0, 1, 1, &
               1, 1, 1, &
               1, 1, 1, &
               2, 1, 1, &
               1, 2, 1, &
               2, 2, 1, &
               1, 1, 2, &
               2, 1, 2, &
               1, 2, 2, &
               2, 2, 2  ])

  !mux_state
  itest = itest + 1
  test_array(itest) = testify("mux_state")
  call test_array(itest)%assert(&
       [mux_state(12_i64, 1_i64, 1_i64), &
        mux_state(12_i64, 12_i64, 1_i64), &
        mux_state(6_i64, 3_i64, 6_i64)], &
       1_i64*[1, 12, 33])

  !demux_state
  itest = itest + 1
  test_array(itest) = testify("demux_state")
  call demux_state(1_i64, 12_i64, ib1, ik1)
  call demux_state(12_i64, 12_i64, ib2, ik2)
  call demux_state(33_i64, 6_i64, ib3, ik3)
  call test_array(itest)%assert(&
       [ib1, ik1, ib2, ik2, ib3, ik3], &
       1_i64*[1, 1, 12, 1, 3, 6])

  !coarse_grained
  itest = itest + 1
  test_array(itest) = testify("coarse_grained")
  call test_array(itest)%assert(&
       [coarse_grained(1_i64, 1_i64*[2, 2, 2], 1_i64*[4, 2, 2]), &
        coarse_grained(2_i64, 1_i64*[2, 2, 2], 1_i64*[4, 3, 4]), &
        coarse_grained(3_i64, 1_i64*[2, 2, 2], 1_i64*[4, 3, 4]), &
        coarse_grained(4_i64, 1_i64*[2, 2, 2], 1_i64*[4, 1, 1]), &
        coarse_grained(1_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(2_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(4_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(5_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(6_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(7_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(9_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10]), &
        coarse_grained(10_i64, 1_i64*[5, 5, 5], 1_i64*[10, 10, 10])], &
        1_i64*[1, 3, 3, 1, 1, 1, 6, 6, 6, 6, 1, 1])

  !unique
  itest = itest + 1
  test_array(itest) = testify("unique")
  call test_array(itest)%assert(&
       [unique(1_i64*[4, 5, 1, 3, 3, 1, 4]), &
        unique(1_i64*[0, 0, 0]), &
        unique(1_i64*[5, 5, 4, 3, 2, 1, 1, 0, 0, -1]), &
        unique(1_i64*[1, 2, 2, 4, 4, 5, 5])], &
       1_i64*[[4, 5, 1, 3], [0], [5, 4, 3, 2, 1, 0, -1], [1, 2, 4, 5]])

  !expm1
  itest = itest + 1
  test_array(itest) = testify("expm1 (numpy precision equivalence)")
  call test_array(itest)%assert(&
       [expm1(pi), expm1(1.0e-8_r64), expm1(1.0e-10_r64)], &
       [22.140692632779267_r64, 1.0000000050000001e-08_r64, 1.0000000050000000e-08_r64], &
       tol = 1.0e-8_r64)

  !Bose
  itest = itest + 1
  test_array(itest) = testify("Bose")
  call test_array(itest)%assert(&
       [1.0e-10*Bose(1.0e-6_r64, 1.0e8_r64), &
        Bose(1.0_r64, 1.0e-2_r64)], &
       [1.0e-10*kB*1.0e8_r64/1.0e-6_r64, 0.0_r64], &
       tol = 1.0e-10_r64)
  
  !Fermi
  itest = itest + 1
  test_array(itest) = testify("Fermi")
  call test_array(itest)%assert(&
       [Fermi(1.0_r64, 1.0_r64, 300.0_r64), &
        Fermi(1.0_r64, 1.0_r64, 1.0e-2_r64)], &
       [0.5_r64, 0.5_r64], &
       tol = 1.0e-8_r64)

  !Interpolate
  itest = itest + 1
  test_array(itest) = testify("Linear interpolation")
  wvmesh = [4, 4, 4]
  mesh_ref_array = [3, 3, 3]
  nk_coarse = product(wvmesh)
  ninterp = 4

  allocate(widc(nk_coarse, 6), idc(nk_coarse, 9), &
       ksint(nk_coarse, 3), f_coarse(nk_coarse), f_interp(ninterp), ik_interp(ninterp))
  do ik = 1, nk_coarse
     call demux_vector(ik, ksint(ik, :), wvmesh, 0_i64)
  end do
  call precompute_interpolation_corners_and_weights(wvmesh, &
       mesh_ref_array, ksint, idc, widc)

  call random_number(f_coarse)
  f_coarse = 2.0_r64*f_coarse - 1.0_r64

  ik_interp = [1, 2, 3, 4]

  do ik = 1, ninterp
     call interpolate_using_precomputed(idc(ik_interp(ik), :), widc(ik_interp(ik), :), &
          f_coarse, f_interp(ik))
  end do
  
  call test_array(itest)%assert(&
       f_interp, &
       [f_coarse(1), (2.0_r64*f_coarse(1) + f_coarse(2))/3.0_r64, (f_coarse(1) + 2.0_r64*f_coarse(2))/3.0_r64, f_coarse(2)], &
       tol = 1.0e-12_r64)

  !TODO Bilinear

  !TODO Trilinear
  
  !Pade_coeffs & Pade_continued
  itest = itest + 1
  test_array(itest) = testify("Pade approximant")
  allocate(im_axis(10), real_func(10))
  call linspace(im_axis, 0.0_r64, 1.0_r64, 10_i64)
  real_func = 1.0_r64/(im_axis - 1.0_r64)
  call test_array(itest)%assert(&
       Pade_continued(oneI*im_axis, real_func, [0.0_r64, 0.5_r64, 1.0_r64]), &
       [-1.0_r64 + 0.0_r64*oneI, -0.8_r64 + 0.4_r64*oneI, -0.5_r64 + 0.5_r64*oneI], &
       tol = 1.0e-10_r64)

  !.umklapp.
  q1(:, 1) = [0.5_r64, 0.0_r64, 0.0_r64]
  q1(:, 2) = [0.5_r64, 0.5_r64, 0.0_r64]
  q1(:, 3) = [0.5_r64, 0.5_r64, 0.5_r64]
  q1(:, 4) = [0.1_r64, 0.2_r64, 0.0_r64]
  q2(:, 1) = [0.5_r64, 0.0_r64, 0.0_r64]
  q2(:, 2) = [0.5_r64, 0.5_r64, 0.0_r64]
  q2(:, 3) = [0.5_r64, 0.5_r64, 0.5_r64]
  q2(:, 4) = -[0.5_r64, 0.5_r64, 0.1_r64]
  q3(:, 1) = [0.0_r64, 0.0_r64, 0.0_r64]
  q3(:, 2) = [0.0_r64, 0.0_r64, 0.0_r64]
  q3(:, 3) = [0.0_r64, 0.0_r64, 0.0_r64]
  q3(:, 4) = [0.6_r64, 0.7_r64, 0.9_r64]
  
  itest = itest + 1
  test_array(itest) = testify(".umklapp.")

  call test_array(itest)%assert( &
       [ &
       q1(:, 1) .umklapp. q2(:, 1), &
       q1(:, 2) .umklapp. q2(:, 2), &
       q1(:, 3) .umklapp. q2(:, 3), &
       q1(:, 4) .umklapp. q2(:, 4)], &
       [q3(:, 1), q3(:, 2), q3(:, 3), q3(:, 4)], tol = 1.0e-12_r64)

  !.umklapp. elemental
  itest = itest + 1
  test_array(itest) = testify("elemental .umklapp.")  
  call test_array(itest)%assert(pack(q1 .umklapp. q2, .true.), reshape(q3, [size(q3)]))

  !shrink int
  itest = itest + 1
  test_array(itest) = testify("shrink integer array")
  allocate(array_of_ints(5))
  array_of_ints = [1, 2, 3, 4, 5]*1_i64
  call shrink(array_of_ints, 2_i64)
  call test_array(itest)%assert(array_of_ints, [1, 2]*1_i64)

  !shrink real
  itest = itest + 1
  test_array(itest) = testify("shrink real array")
  allocate(array_of_reals(5))
  array_of_reals = [1, 2, 3, 4, 5]*1.0_r64
  call shrink(array_of_reals, 2_i64)
  call test_array(itest)%assert(array_of_reals, [1, 2]*1.0_r64)
  
  tests_all = testify(test_array)
  call tests_all%report
  
  if(tests_all%get_status() .eqv. .false.) error stop -1
end program test_misc

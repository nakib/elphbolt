program test_misc

  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use params, only: pi, kB, oneI, twopi
  use misc, only: int_div, expi, trace, kronecker, sort, cross_product, &
       twonorm, binsearch, mux_vector, demux_vector, interpolate, coarse_grained, &
       unique, linspace, compsimps, mux_state, demux_state, demux_mesh, expm1, &
       Fermi, Bose, Pade_continued, precompute_interpolation_corners_and_weights, &
       interpolate_using_precomputed, operator(.umklapp.), shrink, Hilbert_transform, &
       fft_next_pow2, interpolator_1d, permutations, lex_less_2d, lex_less_1d, &
       map_triplet_full_to_reduced, Jacobian

  implicit none

  !   external :: pgbeg, pgenv, pgline, pgend

  integer :: itest, dim, nb, nk
  integer, parameter :: num_tests = 49
  type(testify) :: test_array(num_tests), tests_all
  integer(i64) :: index, quotient, remainder, int_array(5), v1(3), v2(3), &
       v1_muxed, v2_muxed, ik, ik1, ik2, ik3, ib1, ib2, ib3, wvmesh(3), &
       mesh_ref_array(3), nk_coarse, ninterp, N, i, j, l, a(3,2), b(3,2), c(3), d(3), &
       nbands, mesh_size(3), kmesh(3)
  integer(i64), allocatable :: index_mesh_0(:, :), index_mesh_1(:, :), &
       ksint(:, :), idc(:, :), ik_interp(:), array_of_ints(:), perm(:, :), &
       lambda1_list(:), lambda2_list(:), M(:, :, :), ind_even(:), ind_odd(:), indexlist(:)
  real(r64) :: pauli1(2, 2), ipauli2(2, 2), pauli3(2, 2), &
       real_array(5), result, q1(3, 4), q2(3, 4), q3(3, 4), &
       lattvecs(3,3), b_matrix(3,3), k_c(3), k_f(3), kc1
  real(r64), allocatable :: integrand(:), domain(:), im_axis(:), real_func(:), &
       widc(:, :), f_coarse(:), f_interp(:), array_of_reals(:), x_axis(:)
  real(r64), allocatable :: hfx1_even(:), hfx1_odd(:), hfx2_even(:), hfx2_odd(:), &
       x_even(:), x_odd(:), xmin, xmax, & 
       f(:, :, :), gradf(:, :, :, :), correct_gradf(:, :, :, :)
  integer(i64) :: n_even, n_odd, n_hilb

  logical :: lex_order_2d, lex_order_1d, blocks

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

  !permutations(N)
  itest = itest + 1
  test_array(itest) = testify("permutations of 3 elements")
  N = 3
  perm = permutations(N)
  call test_array(itest)%assert(reshape(perm, [6*3]), [ &
       1, 2, 3,  1, 3, 2,  3, 1, 2,  3, 2, 1,  2, 3, 1,  2, 1, 3 &
       ]*1_i64)

  itest =  itest + 1
  test_array(itest) = testify("permutations of 4 elements")
  N = 4
  perm = permutations(N)
  call test_array(itest)%assert(reshape(perm, [24*4]), [ &
       1, 2, 3, 4,  1, 2, 4, 3,  1, 4, 2, 3,  4, 1, 2, 3,  4, 1, 3, 2,  1, 4, 3, 2, &
       1, 3, 4, 2,  1, 3, 2, 4,  3, 1, 2, 4,  3, 1, 4, 2,  3, 4, 1, 2,  4, 3, 1, 2, &
       4, 3, 2, 1,  3, 4, 2, 1,  3, 2, 4, 1,  3, 2, 1, 4,  2, 3, 1, 4,  2, 3, 4, 1, &
       2, 4, 3, 1,  4, 2, 3, 1,  4, 2, 1, 3,  2, 4, 1, 3,  2, 1, 4, 3,  2, 1, 3, 4 &
       ]*1_i64)

  !lex_less_2d(a, b)
  itest = itest + 1
  test_array(itest) = testify("lex_less_2d: a smaller than b")
  a = reshape([0, 1,  1, 0,  2, 0]*1_i64, [3, 2])
  b = reshape([0, 1,  1, 0,  2, 1]*1_i64, [3, 2])
  lex_order_2d = lex_less_2d(a, b)
  call test_array(itest)%assert(lex_order_2d, .true.)

  itest = itest + 1
  test_array(itest) = testify("lex_less_2d: a equal to b")
  a = reshape([0, 1,  1, 0,  2, 1]*1_i64, [3, 2])
  b = reshape([0, 1,  1, 0,  2, 1]*1_i64, [3, 2])
  lex_order_2d = lex_less_2d(a, b)
  call test_array(itest)%assert(lex_order_2d, .false.)

  itest = itest + 1
  test_array(itest) = testify("lex_less_2d: a greater than b")
  a = reshape([1, 1,  2, 2,  2, 1]*1_i64, [3, 2])
  b = reshape([1, 1,  2, 1,  2, 2]*1_i64, [3, 2])
  lex_order_2d = lex_less_2d(a, b)
  call test_array(itest)%assert(lex_order_2d, .false.)

  !lex_less_1d(c, d)
  itest = itest + 1
  test_array(itest) = testify("lex_less_1d: c smaller than d")
  c = reshape([0, 1, 2]*1_i64, [3])
  d = reshape([0, 1, 3]*1_i64, [3])
  lex_order_1d = lex_less_1d(c, d)
  call test_array(itest)%assert(lex_order_1d, .true.)

  ! map_triplet_full_to_reduced(ilambda2, ilambda1) = canonical representative of triplet (lambda1, lambda2, lambda3)
  itest = itest + 1
  test_array(itest) = testify("map_triplet_full_to_reduced: M(1, 1) should be [1, 1, 1]")
  nbands = 2
  mesh_size = [2, 2, 2]*1_i64
  allocate(lambda1_list(2), lambda2_list(16))
  lambda1_list = [1, 2]*1_i64
  lambda2_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]*1_i64
  call map_triplet_full_to_reduced(nbands, mesh_size, lambda1_list, lambda2_list, M)
  call test_array(itest)%assert(reshape(M(:, 1, 1), [3]), [1, 1, 1]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("map_triplet_full_to_reduced: M(9, 1) should be [1, 9, 9]")
  call test_array(itest)%assert(reshape(M(:, 9, 1), [3]), [1, 9, 9]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("map_triplet_full_to_reduced: M(16, 2) should be [2, 15, 16]")
  call test_array(itest)%assert(reshape(M(:, 16, 2), [3]), [2, 15, 16]*1_i64)

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

  ! Next power of 2
  itest = itest + 1
  test_array(itest) = testify("next power of 2")
  call test_array(itest)%assert(4096_i64, fft_next_pow2(3000_i64))

  ! Hilbert transform tests (H)
  ! fx1 -> function 1, fx2 -> function 2
  ! hfx1_even stores hilbert transform calculated for fx1, and for even number
  ! of points
  ! New unit test
  xmin = -100.0
  xmax = 100.0
  n_even = 4000
  allocate(x_even(n_even))
  call linspace(x_even, xmin, xmax, n_even)
  n_odd = n_even + 1
  allocate(x_odd(n_odd))
  call linspace(x_odd, xmin, xmax, n_odd)
  ! ind_even are indices to compare in case of even number of points
  ! ind_odd are indices to compare in case of odd number of points
  allocate(ind_even(6),ind_odd(5))
  ind_even = [1900, 1950, 2000, 2001, 2251, 2201]
  ind_odd = [1900, 1950, 2001, 2250, 2200] 

  itest = itest + 1
  test_array(itest) = testify("Hilbert transform: f(x) = 1/(1 + x^4), even points")
  allocate(hfx1_even(n_even))
  call Hilbert_transform(fx1(x_even), hfx1_even)
  call test_array(itest)%assert(hfx1_even(ind_even), hfx1(x_even(ind_even)), &
       tol = 7e-4_r64)

  itest = itest + 1
  test_array(itest) = testify("Hilbert transform: f(x) = sin(x)/(1 + x^2), even points")
  allocate(hfx2_even(n_even))
  call Hilbert_transform(fx2(x_even), hfx2_even)
  call test_array(itest)%assert(hfx2_even(ind_even), hfx2(x_even(ind_even)), &
       tol = 3e-5_r64)

  itest = itest + 1
  test_array(itest) = testify("Hilbert transform: f(x) = 1/(1 + x^4), odd points")
  allocate(hfx1_odd(n_odd))
  call Hilbert_transform(fx1(x_odd), hfx1_odd)
  call test_array(itest)%assert(hfx1_odd(ind_odd), hfx1(x_odd(ind_odd)), &
       tol = 7e-4_r64)

  itest = itest + 1
  test_array(itest) = testify("Hilbert transform: f(x) = sin(x)/(1 + x^2), odd points")
  allocate(hfx2_odd(n_odd))
  call Hilbert_transform(fx2(x_odd), hfx2_odd)
  call test_array(itest)%assert(hfx2_odd(ind_odd), hfx2(x_odd(ind_odd)), &
       tol = 3e-5_r64)

  !$!   ! Old unit test
  !$!   xmin = -30.0
  !$!   xmax = 30.0
  !$!   n_even = 4000
  !$!   n_odd = 4001
  !$!   ! ind_even are indices to compare in case of even number of points
  !$!   ! ind_odd are indices to compare in case of odd number of points
  !$!   allocate(ind_even(6),ind_odd(5))
  !$!   ind_even = [801, 1201, 1601, 2001, 2401, 2801]
  !$!   ind_odd = [889, 1333, 1777, 2221, 2665]
  !$! 
  !$!   itest = itest + 1
  !$!   test_array(itest) = testify("Hilbert transform: f(x) = 1/(1 + x^2), even points")
  !$!   allocate(x_even(n_even), hfx1_even(n_even))
  !$!   call linspace(x_even, xmin, xmax, n_even)
  !$!   call Hilbert_transform(fx1(x_even), hfx1_even)
  !$!   call test_array(itest)%assert(hfx1_even(ind_even), hfx1(x_even(ind_even)), &
  !$!        tol = 2e-4_r64)
  !$! 
  !$!   itest = itest + 1
  !$!   test_array(itest) = testify("Hilbert transform: f(x) = sin(x)/(1 + x^2), even points")
  !$!   allocate(hfx2_even(n_even))
  !$!   call Hilbert_transform(fx2(x_even), hfx2_even)
  !$!   call test_array(itest)%assert(hfx2_even(ind_even), hfx2(x_even(ind_even)), &
  !$!        tol = 1e-4_r64)
  !$! 
  !$!   itest = itest + 1
  !$!   test_array(itest) = testify("Hilbert transform: f(x) = 1/(1 + x^2), odd points")
  !$!   allocate(x_odd(n_odd), hfx1_odd(n_odd))
  !$!   call linspace(x_odd, xmin, xmax, n_odd)
  !$!   call Hilbert_transform(fx1(x_odd), hfx1_odd)
  !$!   call test_array(itest)%assert(hfx1_odd(ind_odd), hfx1(x_odd(ind_odd)), &
  !$!        tol = 4e-4_r64)
  !$! 
  !$!   itest = itest + 1
  !$!   test_array(itest) = testify("Hilbert transform: f(x) = sin(x)/(1 + x^2), odd points")
  !$!   allocate(hfx2_odd(n_odd))
  !$!   call Hilbert_transform(fx2(x_odd), hfx2_odd)
  !$!   call test_array(itest)%assert(hfx2_odd(ind_odd), hfx2(x_odd(ind_odd)), &
  !$!        tol = 1e-5_r64)


  ! 1D Interpolation
  itest = itest + 1
  test_array(itest) = testify("1D Interpolation")
  ! Testing for function f(x) = 5x - 1
  call test_array(itest)%assert([9, 19, 29]*1.0_r64, interpolator_1d([2, 4, 6]*1.0_r64, &
       [1, 3, 5, 7]*1.0_r64, [4, 14, 24, 34]*1.0_r64))


  ! Jacobian test 

  nb = 1 
  allocate(indexlist(1))
  indexlist = 1_i64 
  lattvecs = reshape( [ ( (merge(1.0, 0.0, i==j), j=1,3), i=1,3 ) ], [3, 3] )
  kmesh = [3, 1, 1]
  deallocate(array_of_reals)
  allocate(f(kmesh(1), nb, 3), gradf(kmesh(1), nb, 3, 3), correct_gradf(kmesh(1), nb, 3, 3), &
       array_of_reals(3))
  f = 0.0_r64
  do i = 1, kmesh(1)
     f(i, 1, 1) = (i**2)*1.0_r64*twopi*2
     f(i, 1, 2) = (i**3)*1.0_r64*twopi*2
     f(i, 1, 3) = 2.0_r64
     array_of_reals = [-15.0_r64, 24.0_r64, -9.0_r64]
     correct_gradf(i, 1, 1, 1) =  array_of_reals(i)
     array_of_reals = [-57_r64, 78_r64,-21_r64]
     correct_gradf(i, 1, 1, 2) =  array_of_reals(i)
  end do
  itest = itest + 1
  test_array(itest) = testify("Jacobian test 1d mesh 3x1x1")

  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .false.)
  call test_array(itest)%assert(reshape(correct_gradf, [size(correct_gradf)]), reshape(gradf, [size(gradf)]), tol = 3e-5_r64)

  kmesh = [3, 3, 1]*1_i64
  deallocate(f, gradf, array_of_reals, correct_gradf)
  allocate(f(product(kmesh), nb, 3), gradf(product(kmesh), nb, 3, 3), correct_gradf(product(kmesh), nb, 3, 3), &
       array_of_reals(3))
  correct_gradf = 0.0_r64

  do i = 1, kmesh(1)
     do j = 1, kmesh(2)
        nk = mux_vector([i*1_i64, j*1_i64, 1*1_i64], kmesh, 1_i64)          
        f(nk, 1, 1) = (i**2)*1.0_r64*twopi*2
        f(nk, 1, 2) = (i**3)*1.0_r64*twopi*2
        array_of_reals = [-15.0_r64, 24.0_r64, -9.0_r64]
        correct_gradf(nk, 1, 1, 1) =  array_of_reals(i)
        array_of_reals = [-57_r64, 78_r64,-21_r64]
        correct_gradf(nk, 1, 1, 2) =  array_of_reals(i)
     end do
  end do
  itest = itest + 1

  test_array(itest) = testify("Jacobian test 2d mesh 3x3x1")
  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .false.)
  call test_array(itest)%assert(reshape(correct_gradf, [size(correct_gradf)]), reshape(gradf, [size(gradf)]), tol = 3e-5_r64)

  do i = 1, kmesh(1)
     do j = 1, kmesh(2)
        nk = mux_vector([i*1_i64, j*1_i64, 1*1_i64], kmesh, 1_i64)          
        f(nk, 1, 3) = (sin(pi*(j-1)/3.0_r64)**2)*twopi*2
        array_of_reals = [0.0_r64, 9.0/4.0*1.0_r64, -9.0/4.0*1.0_r64]
        correct_gradf(nk, 1, 2, 3) =  array_of_reals(j)
     end do
  end do
  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .false.)

  itest = itest + 1
  test_array(itest) = testify("Jacobian test 2d mesh 3x3x1, with additional component")
  call test_array(itest)%assert(reshape(correct_gradf, [size(correct_gradf)]), reshape(gradf, [size(gradf)]), tol = 3e-5_r64)

  kmesh = [3, 3, 4]*1_i64
  deallocate(f, gradf, array_of_reals, correct_gradf)
  allocate(f(product(kmesh), nb, 3), gradf(product(kmesh), nb, 3, 3), correct_gradf(product(kmesh), nb, 3, 3), array_of_reals(3))
  correct_gradf = 0.0_r64 

  do i = 1, kmesh(1)
     do j = 1, kmesh(2)
        do l = 1, kmesh(3)
           nk = mux_vector([i*1_i64, j*1_i64, l*1_i64], kmesh, 1_i64) 
           f(nk, 1, 1) = (i**2)*1.0_r64*twopi*2 + (cos(pi*(l-1)/4.0_r64)**2)*twopi*2
           f(nk, 1, 2) = (i**3)*1.0_r64*twopi*2 + (cos(pi*(l-1)/4.0_r64)**2)*twopi*2        
           f(nk, 1, 3) = (sin(pi*(j-1)/3.0_r64)**2)*twopi*2
           array_of_reals = [0.0_r64, 9.0/4.0*1.0_r64, -9.0/4.0*1.0_r64]
           correct_gradf(nk, 1, 2, 3) =  array_of_reals(j)
           array_of_reals = [-15.0_r64, 24.0_r64, -9.0_r64]
           correct_gradf(nk, 1, 1, 1) =  array_of_reals(i)
           array_of_reals = [-57_r64, 78_r64,-21_r64]
           correct_gradf(nk, 1, 1, 2) =  array_of_reals(i)
        end do
     end do
  end do

  deallocate(array_of_reals)
  allocate(array_of_reals(4))
  do i = 1, kmesh(1)
     do j = 1, kmesh(2)
        do l = 1, kmesh(3)
           nk = mux_vector([i*1_i64, j*1_i64, l*1_i64], kmesh, 1_i64) 
           array_of_reals = [0.0_r64, -4.0_r64, 0.0_r64, 4.0_r64]
           correct_gradf(nk, 1, 3, 1) =  array_of_reals(l)
           correct_gradf(nk, 1, 3, 2) =  array_of_reals(l)
        end do
     end do
  end do
  itest = itest + 1

  test_array(itest) = testify("Jacobian test 3d mesh 3x3x4")
  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .false.)
  call test_array(itest)%assert(reshape(correct_gradf, [size(correct_gradf)]), reshape(gradf, [size(gradf)]), tol = 3e-5_r64)

  deallocate(indexlist)
  allocate(indexlist(600), x_axis(600))
  kmesh = [1000, 1, 1]*1_i64 
  deallocate(f, gradf, array_of_reals, correct_gradf)
  allocate(f(size(indexlist), nb, 3), gradf(size(indexlist), nb, 3, 3), array_of_reals(3), correct_gradf(size(indexlist), nb, 3, 3))
  array_of_reals = [0.5_r64, 0.5_r64, 1.0_r64]
  b_matrix = twopi*lattvecs
  do ik1 = 1, size(indexlist)
     if (ik1 <= size(indexlist)/2) then 
        kc1 = twopi* real(ik1-1)/kmesh(1)
        nk = mux_vector([ik1*1_i64, 1_i64, 1_i64], kmesh, 1_i64)
     else 
        kc1 = twopi* real(ik1+300-1)/kmesh(1)
        nk = mux_vector([(ik1+300)*1_i64, 1_i64, 1_i64], kmesh, 1_i64)
     end if
     indexlist(ik1) = nk
     f(ik1, 1, 1) = cos(kc1*array_of_reals(1))**2 
     correct_gradf(ik1, 1, 1, 1) = -sin(2*kc1*array_of_reals(1))*array_of_reals(1)
     x_axis(ik1) = kc1
  end do

  itest = itest + 1
  test_array(itest) = testify("Jacobian for blocks = .true. ")

  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .true.)
  gradf(300, 1, 1, 1) = correct_gradf(300, 1, 1, 1)
  gradf(301, 1, 1, 1) = correct_gradf(301, 1, 1, 1)
  gradf(600, 1, 1, 1) = correct_gradf(600, 1, 1, 1)

  call test_array(itest)%assert(gradf(:, 1, 1, 1), correct_gradf(:, 1, 1, 1), tol = 3e-3_r64)

  !   open(unit=10, file='data_graph_ind.dat', status='replace')
  !   do i = 1, size(indexlist)
  !     write(10, *) x_axis(i), gradf(i, 1, 1, 1)
  !   end do
  !   close(10)

  !   open(unit=10, file='correct_graph_ind.dat', status='replace')
  !   do i = 1, size(indexlist)
  !     write(10, *) x_axis(i), correct_gradf(i, 1, 1, 1)
  !   end do
  !   close(10)

  kmesh = [100, 100, 300]*1_i64
  deallocate(f, gradf, array_of_reals, correct_gradf, x_axis)
  allocate(f(product(kmesh), nb, 3), gradf(product(kmesh), nb, 3, 3), correct_gradf(product(kmesh), nb, 3, 3), array_of_reals(3), &
       x_axis(kmesh(1)))
  array_of_reals = [0.5_r64, 0.5_r64, 1.0_r64]
  lattvecs = reshape([0, 1, 1, 1, 0, 1, 1, 1, 0], [3, 3])/2.0_r64
  b_matrix = twopi*reshape([-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3])
  do ik1 = 1, kmesh(1)
     do ik2 = 1, kmesh(2)
        do ik3 = 1, kmesh(3)
           k_f = [real(ik1-1, kind=r64), real(ik2-1, kind=r64),real(ik3-1, kind=r64)]/kmesh
           nk = mux_vector([ik1*1_i64, ik2*1_i64, ik3*1_i64], kmesh, 1_i64) 
           k_c = matmul(b_matrix, k_f)
           do i = 1, 3
              f(nk, 1, i) = cos(i*DOT_PRODUCT(k_c, array_of_reals))**2 
              do j = 1, 3
                 correct_gradf(nk, 1, j, i) = -i*sin(2*i*DOT_PRODUCT(k_c, array_of_reals))*array_of_reals(j)
              end do
           end do
        end do
     end do
     x_axis(ik1) = k_c(1)
  end do

  itest = itest + 1
  test_array(itest) = testify("Jacobian for fcc crystal test 3d mesh 100x100x300")

  call Jacobian(f, gradf, lattvecs, kmesh, indexlist, blocks = .false.)
  call test_array(itest)%assert(reshape(correct_gradf, [size(correct_gradf)]), reshape(gradf, [size(gradf)]), tol = 8e-2_r64)


  !   do i = 1, 3
  !      do j = 1, 3
  !           print *, "Max value Jij, i = ",i, "j =",j , maxval(abs(gradf(:, 1,j,i) - correct_gradf(:, 1,j,i)))
  !      end do
  ! end do

  !   open(unit=10, file='data_graph.dat', status='replace')
  !   do i = 1, kmesh(1)
  !     write(10, *) x_axis(i), gradf(i, 1, 1, 3)
  !   end do
  !   close(10)

  !   open(unit=10, file='correct_graph.dat', status='replace')
  !   do i = 1, kmesh(1)
  !     write(10, *) x_axis(i), correct_gradf(i, 1, 1, 3)
  !   end do
  !   close(10)

  tests_all = testify(test_array)              
  call tests_all%report                                

  if(tests_all%get_status() .eqv. .false.) error stop -1
  itest = itest + 1
contains    

  ! Some reference functions and their Hilbert transforms:
  !$!   pure elemental real(r64) function fx1(x)  ! Old
  !$!     real(r64), intent(in) :: x
  !$! 
  !$!     fx1 = 1/(1.0_r64 + x**2)
  !$!   end function fx1
  !$! 
  !$!   pure elemental real(r64) function hfx1(x)
  !$!     real(r64), intent(in) :: x
  !$! 
  !$!     hfx1 = x/(1.0_r64 + x**2)
  !$!   end function hfx1

  pure elemental real(r64) function fx1(x)
    real(r64), intent(in) :: x

    fx1 = 1/(1.0_r64 + x**4)
  end function fx1

  pure elemental real(r64) function hfx1(x)
    real(r64), intent(in) :: x

    hfx1 = x*(1.0_r64 + x**2)/(1.0_r64 + x**4)/sqrt(2.0_r64)
  end function hfx1

  pure elemental real(r64) function fx2(x)
    real(r64), intent(in) :: x

    fx2 = sin(x)/(1.0_r64 + x**2)
  end function fx2

  pure elemental real(r64) function hfx2(x)
    real(r64), intent(in) :: x

    hfx2 = (exp(-1.0_r64) - cos(x))/(1.0_r64 + x**2)
  end function hfx2


end program test_misc

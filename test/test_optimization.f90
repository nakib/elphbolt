program optimization_test
  use optimization, only: find_correction, coeffs_for_optimization, global_coeffs, global_coeffs_opp_B, fun
  use precision, only: r128, r64, i64
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry, symmetrize_3x3_tensor
  use misc, only: eye, compute_eigenvalues, F_norm
  use bte_module, only: fill_coeff, transport_coeffs 
  use testify_m, only : testify

  implicit none

  type(crystal) :: crys
  type(symmetry) :: sym
  type(transport_coeffs) :: trans
  real(r64) :: T, corr_threshold, x0(9), x_solution(9), dev, corr, x_solution_B(27), x_solution_B_2(18)
  integer(i64) :: i, mesh(3) 
  integer :: exit_code, itest
  integer, parameter :: num_tests = 5
  type(testify) :: test_array(num_tests), tests_all
 
  print*, '<<module optimization unit tests>>'

  call system('cp ./test/input_SiC.nml ./input.nml')

  !Set up crystal, Zincblende SiC, we will need symmetry elements, it has T_d crystallographic point group
  call crys%initialize(print_flag = .false.)

  mesh = [4_i64, 4_i64, 4_i64] 

  ! Calculate crystal and BZ symmetries
  ! call sym%calculate_symmetries(crys, mesh, print_flag = .false.)

  ! Initialize transport coefficients, we will need them for optimization
  call trans%initialize_ph(1_i64)
  call trans%initialize_el(1_i64)

  ! Test 1: sigma = kappa_el = kappa_ph = alpha_el = I, alpha_ph = 0, sigmaS = -I, T = 0, corr_threshold = 20%
  ! Output: x = I, dev = 200%, corr = 0%, because correction in higher than corr_treshold, so it is not applied

  itest = 1

  T = 1.0_r64
  trans%ph_kappa(1, :, :) = 1.0_r64 * eye(3_i64)
  trans%el_sigma(1, :, :) = 1.0_r64 * eye(3_i64)
  trans%el_sigmaS(1, :, :) = -1.0_r64 * eye(3_i64)
  trans%el_kappa0(1, :, :) = 1.0_r64 * eye(3_i64)
  trans%el_alphabyT(1, :, :) = 1.0_r64 * eye(3_i64)/T
  trans%ph_alphabyT(1, :, :) = 0.0_r64 * eye(3_i64)/T
   
  call fill_coeff(trans, T, global_coeffs)
  ! allocate(global_coeffs%crotations(3, 3, sym%nsymm))
  ! global_coeffs%crotations = sym%crotations

  ! Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, mesh, Bfield = [2.0_r64, 0.0_r64, 1.0_r64], print_flag = .false.)

  ! deallocate(global_coeffs%crotations)
  allocate(global_coeffs%crotations(3, 3, sym%nsymm_Bfield))
  global_coeffs%crotations = sym%crotations_Bfield

  ! Set global_coeffs_opp_B (Global coef -B)
! alpha_el
global_coeffs_opp_B%el_alpha = transpose(reshape([-43056.616546941281,       -679.73047841474624,       -9.1190372299476971E-013, &
                                         679.73047841474647,       -43056.616546941281,       -1359.4609568294959,     &
                                         8.2071335069529272E-013,   1359.4609568294945,       -43056.616546941288]*1.0_r64, [3,3]))

! alpha_ph (all zeros)
global_coeffs_opp_B%ph_alpha = 0.0_r64

! sigma
global_coeffs_opp_B%el_sigma = transpose(reshape([930078.19538570102,        15365.869642558719,        2.9180919135832629E-012, &
                                        -15365.869642558684,       930078.19538570102,        30731.739285117419,     &
                                        -1.1672367654333052E-011,  -30731.739285117397,       930078.19538570137]*1.0_r64, [3,3]))

! sigmaS
global_coeffs_opp_B%el_sigmaS = 0.0_r64



! el_kappa0
global_coeffs_opp_B%el_kappa0 = transpose(reshape([8.7169932921567774,       0.17317697833760359,       3.6189438642889487E-016, &
                                         -0.17317697833760320,      8.7169932921567757,       0.34635395667520641,   &
                                         -1.6449744837677038E-016, -0.34635395667520652,       8.7169932921567757]*1.0_r64, [3,3]))

! ph_kappa (all zeros)
global_coeffs_opp_B%ph_kappa = 0.0_r64
global_coeffs%T = 300.0_r64

! ============================================================
! Set global_coeffs (Global coef B)
! ============================================================

! alpha_el
global_coeffs%el_alpha = transpose(reshape([-43056.616546941281,        679.73047841474965,      -5.4714223379686175E-013, &
                                  -679.73047841474511,       -43056.616546941281,        1359.4609568294948,      &
                                  -5.4714223379686175E-013,  -1359.4609568294952,       -43056.616546941288]*1.0_r64, [3,3]))

! alpha_ph (all zeros)
global_coeffs%ph_alpha = 0.0_r64

! sigma
global_coeffs%el_sigma = transpose(reshape([930078.19538570102,       -15365.869642558700,        1.1672367654333052E-011, &
                                  15365.869642558733,        930078.19538570091,       -30731.739285117408,      &
                                  7.2952297839581577E-012,   30731.739285117419,        930078.19538570137]*1.0_r64, [3,3]))

! sigmaS
global_coeffs%el_sigmaS = transpose(reshape([-147.57328106271726,        3.8920964548355510,       -4.2111346784453218E-015, &
                                   -3.8920964548355550,       -147.57328106271726,        7.7841929096711100,      &
                                    2.1055673392226609E-015,  -7.7841929096711127,       -147.57328106271723]*1.0_r64, [3,3]))

! el_kappa0
global_coeffs%el_kappa0 = transpose(reshape([8.7169932921567792,      -0.17317697833760365,       -3.2899489675354077E-017, &
                                   0.17317697833760365,      8.7169932921567774,        -0.34635395667520624,      &
                                  -3.6189438642889487E-016,  0.34635395667520696,        8.7169932921567757]*1.0_r64, [3,3]))

! ph_kappa (all zeros)
global_coeffs%ph_kappa = 0.0_r64
    
    
    
    ! corr_threshold = 20.0_r64
    call find_correction(x_solution_B_2, dev, corr)

    print *, "Solution M_I", x_solution_B_2(1:9)

    print *, "Solution M_J", x_solution_B_2(10:18)

    print*, "corr =", corr
    print*, "dev =", dev

!   test_array(itest) = testify("Test 1: sigma = kappa_el = kappa_ph = alpha_el = I, alpha_ph = 0, sigmaS = -I, T = 1, &
!                               corr_threshold = 200%")

!   call test_array(itest)%assert([reshape(1.0_r64* eye(3_i64), [9]), 200.0_r64, 0.0_r64], [x_solution, dev, corr], tol = 1e-3_r64)

!   ! Test 2: The same coefficients, corr_threshold = 200%
!   ! Output: x = -0.5*I, dev = 50%, corr = 150%  

!   itest = itest + 1
!   corr_threshold = 200.0_r64
!   call find_correction(x_solution, dev, corr, corr_threshold)

!   test_array(itest) = testify("Test 2: The same coefficients, corr_threshold = 200%, no correction applied")
!   call test_array(itest)%assert([reshape(-0.5_r64* eye(3_i64), [9]), 50.0_r64, 150.0_r64], [x_solution, dev, corr], tol = 1e-3_r64)
  
  
!   ! Test 3: sigma = -I, no solution for the system of equations
!   ! Output: x = I, dev = 200%, corr = 0%
!   itest = itest + 1
!   global_coeffs%el_sigma = -global_coeffs%el_sigma
!   call find_correction(x_solution, dev, corr, corr_threshold)
  
!   test_array(itest) = testify("Test 3: sigma = -I, no solution for the system of equations")
!   call test_array(itest)%assert([reshape(1.0_r64* eye(3_i64), [9]), 200.0_r64, 0.0_r64], [x_solution, dev, corr], tol = 1e-3_r64)

!   ! Test 4: sigma = 1100000, kappa_el = 5.3, sigmaS = -3.5, alpha_el = 10500, kappa_ph = alpha_ph = 0, T = 300, corr_threshold = 200%
!   ! Output: x = 0, dev = 100%, corr = 100% 
!   itest = itest + 1

!   global_coeffs%el_sigma = 1100000.0_r64 * eye(3_i64)
!   global_coeffs%el_sigmaS = -3.5_r64 * eye(3_i64)
!   global_coeffs%ph_kappa = 0.0_r64 * eye(3_i64)
!   global_coeffs%el_kappa0 = 5.3_r64 * eye(3_i64)
!   global_coeffs%ph_alpha = 0.0_r64 * eye(3_i64)
!   global_coeffs%el_alpha = 10500.0_r64 * eye(3_i64)
!   global_coeffs%T = 300.0_r64 
!   call find_correction(x_solution, dev, corr, corr_threshold)

!   test_array(itest) = testify("Test 4: sigma = 1100000, kappa_el = 5.3, sigmaS = -3.5, alpha_el = 10500, kappa_ph = alpha_ph = 0, &
!                                T = 300, corr_threshold = 200%")
!   call test_array(itest)%assert([reshape(0.0_r64 * eye(3_i64), [9]), 100.0_r64, 100.0_r64], [x_solution, dev, corr], tol = 1e-3_r64)


!   ! Test 5: change the symmetry of the structure, we will use the symmetry elements of SiC in the presence of B-field along z-axis, which has C_4v point group 
!   ! sigma = diag(10,10,20), kappa_el = diag(10,10,20), sigmaS = diag(1,1,2), alpha_el = diag(1,1,2.1), kappa_ph = diag(1,1,2), alpha_ph = 0, T = 1, corr_threshold = 200%
!   ! Output: x = diag(1,1,1.05), dev = 0%, corr = 5% 
!   itest = itest + 1
!   ! Calculate crystal and BZ symmetries
!   call sym%calculate_symmetries(crys, mesh, Bfield = [0.0_r64, 0.0_r64, 1.0_r64], print_flag = .false.)

!   global_coeffs%el_sigmaS = 1.0_r64*reshape([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0], [3, 3])
!   global_coeffs%el_sigma = 10.0_r64 * global_coeffs%el_sigmaS
!   global_coeffs%ph_kappa = global_coeffs%el_sigmaS
!   global_coeffs%el_kappa0 = 10.0_r64*global_coeffs%el_sigmaS
!   global_coeffs%ph_alpha = 0.0_r64 * eye(3_i64)
!   global_coeffs%el_alpha = reshape(1.0_r64*[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.1], [3, 3])
!   global_coeffs%T = 1.0_r64 

!   deallocate(global_coeffs%crotations)
!   allocate(global_coeffs%crotations(3, 3, sym%nsymm_Bfield))
!   global_coeffs%crotations = sym%crotations_Bfield
!   call find_correction(x_solution, dev, corr, corr_threshold)

!   test_array(itest) = testify("Test 5: C_4v symmetry, diag(10,10,20), kappa_el = diag(10,10,20), sigmaS = diag(1,1,2), &
!                               alpha_el = diag(1,1,2.1), kappa_ph = diag(1,1,2), alpha_ph = 0, T = 1, corr_threshold = 200%")
!   call test_array(itest)%assert([1.0_r64*[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.05], 0.0_r64, 5.0_r64], &
!                                   [x_solution, dev, corr], tol = 1e-3_r64)
  

! ! print*, "x_solution =", x_solution
! !   print*, "corr =", corr
! !   print*, "dev =", dev
  call system('rm  ./input.nml')

!   tests_all = testify(test_array)              
!   call tests_all%report                                

!   if(tests_all%get_status() .eqv. .false.) error stop -1

end program optimization_test
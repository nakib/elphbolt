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
  real(r64) :: T, corr_threshold, x0(9), x_solution(9), dev, corr, x_solution_B(27)
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
  global_coeffs_opp_B%el_alpha(1,:) = [-43049.656081174922_r64, -143.22163919973462_r64, 0.0_r64]
  global_coeffs_opp_B%el_alpha(2,:) = [143.22163919973462_r64, -43049.656081174922_r64, 0.0_r64]
  global_coeffs_opp_B%el_alpha(3,:) = [2.8421709430404007E-014_r64, 0.0_r64, -43049.656081174915_r64]

  ! alpha_ph (all zeros)
  global_coeffs_opp_B%ph_alpha = 0.0_r64

  ! sigma
  global_coeffs_opp_B%el_sigma(1,:) = [928848.77923708223_r64, 8434.5812517018712_r64, 0.0_r64]
  global_coeffs_opp_B%el_sigma(2,:) = [-8434.5812517018712_r64, 928848.77923708223_r64, 0.0_r64]
  global_coeffs_opp_B%el_sigma(3,:) = [0.0_r64, 0.0_r64, 928848.77923708246_r64]

  ! sigmaS
  global_coeffs_opp_B%el_sigmaS(1,:) = [-147.84190844556974_r64, 0.67199942806451363_r64, 2.2204460492503131E-016_r64]
  global_coeffs_opp_B%el_sigmaS(2,:) = [-0.67199942806451374_r64, -147.84190844556974_r64, 0.0_r64]
  global_coeffs_opp_B%el_sigmaS(3,:) = [0.0_r64, 0.0_r64, -147.84190844556974_r64]

  ! el_kappa0
  global_coeffs_opp_B%el_kappa0(1,:) = [8.7458132952234351_r64, -0.17289965476480373_r64, 0.0_r64]
  global_coeffs_opp_B%el_kappa0(2,:) = [0.17289965476480373_r64, 8.7458132952234351_r64, 0.0_r64]
  global_coeffs_opp_B%el_kappa0(3,:) = [0.0_r64, 0.0_r64, 8.7458132952234351_r64]

  ! ph_kappa (all zeros)
  global_coeffs_opp_B%ph_kappa = 0.0_r64

  ! ============================================================
  ! Set global_coeffs (Global coef B)
  ! ============================================================

  ! alpha_el
  global_coeffs%el_alpha(1,:) = [-43049.656081174922_r64, 1216.2393176297628_r64, -2.8421709430404007E-014_r64]
  global_coeffs%el_alpha(2,:) = [-1216.2393176297628_r64, -43049.656081174922_r64, 0.0_r64]
  global_coeffs%el_alpha(3,:) = [2.8421709430404007E-014_r64, 0.0_r64, -43049.656081174915_r64]

  ! alpha_ph (all zeros)
  global_coeffs%ph_alpha = 0.0_r64

  ! sigma
  global_coeffs%el_sigma(1,:) = [928848.77923708223_r64, -22297.158033415610_r64, 0.0_r64]
  global_coeffs%el_sigma(2,:) = [22297.158033415610_r64, 928848.77923708223_r64, 0.0_r64]
  global_coeffs%el_sigma(3,:) = [2.2737367544323206E-013_r64, 0.0_r64, 928848.77923708223_r64]

  ! sigmaS
  global_coeffs%el_sigmaS(1,:) = [-147.84190844556974_r64, 8.4561923377356241_r64, -2.2204460492503131E-016_r64]
  global_coeffs%el_sigmaS(2,:) = [-8.4561923377356241_r64, -147.84190844556974_r64, 0.0_r64]
  global_coeffs%el_sigmaS(3,:) = [-2.2204460492503131E-016_r64, 0.0_r64, -147.84190844556974_r64]

  ! el_kappa0
  global_coeffs%el_kappa0(1,:) = [8.7458132952234351_r64, -0.51925361144001059_r64, -1.3877787807814457E-017_r64]
  global_coeffs%el_kappa0(2,:) = [0.51925361144001059_r64, 8.7458132952234351_r64, 0.0_r64]
  global_coeffs%el_kappa0(3,:) = [0.0_r64, 0.0_r64, 8.7458132952234351_r64]

  ! ph_kappa (all zeros)
  global_coeffs%ph_kappa = 0.0_r64
  print* , "sigma(B) - sigma(-B)^T"
  print *, F_norm(global_coeffs%el_sigma - transpose(global_coeffs_opp_B%el_sigma))

  print* , "sigmaS(B) - alpha(-B)^T"
  print *, F_norm(global_coeffs%el_sigmaS*300.0 - transpose(global_coeffs_opp_B%el_alpha + global_coeffs_opp_B%ph_alpha))

  print* , "kappa(B) - kappa(-B)^T"
  print *, F_norm(global_coeffs%el_kappa0 +  global_coeffs%ph_kappa - transpose(global_coeffs_opp_B%el_kappa0 + &
       global_coeffs_opp_B%ph_kappa))

  ! corr_threshold = 20.0_r64
  call find_correction(x_solution_B, dev, corr)

  ! print *, "Solution", x_solution_B

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

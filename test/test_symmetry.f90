program test_symmetry

  use crystal_module, only: crystal
  use numerics_module, only: numerics
  use testify_m, only : testify
  use precision, only: r128, r64, i64
  use symmetry_module, only: symmetry, symmetrize_3x3_tensor
  use misc, only:   det_3x3

  implicit none

  type(crystal) :: crys
  type(symmetry) :: sym
  integer(i64) :: mesh(3), i, itest
  real(r64) :: Bfield(3), sigma(3, 3)
  integer, parameter :: num_tests = 8
  real, allocatable :: S3x3(:, :)
  real, parameter:: tensor(3,3) = reshape([ 11.0_r64, 12.0_r64, 13.0_r64, &
       21.0_r64, 22.0_r64, 23.0_r64, &
       31.0_r64, 32.0_r64, 33.0_r64], [3,3])
  type(testify) :: test_array(num_tests), tests_all

  print*, '<<module symmetry unit tests>>'

  ! Zincblende SiC in the absence of magnetic field

  itest = 1
  !Copy input file for SiC, 4x4x4 mesh
  call system('cp ./test/input_SiC.nml ./input.nml')

  !Set up crystal
  call crys%initialize

  mesh = [4_i64, 4_i64, 4_i64] 

  ! Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, mesh)

  ! 1st test: order of the crystallographic symmetry group
  test_array(itest) = testify("Zincblende SiC, order of the crystallographic symmetry group")
  call test_array(itest)%assert(sym%nsymm, 24_i64)

  ! 2nd test: order of the FBZ + time symmetry group
  itest = itest + 1
  test_array(itest) = testify("Zincblende SiC, order of the FBZ + time symmetry group")
  call test_array(itest)%assert(sym%nsymm_rot, 48_i64)

  ! 3rd test: symmetrization of a 3x3 tensor for T_d group
  sigma = tensor
  call symmetrize_3x3_tensor(sigma, sym%crotations)

  itest = itest + 1
  test_array(itest) = testify("Symmetrization of 3x3 tensor for T_d group")
  call test_array(itest)%assert(reshape(sigma, [9]), [sigma(1, 1), 0.0_r64, 0.0_r64, &
       0.0_r64, sigma(1, 1), 0.0_r64, &
       0.0_r64, 0.0_r64, sigma(1, 1)], tol = 1e-5_r64)

  ! Zincblende SiC in the presence of magnetic field along z-axis

  Bfield = [0.0_r64, 0.0_r64, 1.0_r64] 

  ! Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, mesh, Bfield)

  ! 4th test: symmetry group with B-field along z-axis
  allocate(S3x3(3, 3))
  S3x3 = 0.0_r64
  do i = 1, sym%nsymm_Bfield
     S3x3 = S3x3 + det_3x3(sym%crotations_Bfield(:, :, i))*sym%crotations_Bfield(:, :, i)
  end do
  S3x3 = S3x3 / real(sym%nsymm_Bfield, r64)

  itest = itest + 1
  test_array(itest) = testify("Zincblende SiC, total symmetry group with B-field along z-axis")
  call test_array(itest)%assert(matmul(S3x3, Bfield), Bfield, tol = 1e-5_r64)  

  ! 5th test: symmetrization of a 3x3 tensor for S_4 group
  sigma = tensor
  call symmetrize_3x3_tensor(sigma, sym%crotations_Bfield)

  itest = itest + 1
  test_array(itest) = testify("Symmetrization of 3x3 tensor for S_4 group")
  call test_array(itest)%assert(reshape(sigma, [9]), [16.5_r64, -4.5_r64, 0.0_r64, &
       4.5_r64, 16.5_r64, 0.0_r64, &
       0.0_r64, 0.0_r64, 33.0_r64], tol = 1e-5_r64)

  ! 6th test: symmetry group with B-field along diagonal (1,1,1)
  Bfield = [1.0_r64, 1.0_r64, 1.0_r64]
  call sym%calculate_symmetries(crys, mesh, Bfield)
  S3x3 = 0.0_r64
  do i = 1, sym%nsymm_Bfield
     S3x3 = S3x3 + det_3x3(sym%crotations_Bfield(:, :, i))*sym%crotations_Bfield(:, :, i)
  end do

  S3x3 = S3x3 / real(sym%nsymm_Bfield, r64)

  itest = itest + 1
  test_array(itest) = testify("Zincblende SiC, total symmetry group with B-field along diagonal")
  call test_array(itest)%assert(matmul(S3x3, Bfield), Bfield, tol = 1e-5_r64)

  ! 7th test: symmetrization of a 3x3 tensor for C_3 group with C_3 axis along (1,1,1)
  sigma = tensor
  call symmetrize_3x3_tensor(sigma, sym%crotations_Bfield) 
  itest = itest + 1
  test_array(itest) = testify("Symmetrization of 3x3 tensor for C_3 group")
  call test_array(itest)%assert(reshape(sigma, [9]), [1, 1, 1, & 
       1, 1, 1, &
       1, 1, 1]*22.0_r64, tol = 1e-5_r64)

  ! 8th test: symmetry group with B-field along (2,0,1)
  Bfield = [2.0_r64, 0.0_r64, 1.0_r64]
  call sym%calculate_symmetries(crys, mesh, Bfield)

  itest = itest + 1
  test_array(itest) = testify("Zincblende SiC, total symmetry group with B-field random")
  call test_array(itest)%assert(sym%nsymm_Bfield, 1_i64)

  ! Deallocate arrays
  deallocate(S3x3)

  !Remove input file 
  call system('rm  ./input.nml')

  tests_all = testify(test_array)              
  call tests_all%report                                

  if(tests_all%get_status() .eqv. .false.) error stop -1

end program test_symmetry

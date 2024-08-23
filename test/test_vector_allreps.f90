program test_vector_allreps

  use precision, only : r64, i64
  use testify_m, only : testify
  use vector_allreps_module, only : &
       vector_allreps, vector_allreps_add, vector_allreps_sub, &
       vector_allreps_print, vector_allreps_change_grid

  implicit none
  
  integer :: itest
  integer, parameter :: num_tests = 25 !17
  type(testify) :: test_array(num_tests), tests_all

  integer(i64) :: imuxed, grid(3), another_grid(3)
  real(r64) :: primitive_vecs(3, 3)
  !integer(i64), parameter :: N = 3
  type(vector_allreps) :: v0, v1, v2, v3, v0p, v1p, v2p

  print*, '<<module vector_allreps unit tests>>'

  primitive_vecs(:, 1) = [-0.5_r64, 0.0_r64, 0.5_r64]
  primitive_vecs(:, 2) = [ 0.0_r64, 0.5_r64, 0.5_r64]
  primitive_vecs(:, 3) = [-0.5_r64, 0.5_r64, 0.0_r64]
  
  grid = [4, 4, 4]

  !Null vector
  v0 = vector_allreps(1_i64, grid, primitive_vecs)
  
  itest = 1
  test_array(itest) = testify("null vector, integer rep")
  call test_array(itest)%assert(&
       v0%int, &
       [0, 0, 0]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("null vector, fractional rep")
  call test_array(itest)%assert(&
       v0%frac, &
       [0.0, 0.0, 0.0]*1.0_r64)

  itest = itest + 1
  test_array(itest) = testify("null vector, Cartesian rep")
  call test_array(itest)%assert(&
       v0%cart, &
       [0.0, 0.0, 0.0]*1.0_r64)

  !A non-null vector
  v1 = vector_allreps(2_i64, grid, primitive_vecs)

  itest = itest + 1
  test_array(itest) = testify("a non-null vector, integer rep")
  call test_array(itest)%assert(&
       v1%int, &
       [1, 0, 0]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("a non-null vector, fractional rep")
  call test_array(itest)%assert(&
       v1%frac, &
       [0.25, 0.0, 0.0]*1.0_r64)

  itest = itest + 1
  test_array(itest) = testify("a non-null vector, Cartesian rep")
  call test_array(itest)%assert(&
       v1%cart, &
       [-0.125, 0.0, 0.125]*1.0_r64)

  !The last vector
  v2 = vector_allreps(product(grid), grid, primitive_vecs)

  itest = itest + 1
  test_array(itest) = testify("another non-null vector, integer rep")
  call test_array(itest)%assert(&
       v2%int, &
       [3, 3, 3]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("another non-null vector, fractional rep")
  call test_array(itest)%assert(&
       v2%frac, &
       [0.75, 0.75, 0.75]*1.0_r64)

  itest = itest + 1
  test_array(itest) = testify("another non-null vector, Cartesian rep")
  call test_array(itest)%assert(&
       v2%cart, &
       [-0.75, 0.75, 0.75]*1.0_r64)

  !Addition with Umklapp
  v3 = vector_allreps_add(v1, v2, grid, primitive_vecs)

  itest = itest + 1
  test_array(itest) = testify("add vectors, muxed index")
  call test_array(itest)%assert(v3%muxed_index, 61_i64)
  
  itest = itest + 1
  test_array(itest) = testify("add vectors, integer rep")
  call test_array(itest)%assert(&
       v3%int, &
       [0, 3, 3]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("add vectors, fractional rep")
  call test_array(itest)%assert(&
       v3%frac, &
       [0.0, 3.0, 3.0]/4.0_r64)

  itest = itest + 1
  test_array(itest) = testify("add vectors, Cartesian rep")
  call test_array(itest)%assert(&
       v3%cart, &
       [-3.0, 6.0, 3.0]/8.0_r64)
  
  !Subtraction with Umklapp
  v3 = vector_allreps_sub(v1, v2, grid, primitive_vecs)

  itest = itest + 1
  test_array(itest) = testify("subtract vectors, muxed index")
  call test_array(itest)%assert(v3%muxed_index, 23_i64)
  
  itest = itest + 1
  test_array(itest) = testify("subtract vectors, integer rep")
  call test_array(itest)%assert(&
       v3%int, &
       [2, 1, 1]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("subtract vectors, fractional rep")
  call test_array(itest)%assert(&
       v3%frac, &
       [0.5, 0.25, 0.25]*1.0_r64)

  itest = itest + 1
  test_array(itest) = testify("subtract vectors, Cartesian rep")
  call test_array(itest)%assert(&
       v3%cart, &
       [-3.0, 2.0, 3.0]/8.0_r64)

  !Test grid change
  another_grid = [12, 12, 12] !3x first grid

  ! zero vector, everything should be unchanged
  v0p = vector_allreps_change_grid(v0, another_grid)

  itest = itest + 1
  test_array(itest) = testify("grid change null vector, muxed index")
  call test_array(itest)%assert(v0p%muxed_index, 1_i64)
  
  itest = itest + 1
  test_array(itest) = testify("grid change null vector, integer rep")
  call test_array(itest)%assert(&
       v0p%int, &
       [0, 0, 0]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("grid change null vector, fractional rep")
  call test_array(itest)%assert(&
       v0p%frac, &
       [0.0, 0.0, 0.0]*1.0_r64)

  itest = itest + 1
  test_array(itest) = testify("grid change null vector, Cartesian rep")
  call test_array(itest)%assert(&
       v0p%cart, &
       [0.0, 0.0, 0.0]*1.0_r64)
  
  ! test the last vector of original grid
  v1 = vector_allreps(product(grid), grid, primitive_vecs)
  v1p = vector_allreps_change_grid(v1, another_grid)

  itest = itest + 1
  test_array(itest) = testify("grid change last vector, muxed index")
  call test_array(itest)%assert(v1p%muxed_index, 1414_i64)
  
  itest = itest + 1
  test_array(itest) = testify("grid change last vector, integer rep")
  call test_array(itest)%assert(&
       v1p%int, &
       [9, 9, 9]*1_i64)

  itest = itest + 1
  test_array(itest) = testify("grid change last vector, fractional rep")
  call test_array(itest)%assert(&
       v1p%frac, &
       v1%frac)

  itest = itest + 1
  test_array(itest) = testify("grid change last vector, Cartesian rep")
  call test_array(itest)%assert(&
       v1p%cart, &
       v1%cart)
  
  tests_all = testify(test_array)
  call tests_all%report

  if(tests_all%get_status() .eqv. .false.) error stop -1
end program test_vector_allreps

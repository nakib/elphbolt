program test_task_manager

  use precision, only : r64, i64
  use testify_m, only : testify
  use task_manager_module, only: task_manager

  implicit none

  integer :: itest
  integer, parameter :: num_tests = 1
  type(testify) :: test_array(num_tests), tests_all

  type(task_manager) :: job1

  itest = 1
  test_array(itest) = testify("num batches and batches range")
  call job1%distribute_load(num_tasks = 100_i64, num_batches = 4_i64)
  call job1%print_report

  print*, job1%get_num_batches()
  print*, job1%get_batch_range(4_i64)

  call job1%distribute_load(num_tasks = 13333_i64, num_batches = 3_i64)
  call job1%print_report

  print*, job1%get_num_batches()
  print*, job1%get_batch_range(3_i64)
end program test_task_manager

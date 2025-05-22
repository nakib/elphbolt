program test_task_manager

  use precision, only : i64
  use testify_m, only : testify
  use task_manager_module, only : task_manager

  implicit none

  integer :: itest
  integer, parameter :: num_tests = 5
  type(testify) :: test_array(num_tests)

  type(task_manager) :: job1
  integer(i64) :: num_batches, start_idx, end_idx, batch_size
  integer(i64), dimension(3) :: batch_info

  print*, '<<module task manager unit tests>>'

  !Equal distribution: 100 tasks split into 4 batches
  itest = 1
  test_array(itest) = testify("equal distribution: 100 tasks, 4 batches")
  call job1%distribute_load(num_tasks = 100_i64, num_batches = 4_i64)
  batch_info = job1%get_batch_range(4_i64)
  num_batches = int(job1%get_num_batches(), i64)
  start_idx   = batch_info(1)
  end_idx     = batch_info(2)
  batch_size  = batch_info(3)
  call test_array(itest)%assert([num_batches, start_idx, end_idx, batch_size], &
       [4, 76, 100, 25]*1_i64)

  !odd distribution: 13333 tasks split into 3 batches
  itest = itest + 1
  test_array(itest) = testify("odd distribution: 13333 tasks, 3 batches")
  call job1%distribute_load(num_tasks = 13333_i64, num_batches = 3_i64)
  batch_info = job1%get_batch_range(3_i64)
  num_batches = int(job1%get_num_batches(), i64)
  start_idx   = batch_info(1)
  end_idx     = batch_info(2)
  batch_size  = batch_info(3)
  call test_array(itest)%assert([num_batches, start_idx, end_idx, batch_size], &
       [3, 8890, 13333, 4444]*1_i64)

  !Single task: 1 tasks split into 1 batches
  itest = itest + 1
  test_array(itest) = testify("single test: 1 tasks, 1 batches")
  call job1%distribute_load(num_tasks = 1_i64, num_batches = 1_i64)
  batch_info = job1%get_batch_range(1_i64)
  num_batches = int(job1%get_num_batches(), i64)
  start_idx   = batch_info(1)
  end_idx     = batch_info(2)
  batch_size  = batch_info(3)
  call test_array(itest)%assert([num_batches, start_idx, end_idx, batch_size], &
       [1, 1, 1, 1]*1_i64)

  !More batches than tasks: 3 tasks split into 5 batches
  itest = itest + 1
  test_array(itest) = testify("more batches than tasks: 3 tasks, 5 batches")
  call job1%distribute_load(num_tasks = 3_i64, num_batches = 5_i64)
  batch_info = job1%get_batch_range(3_i64)
  num_batches = int(job1%get_num_batches(), i64)
  start_idx   = batch_info(1)
  end_idx     = batch_info(2)
  batch_size  = batch_info(3)
  call test_array(itest)%assert([num_batches, start_idx, end_idx, batch_size], &
       [3, 3, 3, 1]*1_i64)
end program test_task_manager


program test_task_manager

  use precision, only : i64
  use testify_m, only : testify
  use task_manager_module, only : task_manager

  implicit none

  integer :: itest
  integer, parameter :: num_tests = 6
  type(testify) :: test_array(num_tests), tests_all

  type(task_manager) :: job1
  integer(i64) :: num_batches, batch, batch_info(3)

  print*, '<<task manager unit tests>>'

  !Equal distribution: 100 tasks split into 4 batches
  itest = 1
  test_array(itest) = testify("equal distribution: 100 tasks, 4 batches")

  call job1%distribute_load(num_tasks = 100_i64, num_batches = 4_i64, &
       restart_from_batch_record = .false., filename = "test_task_log1.txt")

  call test_array(itest)%assert([job1%get_num_batches(), job1%get_batch_range(4_i64)], &
       [4, 76, 100, 25]*1_i64)

  !Odd distribution: 13333 tasks split into 3 batches
  itest = itest + 1
  test_array(itest) = testify("odd distribution: 13333 tasks, 3 batches")

  call job1%distribute_load(num_tasks = 13333_i64, num_batches = 3_i64, &
       restart_from_batch_record = .false., filename = "test_task_log2.txt")

  call test_array(itest)%assert([job1%get_num_batches(), job1%get_batch_range(3_i64)], &
       [3, 8890, 13333, 4444]*1_i64)

  !Single task: 1 task split into 1 batch
  itest = itest + 1
  test_array(itest) = testify("single test: 1 task, 1 batch")

  call job1%distribute_load(num_tasks = 1_i64, num_batches = 1_i64, &
       restart_from_batch_record = .false., filename = "test_task_log3.txt")

  call test_array(itest)%assert([job1%get_num_batches(), job1%get_batch_range(1_i64)], &
       [1, 1, 1, 1]*1_i64)

  !More batches than tasks: 3 tasks split into 5 batches
  itest = itest + 1
  test_array(itest) = testify("more batches than tasks: 3 tasks, 5 batches")

  call job1%distribute_load(num_tasks = 3_i64, num_batches = 5_i64, &
       restart_from_batch_record = .false., filename = "test_task_log4.txt")

  call test_array(itest)%assert([job1%get_num_batches(), job1%get_batch_range(3_i64)], &
       [3, 3, 3, 1]*1_i64)

  !Write and read batch records: 100 tasks, 4 batches
  itest = itest + 1
  test_array(itest) = testify("write and read batch 2 record for 100 tasks and 4 batches")

  call job1%distribute_load(num_tasks = 100_i64, num_batches = 4_i64, &
       restart_from_batch_record = .false., filename = "test_task_log5.txt")

  call job1%write_record(1_i64)
  call job1%write_record(2_i64)
  call job1%read_record(batch, batch_info)
  call test_array(itest)%assert([batch, batch_info], [2, 26, 50, 25]*1_i64)

  !Write and read batch records: 5 tasks, 10 batches
  itest = itest + 1
  test_array(itest) = testify("write and read batch 3 record for 5 tasks and 10 batches")

  call job1%distribute_load(num_tasks = 5_i64, num_batches = 10_i64, &
       restart_from_batch_record = .false., filename = "test_task_log6.txt")

  call job1%write_record(1_i64)
  call job1%read_record(batch, batch_info)
  call test_array(itest)%assert([batch, batch_info], [1, 1, 1, 1]*1_i64)

  call job1%write_record(3_i64)
  call job1%read_record(batch, batch_info)
  call test_array(itest)%assert([batch, batch_info], [3, 3, 3, 1]*1_i64)

  tests_all = testify(test_array)
  call tests_all%report

  !Delete test_task_log*
  call system('rm test_task_log*')

  if(tests_all%get_status() .eqv. .false.) error stop -1
end program test_task_manager

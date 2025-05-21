module task_manager_module
  !! This module help to manage and distribute a set of tasks across a specified number of batches

  use precision, only: i64

  implicit none

  private
  public :: task_manager

  type :: task_manager
     private

     integer(i64) :: num_batches = 1
     integer(i64), allocatable :: batch_info(:, :)

   contains
     procedure, public :: distribute_load, print_report, get_num_batches, get_batch_range
  end type task_manager

contains

  subroutine distribute_load(self, num_tasks, num_batches)
    !!Partitions a total number of tasks into number of batches
    !!
    !!balancing the work as evenly as possible.
    !!the first few batches can receive one additional task if the 
    !!total number of task is not divisible evenly by the number of batches
    !!the number of batches used is limited to min(num_tasks, num_batches), mean no batch is left empty

    class(task_manager), intent(out) :: self
    integer(i64), intent(in) :: num_tasks, num_batches

    integer(i64) :: batch_size, residual_task, ibatch, start_idx, end_idx

    !this handles cases where batches are more than tasks
    self%num_batches = min(num_tasks, num_batches)
    allocate(self%batch_info(self%num_batches, 3))

    !divide the total number of tasks in num_batches (subtasks)
    batch_size = num_tasks/self%num_batches

    !how many tasks are left over
    residual_task = mod(num_tasks, self%num_batches)

    !the first batch index start is 1 always
    start_idx = 1
    !record the start and the end index and the size of each batch
    do ibatch = 1, self%num_batches
       !first check if this batch should receive one of the extra tasks
       if(ibatch <= residual_task) then
          end_idx = start_idx + batch_size
          !for the batches after the extra ones, assign exactly the batch_size tasks
       else
          end_idx = start_idx + batch_size - 1
       end if

       !save the start, end, and size of the batch
       self%batch_info(ibatch, 1) = start_idx
       self%batch_info(ibatch, 2) = end_idx
       self%batch_info(ibatch, 3) = end_idx - start_idx + 1

       !here the starting point for the next batch
       start_idx = end_idx + 1
    end do
  end subroutine distribute_load

  subroutine print_report(self)
    !!print the resume of the task distribution across different batches
    !!
    !!the total number of batches, the first and last indices for each batch, and number of task in each batch

    class(task_manager), intent(in) :: self

    if(this_image() == 1) then
       print*, "Number of batches:", self%num_batches
       print*, "First index of batch:", self%batch_info(:, 1)
       print*, "Last index of batch:", self%batch_info(:, 2)
       print*, "Size of batch:", self%batch_info(:, 3)
       print*, '========'
    end if
  end subroutine print_report

  pure integer function get_num_batches(self)
    !!returns the number of batches

    class(task_manager), intent(in) :: self

    get_num_batches = self%num_batches
  end function get_num_batches

  pure function get_batch_range(self, batch) result(batch_range)
    !!gives the range and size of a batch as a 3-element array: (start, end, size).

    class(task_manager), intent(in) :: self
    integer(i64), intent(in) :: batch
    integer(i64) :: batch_range(3)

    batch_range = self%batch_info(batch, :)
  end function get_batch_range

end module task_manager_module

module task_manager_module
        !! This module help to manage and distribute a set of tasks across a specified number of batches

  implicit none
  private
  public :: task_manager

  type :: task_manager
     private

     integer :: num_batches = 1
     integer, allocatable :: batch_info(:, :)
  
     contains
       procedure, public :: distribute_load
       procedure, public :: print_report
       procedure, public :: get_num_batches
       procedure, public :: get_batch_range
  end type task_manager

contains

  subroutine distribute_load(self, num_tasks, num_batches)
    class(task_manager), intent(out) :: self
    integer, intent(in) :: num_tasks, num_batches

    integer :: base_size, residual_task, ibatch, start_idx, end_idx

    self%num_batches = num_batches
    allocate(self%batch_info(num_batches, 3))

    !divide the total number of tasks in num_batches (subtasks)
    base_size = num_tasks / num_batches
    residual_task = mod(num_tasks, num_batches)
    start_idx = 1

    !record the start and the end index and the size of each batch
    do ibatch = 1, num_batches
       if (ibatch <= residual_task) then
          end_idx = start_idx + base_size
       else
          end_idx = start_idx + base_size - 1
       end if
       self%batch_info(ibatch, 1) = start_idx
       self%batch_info(ibatch, 2) = end_idx
       self%batch_info(ibatch, 3) = end_idx - start_idx + 1
       start_idx = end_idx + 1
    end do

  end subroutine distribute_load

  subroutine print_report(self)
    class(task_manager), intent(in) :: self

    if (this_image() == 1) then
            !if (.true.) then
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
    integer, intent(in) :: batch
    integer :: batch_range(3)

    batch_range = self%batch_info(batch, :)
  end function get_batch_range

end module task_manager_module

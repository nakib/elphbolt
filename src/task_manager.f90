module task_manager_module
  !! Module containing the data type related to task batching.

  use precision, only: i64

  implicit none

  private
  public :: task_manager

  type :: task_manager
     !! Container for task batching strategy that assign tasks to batches.

     private

     integer(i64) :: num_batches 
     !! Number of batches. 
     integer(i64), allocatable :: batch_info(:, :)
     !! Task batching information. 
     character(len = :), allocatable :: filename_record
     !! Saves batch records using write and read record

   contains

     procedure, public :: distribute_load, print_report, get_num_batches, get_batch_range, read_record, write_record
  end type task_manager

contains

  subroutine distribute_load(self, num_tasks, num_batches, filename)
    !! Partitions a total number of tasks into number of batches
    !! balancing the work as evenly as possible.
    !! We use a shuffling algorithm here.
    !! The number of batches used is limited to min(num_tasks, num_batches), which means no batch is left empty.

    class(task_manager), intent(out) :: self
    integer(i64), intent(in) :: num_tasks, num_batches
    character(len = *), intent(in) :: filename

    integer(i64) :: batch_size, residual_task, ibatch, start_idx, end_idx

    !This handles cases where batches are more than tasks.
    self%num_batches = min(num_tasks, num_batches)
    allocate(self%batch_info(self%num_batches, 3))

    !Clear old data
    if(allocated(self%filename_record)) deallocate(self%filename_record)
    allocate(character(len = len_trim(filename)) :: self%filename_record)
    self%filename_record = trim(filename)

    open(unit = 10, file = self%filename_record, status = 'replace', action = 'write')
    close(10)

    !Divide the total number of tasks in num_batches (subtasks).
    batch_size = num_tasks/self%num_batches

    !How many tasks are left over.
    residual_task = mod(num_tasks, self%num_batches)

    !The first batch index start is 1 always.
    start_idx = 1

    !Record the start and the end index and the size of each batch.
    do ibatch = 1, self%num_batches
       !First check if this batch should receive one of the extra tasks.
       if(ibatch <= residual_task) then
          end_idx = start_idx + batch_size
       else
          !For the batches after the extra ones, assign exactly the batch_size tasks.
          end_idx = start_idx + batch_size - 1
       end if

       !Save the start, end, and size of the batch.
       self%batch_info(ibatch, 1) = start_idx
       self%batch_info(ibatch, 2) = end_idx
       self%batch_info(ibatch, 3) = end_idx - start_idx + 1

       !Here the starting point for the next batch.
       start_idx = end_idx + 1
    end do
  end subroutine distribute_load

  subroutine print_report(self)
    !! Print the resume of the task distribution across different batches.
    !!
    !! Print the total number of batches, the first and last indices for each batch, and number of task in each batch.

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
    !! Returns the number of batches.

    class(task_manager), intent(in) :: self

    get_num_batches = self%num_batches
  end function get_num_batches

  pure function get_batch_range(self, batch) result(batch_range)
    !! Gives the range and size of a batch as a 3-element array: (start, end, size).

    class(task_manager), intent(in) :: self
    integer(i64), intent(in) :: batch
    integer(i64) :: batch_range(3)

    batch_range = self%batch_info(batch, :)
  end function get_batch_range

  subroutine write_record(self, batch_number)
    !! Append a line to the filename for the completed batch with timestamp, batch number, start, end  and size of batch.

    class(task_manager), intent(in) :: self
    integer(i64), intent(in) :: batch_number
    character(len = 30) :: timestamp
    ! Date_time in term of YYYY, MM, DD, HH, MM, SS.
    integer(i64) :: date_time(8)
    integer(i64) :: start_idx, end_idx, batch_size

    if(this_image() == 1) then
       !Formatting timestamp informations.
       call date_and_time(values = date_time)
       write(timestamp, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
            date_time(1), date_time(2), date_time(3), date_time(5), date_time(6), date_time(7)

       ! Append to the batch record file
       open(unit = 11, file = self%filename_record, status = 'old', position = 'append', action = 'write')
       write(11, '(A, 1X, I0, 1X, I0, 1X, I0, 1X, I0)') trim(timestamp), batch_number, &
            self%batch_info(batch_number, 1), self%batch_info(batch_number, 2), self%batch_info(batch_number, 3)
       close(11)

       call system("echo 'record saved' >> filename")
    end if
  end subroutine write_record

  subroutine read_record(self, batch_number, start_idx, end_idx, size_batch)
    !! Reads the last record from the filename record file.

    class(task_manager), intent(in) :: self
    integer(i64), intent(out) :: batch_number, start_idx, end_idx, size_batch

    !Holds the last line read from the record file
    character(len = 512) :: line
    character(len = 32) :: timestamp
    character(len = 20)  :: batch_str

    ! ios is an integer variable that stores the result for the read operation.
    integer :: ios
    integer(i64) :: last_batch_number, last_start, last_end, last_size_batch

    open(unit = 12, file = self%filename_record, status = 'old', action = 'read')
    do
       read(12, '(A)', iostat = ios) line
       ! Here exit when the end of file is reached
       if(ios /= 0) exit
    end do
    close(12)

    print *, "last line = '", trim(line), "'"
    read(line, *) timestamp, last_batch_number, last_start, last_end, last_size_batch

    batch_number = last_batch_number
    start_idx = last_start
    end_idx = last_end
    size_batch = last_size_batch

    print*, "Read from log    :"
    print*, "last batch number:", batch_number
    print*, "last sttart index:", start_idx
    print*, "last end index   :", end_idx
    print*, "last size batch  :", size_batch

    ! Convert batch number to string
    write(batch_str, '(I0)') batch_number

    if(this_image() == 1) then
       call system("echo 'Read last record: Batch " // trim(adjustl(batch_str)) // &
            "from" // trim(self%filename_record) // "' >> filename")
    end if
  end subroutine read_record
end module task_manager_module

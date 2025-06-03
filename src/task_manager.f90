module task_manager_module
  !! Module containing the data type related to task batching.

  use precision, only: i64

  implicit none

  private
  public :: task_manager

  type :: task_manager
     !! Container for task batching strategy that assign tasks to batches.

     private

     integer(i64) :: num_batches = 1 
     !! Number of batches.
     integer(i64), public :: num_finished_batches = 0
     !! Number of completed batches. 
     integer(i64), allocatable :: batch_info(:, :)
     !! Task batching information. 
     character(len = :), allocatable :: filename_record
     !! Saves batch records using write and read record

   contains

     procedure, public :: distribute_load, print_report, &
          get_num_batches, get_batch_range, read_record, write_record
  end type task_manager

contains

  subroutine distribute_load(self, num_tasks, num_batches, filename)
    !! Partitions a total number of tasks into number of batches balancing
    !! the work as evenly as possible. We use a shuffling algorithm here.
    !! The number of batches used is limited to min(num_tasks, num_batches),
    !! which means no batch is left empty.

    class(task_manager), intent(out) :: self
    integer(i64), intent(in) :: num_tasks, num_batches
    character(len = *), intent(in) :: filename

    integer(i64) :: batch_size, residual_task, ibatch, start_idx, end_idx, batch_number
    integer :: ios, unit
    logical :: file_exist
    character(len = 256) :: line, last_line
    character(len = 32) :: timestamp

    !This handles cases where batches are more than tasks.
    self%num_batches = min(num_tasks, num_batches)
    allocate(self%batch_info(self%num_batches, 3))

    !Default to zero for completed batches.
    self%num_finished_batches = 0

    !Check if file exists.
    inquire(file = filename, exist = file_exist)

    !This ensures future batch records can be appended without error.
    if(allocated(self%filename_record)) deallocate(self%filename_record)
    allocate(character(len = len_trim(filename)) :: self%filename_record)
    self%filename_record = trim(filename)

    !If the file exist, read the last line and update num_finished_tasks.
    if(file_exist) then
       open(newunit = unit, file = filename, status = "old", action = "read", position = "rewind", iostat = ios)
       if(ios == 0) then
          last_line = ''
          do

             read(unit, '(A)', iostat = ios) line
             !Exit on error or end of file.
             if(ios /= 0) exit
             !Keep updating with latest line.
             if (len_trim(line) > 0) last_line = line
          end do
          close(unit)

          !Verify that the batch record file contains at least one readable line.
          if (len_trim(last_line) > 0) then
             read(last_line, *, iostat = ios) timestamp, batch_number, start_idx, end_idx, batch_size
             if (ios == 0) then
                self%num_finished_batches = batch_number
             else
                self%num_finished_batches = 0
                if (this_image() == 1) print *, "No batches completed yet."
             end if
          else
             self%num_finished_batches = 0
             batch_number = 0
             if (this_image() == 1) print *, "Batch record file is empty. Starting from batch 0."
          end if


          !Read the last line for the completed batch number.
          read(last_line, *, iostat = ios) timestamp, batch_number, start_idx, end_idx, batch_size
          if(ios == 0) self%num_finished_batches = batch_number
       end if

       if(this_image() == 1) then
          print *, " Last record line : ", trim(last_line)
          print *, " Last batch number: ", batch_number
          print *, " Restart from last batch: " , self%num_finished_batches
       end if
    else

       open(newunit = unit, file = self%filename_record, status = 'replace', action = 'write')
       close(unit)
    end if

    !Divide the total number of tasks in num_batches (subtasks).
    batch_size = num_tasks/self%num_batches

    !How many tasks are left over?
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

       !Set the starting point for the next batch.
       start_idx = end_idx + 1
    end do
  end subroutine distribute_load

  subroutine print_report(self)
    !! Prints task distribution information.

    class(task_manager), intent(in) :: self

    if(this_image() == 1) then
       print*, "Number of batches:", self%num_batches
       print*, "First index of batch:", self%batch_info(:, 1)
       print*, "Last index of batch:", self%batch_info(:, 2)
       print*, "Size of batch:", self%batch_info(:, 3)
       print*, '========'
    end if
  end subroutine print_report

  pure integer(i64) function get_num_batches(self)
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
    !! Records completed batch information.

    class(task_manager), intent(in) :: self
    integer(i64), intent(in) :: batch_number

    character(len = 30) :: timestamp
    integer(i64) :: date_time(8) !date_time in YYYY, MM, DD, HH, MM, SS.
    integer :: unit
    integer(i64) :: start_idx, end_idx, batch_size

    if(this_image() == 1) then
       !Formatting timestamp informations.
       call date_and_time(values = date_time)
       write(timestamp, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
            date_time(1), date_time(2), date_time(3), date_time(5), date_time(6), date_time(7)

       !Append to the batch record file
       open(newunit = unit, file = self%filename_record, status = 'old', &
            position = 'append', action = 'write')
       write(unit, '(A, 1X, I0, 1X, I0, 1X, I0, 1X, I0)') trim(timestamp), &
            batch_number, self%batch_info(batch_number, :)
       close(unit)
    end if
  end subroutine write_record

  subroutine read_record(self, batch_number, batch_info)
    !! Reads the last completed batch information from record file.

    class(task_manager), intent(in) :: self
    integer(i64), intent(out) :: batch_number, batch_info(3)

    character(len = 512) :: line
    character(len = 32) :: timestamp
    integer :: ios, unit

    open(newunit = unit, file = self%filename_record, status = 'old', action = 'read')
    do
       read(unit, '(A)', iostat = ios) line
       ! Here exit when the end of file is reached
       if(ios /= 0) exit
    end do
    close(unit)

    read(line, *) timestamp, batch_number, batch_info(:)
  end subroutine read_record
end module task_manager_module

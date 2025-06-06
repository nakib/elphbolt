module task_manager_module
  !! Module containing the data type related to task batching.

  use precision, only: i64
  use misc, only: exit_with_message

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
    !! the work as evenly as possible.
    !
    ! We use a shuffling algorithm here.
    ! The number of batches used is limited to min(num_tasks, num_batches),
    ! which means no batch is left empty.

    class(task_manager), intent(out) :: self
    integer(i64), intent(in) :: num_tasks, num_batches
    character(len = *), intent(in) :: filename

    integer(i64) :: batch_size, residual_tasks, ibatch, start_idx, end_idx, batch_number
    integer :: ios, ios2, unit
    logical :: file_exists
    character(len = 256) :: line, last_line
    character(len = 32) :: timestamp

    !Set num_batches handling the case where there are more batches than tasks.
    self%num_batches = min(num_tasks, num_batches)

    !Allocate batch_info array
    allocate(self%batch_info(self%num_batches, 3))

    !Allocate the record filename
    allocate(character(len = len_trim(filename)) :: self%filename_record)

    !Set record filename
    self%filename_record = trim(filename)

    !Do all file creation business with just one image.
    if(this_image() == 1) then
       !Check if file exists.
       inquire(file = filename, exist = file_exists)

       !If the file exists, read the last line and update num_finished_tasks.
       if(file_exists) then
          open(newunit = unit, file = filename, status = "old", action = "read", &
               position = "rewind", iostat = ios)

          if(ios == 0) then !Was able to open the existing file
             !Defaults
             last_line = ''
             batch_number = 0

             do
                read(unit, '(A)', iostat = ios2) line

                !Exit on error or end of file.
                if(ios2 /= 0) exit

                !Keep updating with latest line.
                if(len_trim(line) > 0) last_line = line
             end do

             !Read last line
             if(len_trim(line) > 0) then !non-empty last line
                last_line = line

                read(last_line, *) &
                     timestamp, batch_number, start_idx, end_idx, batch_size

                !Here assert that the data in the file makes sense
                if(batch_number < 0 .or. start_idx < 0 &
                     .or. end_idx < 0 .or. batch_size < 0) then
                   close(unit)

                   call exit_with_message('Meaningless data in job record file. Exiting.')
                end if

                self%num_finished_batches = batch_number

                print *, " Last record line : ", trim(last_line)
                print *, " Last batch number: ", batch_number
                print *, " Restart from last batch: ", self%num_finished_batches
             else !empty last line
                print *, "Batch record file is empty."
             end if

             self%num_finished_batches = batch_number
          else !Error opening file
             call exit_with_message(&
                  'Could not open batch record file in distribute_load. Exiting.')
          end if
       else !file does not exist, so create it
          open(newunit = unit, file = self%filename_record, status = 'replace')
          write(unit, '(A, 1X, I0, 1X, I0, 1X, I0, 1X, I0)') "Start-marker", &
               0, 0, 0, 0
       end if

       close(unit)
    end if !only on image 1

    !Here broadcast self%num_finished_batches to all other images
    call co_broadcast(self%num_finished_batches, source_image = 1)

    !Divide the total number of tasks in num_batches (subtasks).
    batch_size = num_tasks/self%num_batches

    !How many tasks are left over?
    residual_tasks = mod(num_tasks, self%num_batches)

    !The first batch index start is 1 always.
    start_idx = 1

    !Record the start and the end index and the size of each batch.
    do ibatch = 1, self%num_batches
       !First check if this batch should receive one of the extra tasks.
       if(ibatch <= residual_tasks) then
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

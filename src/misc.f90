! Copyright (C) 2020- Nakib Haider Protik <nakib.haider.protik@gmail.com>
! This file is part of elphbolt <https://github.com/nakib/elphbolt>.
!
! elphbolt is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! elphbolt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with elphbolt. If not, see <http://www.gnu.org/licenses/>.

module misc
  !! Module containing miscellaneous math and numerics related functions and subroutines.

  use params, only: dp, k8, kB
  
  implicit none
  
  public
  private :: sort_int, sort_real

  type timer
     !! Container for timing related data and procedures.

     integer(k8) :: rate = 0
     integer(k8) :: start = -1
     integer(k8) :: end = -1
     character(len = 1024) :: event

   contains

     procedure :: start_timer, end_timer
  end type timer

  interface sort
     module procedure :: sort_int, sort_real
  end interface sort
  
contains

  subroutine start_timer(self, event)
    !! Start/Reset the timer. This is a blocking call.
    !! Only image 1 can modify timing information.

    class(timer), intent(out) :: self
    character(len = *), intent(in) :: event

    sync all
    if(this_image() == 1) then
       !Set the clock rate
       call system_clock(count_rate = self%rate)

       !Clock in
       call system_clock(count = self%start)

       !(Re)set the event name
       self%event = event
    end if
  end subroutine start_timer

  subroutine end_timer(self, event)
    !! End the timer and print the elapsed time. This is a blocking call.
    !! Only image 1 can modify timing information.

    class(timer), intent(inout) :: self
    character(len = *), intent(in) :: event

    !Local variable
    real(dp) :: time_elapsed

    sync all
    if(this_image() == 1) then
       !Clock in
       call system_clock(count = self%end)

       !Check the event name and if clock-in happened
       if((event /= self%event) .or. (self%start == -1_k8)) then
          call exit_with_message('Clock-in event does not match this clock-out event.')
       end if

       !Calculate and print time taken for this event
       time_elapsed = dble(self%end - self%start)/self%rate/3600.0_dp !hours
       !time_elapsed = dble(end - start)/rate/3600.0_dp*num_images() !cpu-hours
       write(*, "(A)") ".............."
       write(*, "(A, A, 1E16.8, A)") "| Timing info:", trim(event), time_elapsed, " hr"
       !write(*, "(A, A, 1E16.8, A)") "| Timing info:", trim(event), time_elapsed, " cpu-hrs"
       write(*, "(A)") ".............."
    end if
  end subroutine end_timer
  
  subroutine linspace(grid, min, max, num)
    !! Create equidistant grid.

    real(dp), allocatable, intent(out) :: grid(:)
    real(dp), intent(in) :: min, max
    integer(k8), intent(in) :: num

    !Local variables
    integer(k8) :: i
    real(dp) :: spacing

    !Allocate grid array
    allocate(grid(num))

    !Calculate grid spacing
    spacing = (max - min)/dble(num - 1)

    !Calculate grid
    do i = 1, num
       grid(i) = min + (i - 1)*spacing
    end do
  end subroutine linspace
  
  subroutine exit_with_message(message)
    !! Exit with error message.

    character(len = *), intent(in) :: message

    if(this_image() == 1) then
       write(*, "(A)") trim(message)
       stop
    end if
  end subroutine exit_with_message

  subroutine print_message(message)
    !! Print message.
    
    character(len = *), intent(in) :: message

    if(this_image() == 1) write(*, "(A)") trim(message)
  end subroutine print_message
  
  subroutine write2file_rank1_real(filename, data)
    !! Write rank-1 data to file.

    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: data(:)

    integer(k8) :: ik, nk

    nk = size(data(:))

    if(this_image() == 1) then
       open(1, file = trim(filename), status = "replace")
       do ik = 1, nk
          write(1, "(E20.10)") data(ik)
       end do
       close(1)
    end if
    sync all
  end subroutine write2file_rank1_real
  
  subroutine write2file_rank2_real(filename, data)
    !! Write rank-2 data to file.

    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: data(:,:)

    integer(k8) :: ik, nk
    character(len = 1024) :: numcols

    nk = size(data(:, 1))
    write(numcols, "(I0)") size(data(1, :))

    if(this_image() == 1) then
       open(1, file = trim(filename), status = "replace")
       do ik = 1, nk
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               data(ik, :)
       end do
       close(1)
    end if
    sync all
  end subroutine write2file_rank2_real

  subroutine write2file_rank3_real(filename, data)
    !! Write rank-3 data to file.

    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: data(:,:,:)

    integer(k8) :: ik, nk
    character(len = 1024) :: numcols

    nk = size(data(:, 1, 1))
    write(numcols, "(I0)") size(data(1, :, 1))*size(data(1, 1, :))

    if(this_image() == 1) then
       open(1, file = trim(filename), status = "replace")
       do ik = 1, nk
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               data(ik, :, :)
       end do
       close(1)
    end if
  end subroutine write2file_rank3_real

  subroutine write2file_response(filename, data, bandlist)
    !! Write list of vectors to band/branch resolved files.

    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: data(:,:,:)
    integer(k8), intent(in), optional :: bandlist(:)

    !Local variables
    integer(k8) :: ib, ibstart, ibend, nb, ik, nk, dim
    character(len = 1) :: numcols
    character(len = 1024) :: bandtag
    real(dp), allocatable :: aux(:,:)

    if(this_image() == 1) then
       nk = size(data(:, 1, 1))
       if(present(bandlist)) then
          nb = size(bandlist)
          ibstart = bandlist(1)
          ibend = bandlist(nb)
       else
          nb = size(data(1, :, 1))
          ibstart = 1
          ibend = nb
       end if
       dim = size(data(1, 1, :))
       write(numcols, "(I0)") dim

       !Band/branch summed
       open(1, file = trim(filename//"tot"), status = "replace")
       allocate(aux(nk, 3))
       aux = sum(data, dim = 2)
       do ik = 1, nk
          write(1, "(3(1E20.10),x)") aux(ik, :)
       end do
       close(1)

       !Band/branch resolved
       do ib = ibstart, ibend
          write(bandtag, "(I0)") ib
          open(2, file = trim(filename//bandtag), status = "replace")
          do ik = 1, nk
             write(2, "(3(1E20.10),x)") data(ik, ib, :)
          end do
          close(2)
       end do
    end if
    sync all
  end subroutine write2file_response

  subroutine readfile_response(filename, data, bandlist)
    !! Read list of vectors to band/branch resolved files.

    character(len = *), intent(in) :: filename
    real(dp), intent(out) :: data(:,:,:)
    integer(k8), intent(in), optional :: bandlist(:)

    !Local variables
    integer(k8) :: ib, ibstart, ibend, nb, ik, nk, dim
    character(len = 1) :: numcols
    character(len = 1024) :: bandtag
    
    nk = size(data(:, 1, 1))
    if(present(bandlist)) then
       nb = size(bandlist)
       ibstart = bandlist(1)
       ibend = bandlist(nb)
    else
       nb = size(data(1, :, 1))
       ibstart = 1
       ibend = nb
    end if
    dim = size(data(1, 1, :))
    write(numcols, "(I0)") dim

    !Band/branch resolved
    do ib = ibstart, ibend
       write(bandtag, "(I0)") ib
       open(1, file = trim(filename//bandtag), status = "old")
       do ik = 1, nk
          read(1, *) data(ik, ib, :)
       end do
       close(1)
    end do
    sync all
  end subroutine readfile_response
  
  subroutine append2file_transport_tensor(filename, it, data, bandlist)
    !! Append 3x3 tensor to band/branch resolved files.

    character(len = *), intent(in) :: filename
    integer(k8), intent(in) :: it
    real(dp), intent(in) :: data(:,:,:)
    integer(k8), intent(in), optional :: bandlist(:)

    !Local variables
    integer(k8) :: ib, nb, ibstart, ibend
    character(len = 1) :: numcols
    character(len = 1024) :: bandtag

    if(this_image() == 1) then
       if(present(bandlist)) then
          nb = size(bandlist)
          ibstart = bandlist(1)
          ibend = bandlist(nb)
       else
          nb = size(data(:, 1, 1))
          ibstart = 1
          ibend = nb
       end if

       write(numcols, "(I0)") 9

       !Band/branch summed
       if(it == 0) then
          open(1, file = trim(filename//"tot"), status = "replace")
       else
          open(1, file = trim(filename//"tot"), access = "append", status = "old")
       end if
       write(1, "(I3, " //trim(adjustl(numcols))//"E20.10)") &
            it, sum(data, dim = 1)
       close(1)

       !Band/branch resolved
       do ib = ibstart, ibend
          write(bandtag, "(I0)") ib
          if(it == 0) then
             open(2, file = trim(filename//bandtag), status = "replace")
          else
             open(2, file = trim(filename//bandtag), access = "append", status = "old")
          end if
          write(2, "(I3, " //trim(adjustl(numcols))//"E20.10)") &
               it, data(ib, :, :)
          close(2)
       end do
    end if
    sync all
  end subroutine append2file_transport_tensor

  subroutine write2file_spectral_tensor(filename, data, bandlist)
    !! Append 3x3 spectral transport tensor to band/branch resolved files.

    character(len = *), intent(in) :: filename
    real(dp), intent(in) :: data(:,:,:,:)
    integer(k8), intent(in), optional :: bandlist(:)

    !Local variables
    integer(k8) :: ie, ne, ib, nb, ibstart, ibend
    character(len = 1) :: numcols
    character(len = 1024) :: bandtag
    real(dp) :: aux(3,3)
    
    if(this_image() == 1) then
       !Number of energy points on grid
       ne = size(data(1, 1, 1, :))
       
       !Number of bands/branches and bounds
       if(present(bandlist)) then
          nb = size(bandlist)
          ibstart = bandlist(1)
          ibend = bandlist(nb)
       else
          nb = size(data(:, 1, 1, 1))
          ibstart = 1
          ibend = nb
       end if

       write(numcols, "(I0)") 9

       !Band/branch summed
       open(1, file = trim(filename//"tot"), status = "replace")
       do ie = 1, ne
          aux = 0.0_dp
          do ib = ibstart, ibend
             aux(:,:) = aux(:,:) + data(ib, :, :, ie)
          end do
          write(1, "("//trim(adjustl(numcols))//"E20.10)") aux
       end do
       close(1)

       !Band/branch resolved
       do ib = ibstart, ibend
          write(bandtag, "(I0)") ib
          open(2, file = trim(filename//bandtag), status = "replace")
          do ie = 1, ne
             write(2, "("//trim(adjustl(numcols))//"E20.10)") &
                  data(ib, :, :, ie)
          end do
          close(2)
       end do
    end if
    sync all
  end subroutine write2file_spectral_tensor

  subroutine int_div(num, denom, q, r)
    !! Quotient(q) and remainder(r) of the integer division num/denom.

    integer(k8), intent(in) :: num, denom
    integer(k8), intent(out) :: q, r

    q = num/denom
    r = mod(num, denom)
  end subroutine int_div

  subroutine distribute_points(npts, chunk, istart, iend, num_active_images)
    !! Distribute points among processes

    integer(k8), intent(in) :: npts
    integer(k8), intent(out) :: chunk, istart, iend, num_active_images

    !Number of active images
    num_active_images = min(npts, num_images())
    !Number of points per process
    chunk = ceiling(dble(npts)/num_images())
    !Start index
    istart = (this_image()-1)*chunk + 1
    !End index
    iend = min(chunk*this_image(), npts)
    !Update chunk
    if(istart < iend) then
       chunk = iend - istart + 1
    else
       chunk = 1
    end if
  end subroutine distribute_points
  
  function cross_product(A, B)
    !! Cross product of A and B.

    real(dp), intent(in) :: A(3), B(3)
    real(dp) :: cross_product(3)

    cross_product(1) = A(2)*B(3) - A(3)*B(2)
    cross_product(2) = A(3)*B(1) - A(1)*B(3)
    cross_product(3) = A(1)*B(2) - A(2)*B(1)
  end function cross_product

  pure complex(dp) function expi(x)
    !! Calculate exp(i*x) = cos(x) + isin(x)

    real(dp), intent(in) :: x

    expi = cmplx(cos(x), sin(x), dp)
  end function expi

  pure real(dp) function twonorm(v)
    !! 2-norm of a vector

    real(dp), intent(in) :: v(:)
    integer(k8) :: i, s

    s = size(v)
    twonorm = 0.0_dp
    do i = 1, s
       twonorm = v(i)**2 + twonorm
    end do
    twonorm = sqrt(twonorm)
  end function twonorm

  pure real(dp) function trace(mat)
    !! Trace of square matrix

    real(dp), intent(in) :: mat(:,:)
    integer(k8) :: i, dim

    dim = size(mat(:, 1))
    trace = 0.0_dp
    do i = 1, dim
       trace = trace + mat(i, i)
    end do
  end function trace

  subroutine sort_int(list)
    !! Swap sort list of integers

    integer(k8), intent(inout) :: list(:)
    integer(k8) :: i, j, n
    integer(k8) :: aux, tmp
    
    n = size(list)

    do i = 1, n
       aux = list(i)
       do j = i + 1, n
          if (aux > list(j)) then
             tmp = list(j)
             list(j) = aux
             list(i) = tmp
             aux = tmp
          end if
       end do
    end do
  end subroutine sort_int

  subroutine sort_real(list)
    !! Swap sort list of reals
    
    real(dp), intent(inout) :: list(:)
    real(dp) :: aux, tmp
    integer(k8) :: i, j, n

    n = size(list)

    do i = 1, n
       aux = list(i)
       do j = i + 1, n
          if (aux > list(j)) then
             tmp = list(j)
             list(j) = aux
             list(i) = tmp
             aux = tmp
          end if
       end do
    end do
  end subroutine sort_real

  subroutine binsearch(array, e, m)
    !! Binary search in a list of integers and return index.
    
    integer(k8), intent(in) :: array(:), e
    integer(k8), intent(out) :: m
    integer(k8) :: a, b, mid

    a = 1
    b = size(array)
    m = (b + a)/2
    mid = array(m)

    do while(mid /= e)
       if(e > mid) then
          a = m + 1
       else if(e < mid) then
          b = m - 1
       end if
       if(a > b) then
          m = -1
          exit
       end if
       m = (b + a)/2
       mid = array(m)
    end do
  end subroutine binsearch

  subroutine compsimps(f, h, s)
    !! Composite Simpson's rule for real function
    !! f integrand
    !! h integration variable spacing
    !! s result

    real(dp), intent(in) :: h, f(:)
    real(dp), intent(out) :: s

    !Local variables
    integer(kind=4) :: i, numint, n
    real(dp) :: a, b

    n = size(f)
    
    s = 0.0_dp

    a = f(1)

    !If n is odd then number of intervals is even, carry on.
    !Otherwise, do trapezoidal rule in the last interval.
    if(mod(n, 2) /= 0) then
       numint = n - 1
       b = f(n)
    else
       !Note: Number of sample points is even, so
       !I will do trapezoidal rule in the last interval.
       numint = n - 2
       b = f(n - 1)
    end if

    s = s + a + b

    !even sites
    do i = 2, numint, 2
       s = s + 4.0_dp*f(i)
    end do

    !odd sites
    do i = 3, numint, 2
       s = s + 2.0_dp*f(i)
    end do
    s = s*h/3.0_dp

    if(mod(n, 2) == 0) then
       !trapezoidal rule
       s = s + 0.5_dp*(f(n) + f(n - 1))*h
    end if
  end subroutine compsimps
  
  function mux_vector(v, mesh, base)
    !! Multiplex index of a single wave vector.
    !! v is the demultiplexed triplet of a wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    integer(k8), intent(in) :: v(3), mesh(3), base
    integer(k8) :: mux_vector

    if(base < 0 .or. base > 1) then
       call exit_with_message("Base has to be either 0 or 1 in misc.f90:mux_vector")
    end if

    if(base == 0) then
       mux_vector = (v(3)*mesh(2) + v(2))*mesh(1) + v(1) + 1
    else
       mux_vector = ((v(3) - 1)*mesh(2) + (v(2) - 1))*mesh(1) + v(1)
    end if
  end function mux_vector

  subroutine demux_vector(i, v, mesh, base)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.
    
    integer(k8), intent(in) :: i, mesh(3), base
    integer(k8), intent(out) :: v(3)
    integer(k8) :: aux

    if(base < 0 .or. base > 1) then
       call exit_with_message("Base has to be either 0 or 1 in misc.f90:demux_vector")
    end if

    call int_div(i - 1, mesh(1), aux, v(1))
    call int_div(aux, mesh(2), v(3), v(2))
    if(base == 1) then
       v(1) = v(1) + 1
       v(2) = v(2) + 1
       v(3) = v(3) + 1 
    end if
  end subroutine demux_vector
  
  subroutine demux_mesh(index_mesh, nmesh, mesh, base, indexlist)
    !! Demultiplex all wave vector indices 
    !! (optionally, from a list of indices).
    !! Internally uses demux_vector.

    integer(k8), intent(in) :: nmesh, mesh(3), base
    integer(k8), optional, intent(in) :: indexlist(nmesh)
    integer(k8), intent(out) :: index_mesh(3, nmesh)
    integer(k8) :: i

    do i = 1, nmesh !over total number of wave vectors
       if(present(indexlist)) then
          call demux_vector(indexlist(i), index_mesh(:, i), mesh, base)
       else
          call demux_vector(i, index_mesh(:, i), mesh, base)
       end if
    end do
  end subroutine demux_mesh

  pure integer(k8) function mux_state(nbands, iband, ik)
    !! Multiplex a (band index, wave vector index) pair into a state index 
    !!
    !! nbands is the number of bands
    !! iband is the band index
    !! ik is the wave vector index
    
    integer(k8), intent(in) :: nbands, ik, iband 

    mux_state = (ik - 1)*nbands + iband
  end function mux_state

  subroutine demux_state(m, nbands, iband, ik)
    !! Demultiplex a state index into (band index, wave vector index) pair
    !!
    !! m is the multiplexed state index
    !! nbands is the number of bands
    !! iband is the band index
    !! ik is the wave vector index
    
    integer(k8), intent(in) :: m, nbands
    integer(k8), intent(out) :: ik, iband 

    iband = modulo(m - 1, nbands) + 1
    ik = int((m - 1)/nbands) + 1
  end subroutine demux_state

  pure real(dp) function Bose(e, T)
    !! e Energy in eV
    !! T temperature in K

    real(dp), intent(in) :: e, T
    
    Bose = 1.0_dp/(exp(e/kB/T) - 1.0_dp)
  end function Bose

  pure real(dp) function Fermi(e, chempot, T)
    !! e Energy in eV
    !! chempot Chemical potential in eV
    !! T temperature in K

    real(dp), intent(in) :: e, chempot, T

    Fermi = 1.0_dp/(exp((e - chempot)/kB/T) + 1.0_dp)
  end function Fermi

  subroutine interpolate(coarsemesh, refinement, f, q, interpolation)
    !! Subroutine to perform BZ interpolation.
    !!
    !! coarsemesh The coarse mesh.
    !! refinement The mesh refinement factor.
    !! f The coarse mesh function to be interpolated.
    !! q The 0-based index vector where to evaluate f.
    !! interpolation The result
    
    integer(k8), intent(in) :: coarsemesh(3), q(3), refinement(3)
    real(dp), intent(in) :: f(:)
    real(dp), intent(out) :: interpolation
    
    integer(k8) :: info, r0(3), r1(3), ipol, mode, count
    integer(k8), allocatable :: pivot(:)
    integer(k8) :: i000, i100, i010, i110, i001, i101, i011, i111, equalpol
    real(dp) :: x0, x1, y0, y1, z0, z1, x, y, z, v(2), v0(2), v1(2)
    real(dp), allocatable :: T(:, :), c(:)
    real(dp) :: aux

    !External procedures
    external :: dgesv
    
    aux = 0.0_dp
    equalpol = 0_k8

    !Find on the coarse mesh the two diagonals.
    r0 = modulo(floor(q/dble(refinement)), coarsemesh)
    r1 = modulo(ceiling(q/dble(refinement)), coarsemesh)

    mode = 0
    do ipol = 1, 3
       if(r1(ipol) == r0(ipol)) then
          mode = mode + 1
       end if
    end do
    
    !mode = 0: 3d interpolation
    !mode = 1: 2d interpolation
    !mode = 2: 1d interpolation
    if(mode == 0) then !3d
       allocate(pivot(8), T(8, 8), c(8))

       !Fine mesh point
       x =  q(1)/dble(refinement(1)*coarsemesh(1))
       y =  q(2)/dble(refinement(2)*coarsemesh(2))
       z =  q(3)/dble(refinement(3)*coarsemesh(3))

       !Coarse mesh walls
       x0 = floor(q(1)/dble(refinement(1)))/dble(coarsemesh(1))
       y0 = floor(q(2)/dble(refinement(2)))/dble(coarsemesh(2))
       z0 = floor(q(3)/dble(refinement(3)))/dble(coarsemesh(3))
       x1 = ceiling(q(1)/dble(refinement(1)))/dble(coarsemesh(1))
       y1 = ceiling(q(2)/dble(refinement(2)))/dble(coarsemesh(2))
       z1 = ceiling(q(3)/dble(refinement(3)))/dble(coarsemesh(3))

       !Coarse mesh corners
       i000 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
       i100 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
       i010 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
       i110 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r1(1)+1
       i001 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
       i101 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
       i011 = (r1(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
       i111 = (r1(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r1(1)+1

       !Evaluate functions at the corners and form rhs    
       c = (/ f(i000), f(i100), f(i010), f(i110),&
            f(i001), f(i101), f(i011), f(i111) /)

       !Form the transformation matrix T
       T(1,:) = (/1.0_dp, x0, y0, z0, x0*y0, x0*z0, y0*z0, x0*y0*z0/)
       T(2,:) = (/1.0_dp, x1, y0, z0, x1*y0, x1*z0, y0*z0, x1*y0*z0/)
       T(3,:) = (/1.0_dp, x0, y1, z0, x0*y1, x0*z0, y1*z0, x0*y1*z0/)
       T(4,:) = (/1.0_dp, x1, y1, z0, x1*y1, x1*z0, y1*z0, x1*y1*z0/)
       T(5,:) = (/1.0_dp, x0, y0, z1, x0*y0, x0*z1, y0*z1, x0*y0*z1/)
       T(6,:) = (/1.0_dp, x1, y0, z1, x1*y0, x1*z1, y0*z1, x1*y0*z1/)
       T(7,:) = (/1.0_dp, x0, y1, z1, x0*y1, x0*z1, y1*z1, x0*y1*z1/)
       T(8,:) = (/1.0_dp, x1, y1, z1, x1*y1, x1*z1, y1*z1, x1*y1*z1/)

       !Solve Ta = c for a,
       !where c is an array containing the function values at the 8 corners.
       call dgesv(8,1,T,8,pivot,c,8,info)

       !Approximate f(x,y,z) in terms of a.
       aux = c(1) + c(2)*x + c(3)*y + c(4)*z +&
            c(5)*x*y + c(6)*x*z + c(7)*y*z + c(8)*x*y*z
    else if(mode == 1) then !2d
       allocate(pivot(4), T(4, 4), c(4))
       
       count = 1
       do ipol = 1, 3
          if(r1(ipol) .eq. r0(ipol)) then
             equalpol = ipol
          else
             v(count) = q(ipol)/dble(refinement(ipol)*coarsemesh(ipol))
             v0(count) = floor(q(ipol)/dble(refinement(ipol)))/dble(coarsemesh(ipol))
             v1(count) = ceiling(q(ipol)/dble(refinement(ipol)))/dble(coarsemesh(ipol))
             count = count+1 
          end if
       end do

       i000 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
       if(equalpol .eq. 1) then !1st 2 subindices of i are y,z
          i010 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
          i100 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
          i110 = (r1(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
       else if(equalpol .eq. 2) then !x,z
          i010 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
          i100 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
          i110 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
       else !x,y
          i010 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
          i100 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
          i110 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r1(1)+1
       end if

       c = (/ f(i000), f(i010), f(i100), f(i110)/)

       T(1,:) = (/1.0_dp, v0(1), v0(2), v0(1)*v0(2)/)
       T(2,:) = (/1.0_dp, v0(1), v1(2), v0(1)*v1(2)/)
       T(3,:) = (/1.0_dp, v1(1), v0(2), v1(1)*v0(2)/)
       T(4,:) = (/1.0_dp, v1(1), v1(2), v1(1)*v1(2)/)

       call dgesv(4,1,T,4,pivot,c,4,info)

       aux = c(1) + c(2)*v(1) + c(3)*v(2) + c(4)*v(1)*v(2)
    else !1d
       do ipol = 1, 3
          if(r1(ipol) /= r0(ipol)) then
             x =  q(ipol)/dble(refinement(ipol)*coarsemesh(ipol))
             x0 = floor(q(ipol)/dble(refinement(ipol)))/dble(coarsemesh(ipol))
             x1 = ceiling(q(ipol)/dble(refinement(ipol)))/dble(coarsemesh(ipol))

             i000 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
             if(ipol .eq. 1) then
                i100 = (r0(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r1(1)+1
             else if(ipol .eq. 2) then
                i100 = (r0(3)*coarsemesh(2)+r1(2))*coarsemesh(1)+r0(1)+1
             else
                i100 = (r1(3)*coarsemesh(2)+r0(2))*coarsemesh(1)+r0(1)+1
             end if
             aux = f(i000) + (x - x0)*(f(i100) - f(i000))/(x1 - x0)
          end if
       end do
    end if
    
    interpolation = aux
  end subroutine interpolate

  subroutine welcome
    !! Subroutine to print a pretty banner.
    
    if(this_image() == 1) then
       write(*,'(A75)') "+-------------------------------------------------------------------------+"
       write(*,'(A75)') "| \                                                                       |"
       write(*,'(A75)') "|  \                                                                      |"
       write(*,'(A75)') "|   \   \                                                                 |"
       write(*,'(A75)') "|    \   \                                                                |"
       write(*,'(A75)') "|   __\   \              _        _    _           _    _                 |"
       write(*,'(A75)') "|   \      \         ___|.|      |.|  | |__   ___ |.|_ / /__              |"
       write(*,'(A75)') "|    \    __\       / _ \.|   _  |.|_ | '_ \ / _ \|.|_  ___/              |"
       write(*,'(A75)') "|     \  \         |  __/.| |/ \_|/  \| |_) : (_) |.|/ /__                |"
       write(*,'(A75)') "|      \ \          \___|_|/|__/ |   /| ___/ \___/|_|\___/                |"
       write(*,'(A75)') "|       \ \                /|                                             |"
       write(*,'(A75)') "|        \\                \|                                             |"
       write(*,'(A75)') "|         \\                '                                             |"
       write(*,'(A75)') "|          \                                                              |"
       write(*,'(A75)') "|           \                                                             |"
       write(*,'(A75)') "| A solver for the coupled electron-phonon Boltzmann transport equations. |"
       write(*,'(A75)') "| Copyright (C) 2020- Nakib Haider Protik.                                |"
       write(*,'(A75)') "|                                                                         |"
       write(*,'(A75)') "| This is a 'free as in freedom'[*] software, distributed under the GPLv3.|"
       write(*,'(A75)') "| [*] https://www.gnu.org/philosophy/free-sw.en.html                      |"
       write(*,'(A75)') "+-------------------------------------------------------------------------+" 
       print*, ' '

       write(*, '(A, I5)') 'Number of coarray images = ', num_images()
    end if
  end subroutine welcome
  
  subroutine subtitle(text)
    !! Subroutine to print a subtitle.

    character(len = *), intent(in) :: text
    integer(k8) :: length
    character(len = 75) :: string2print

    length = len(text)
    
    string2print = '___________________________________________________________________________'
    if(this_image() == 1) write(*,'(A75)') string2print
    string2print(75 - length + 1 : 75) = text
    if(this_image() == 1) write(*,'(A75)') string2print
  end subroutine subtitle
end module misc

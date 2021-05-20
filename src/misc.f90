module misc
  !! Module containing miscellaneous math and numerics related functions and subroutines.

  use params, only: dp, k8, kB

  implicit none

  public
  private :: sort_int, sort_real

  interface sort
     module procedure :: sort_int, sort_real
  end interface sort
  
contains

  subroutine exit_with_message(message)
    !! Exit with error message.

    character(len = *), intent(in) :: message

    if(this_image() == 1) then
       print*, trim(message)
       stop
    end if
  end subroutine exit_with_message

  subroutine print_message(message)
    !! Print message.
    
    character(len = *), intent(in) :: message

    if(this_image() == 1) print*, trim(message)
  end subroutine print_message

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
  end subroutine write2file_rank2_real

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

    expi = cmplx(cos(x), sin(x))
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
    real(kind=8) :: aux, tmp
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
    
    integer(k8), intent(in) :: coarsemesh(3), q(3), refinement
    real(dp), intent(in) :: f(:)
    real(dp), intent(out) :: interpolation
    
    integer(k8) :: info, r0(3), r1(3), ipol, mode, count
    integer(k8), allocatable :: pivot(:)
    integer(k8) :: i000, i100, i010, i110, i001, i101, i011, i111, equalpol
    real(dp) :: x0, x1, y0, y1, z0, z1, x, y, z, v(2), v0(2), v1(2)
    real(dp), allocatable :: T(:, :), c(:)
    real(dp) :: aux

    aux = 0.0_dp

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
       x =  q(1)/dble(refinement*coarsemesh(1))
       y =  q(2)/dble(refinement*coarsemesh(2))
       z =  q(3)/dble(refinement*coarsemesh(3))

       !Coarse mesh walls
       x0 = floor(q(1)/dble(refinement))/dble(coarsemesh(1))
       y0 = floor(q(2)/dble(refinement))/dble(coarsemesh(2))
       z0 = floor(q(3)/dble(refinement))/dble(coarsemesh(3))
       x1 = ceiling(q(1)/dble(refinement))/dble(coarsemesh(1))
       y1 = ceiling(q(2)/dble(refinement))/dble(coarsemesh(2))
       z1 = ceiling(q(3)/dble(refinement))/dble(coarsemesh(3))

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
             v(count) = q(ipol)/dble(refinement*coarsemesh(ipol))
             v0(count) = floor(q(ipol)/dble(refinement))/dble(coarsemesh(ipol))
             v1(count) = ceiling(q(ipol)/dble(refinement))/dble(coarsemesh(ipol))
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
             x =  q(ipol)/dble(refinement*coarsemesh(ipol))
             x0 = floor(q(ipol)/dble(refinement))/dble(coarsemesh(ipol))
             x1 = ceiling(q(ipol)/dble(refinement))/dble(coarsemesh(ipol))

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
       write(*,'(A75)') "|   \      \         ___| |      | |  |.|__   ___ | |  / /_               |"
       write(*,'(A75)') "|    \    __\       / _ \ |   _  | |  |.'_ \ / _ \| | / __/               |"
       write(*,'(A75)') "|     \  \         |  __/ | |/ \_|/ \ |.|_) | (_) | |/ /_                 |"
       write(*,'(A75)') "|      \ \          \___|_|/|__/ |   ||..__/ \___/|_|\__/                 |"
       write(*,'(A75)') "|       \ \                /|                                             |"
       write(*,'(A75)') "|        \\                \|                                             |"
       write(*,'(A75)') "|         \\                                                              |"
       write(*,'(A75)') "|          \                                                              |"
       write(*,'(A75)') "|           \                                                             |"
       write(*,'(A75)') "| A solver for the coupled electron-phonon Boltzmann transport equations. |"
       write(*,'(A75)') "| Copyright (C) 2020- Nakib Haider Protik.                                |"
       write(*,'(A75)') "| This is a 'free as in freedom' software distributed under the GPLv3.    |"
       write(*,'(A75)') "+-------------------------------------------------------------------------+" 
       print*, ' '
       print*, 'Number of compute processes = ', num_images()
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

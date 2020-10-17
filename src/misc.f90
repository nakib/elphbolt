module misc
  !! Module containing miscellaneous math and numerics related functions and subroutines.

  use params, only: dp, k4

  implicit none

contains

  subroutine int_div(num,denom,q,r)
    !! Quotient(q) and remainder(r) of the integer division num/denom.

    integer(k4), intent(in) :: num, denom
    integer(k4), intent(out) :: q, r

    q = num/denom
    r = mod(num, denom)
  end subroutine int_div
  
  function cross_product(A,B)
    !! Cross product of A and B.

    real(dp), intent(in) :: A(3), B(3)
    real(dp) :: cross_product(3)

    cross_product(1) = A(2)*B(3) - A(3)*B(2)
    cross_product(2) = A(3)*B(1) - A(1)*B(3)
    cross_product(3) = A(1)*B(2) - A(2)*B(1)
  end function cross_product

  function mux_vector(v,mesh,base)
    !! Multiplex index of a single wave vector.
    !! v is the demultiplexed triplet of a wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    integer(k4), intent(in) :: v(3), mesh(3), base
    integer(k4) :: mux_vector 

    !TODO catch error
    !if(base < 0 .or. base > 1) then
    !   call exit_with_message("Base has to be either 0 or 1 in index.f90:mux_vector")
    !end if

    if(base == 0) then
       mux_vector = (v(3)*mesh(2) + v(2))*mesh(1) + v(1) + 1
    else
       mux_vector = ((v(3) - 1)*mesh(2) + (v(2) - 1))*mesh(1) + v(1)
    end if
  end function mux_vector

  subroutine demux_vector(i,v,mesh,base)
    !! Demultiplex index of a single wave vector.
    !! i is the multiplexed index of a wave vector (always 1-based).
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three.
    !! reciprocal lattice vectors.
    !! base chooses whether v has 0- or 1-based indexing.
    
    integer(k4), intent(in) :: i, mesh(3), base
    integer(k4), intent(out) :: v(3)
    integer(k4) :: aux

    !TODO Catch error
    !if(base < 0 .or. base > 1) then
    !   call exit_with_message("Base has to be either 0 or 1 in index.f90:demux_vector")
    !end if

    call int_div(i - 1, mesh(1), aux, v(1))
    call int_div(aux, mesh(2), v(3), v(2))
    if(base == 1) then
       v(1) = v(1) + 1
       v(2) = v(2) + 1
       v(3) = v(3) + 1 
    end if
  end subroutine demux_vector
  
  subroutine demux_mesh(index_mesh,nmesh,mesh,base,indexlist)
    !! Demultiplex all wave vector indices 
    !! (optionally, from a list of indices).
    !! Internally uses demux_vector.

    integer(k4), intent(in) :: nmesh, mesh(3), base
    integer(k4), optional, intent(in) :: indexlist(nmesh)
    integer(k4), intent(out) :: index_mesh(3, nmesh)
    integer(k4) :: i

    do i = 1, nmesh !over total number of wave vectors
       if(present(indexlist)) then
          call demux_vector(indexlist(i), index_mesh(:, i), mesh, base)
       else
          call demux_vector(i, index_mesh(:, i), mesh, base)
       end if
    end do
  end subroutine demux_mesh
end module misc

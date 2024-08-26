module vector_allreps_module

  use precision, only: i64, r64
  use misc, only: mux_vector, demux_vector, operator(.umklapp.)
  
  implicit none

  private
  public :: vector_allreps, &
       vector_allreps_add, vector_allreps_sub, &
       vector_allreps_change_grid, vector_allreps_print
  
  type vector_allreps
     !! A container for a (3-)vector and relevant arithmetic operations.
     !! It is convenient to use under certain circumstances, for example,
     !! when low level vector arithmetic is repetitive and error prone.
     !! However, the low level method will be more performant.

     !Grid independent representations:
     !Fractional coordinates (fractions of real/reciprocal lattice vectors)
     real(r64) :: frac(3) = 0.0
     !Cartesian coordinates
     real(r64) :: cart(3) = 0.0
     
     !Grid dependent representations:
     !0-based integer triplet with respect to a discretized grid
     integer(i64) :: int(3) = 0
     !Multipled index (1-based) of the integer triplet
     integer(i64) :: muxed_index = -1
  end type vector_allreps

  interface vector_allreps
     module procedure create
  end interface vector_allreps

contains
  
  function create(ind, grid, primitive_vecs) result(vector_obj)
    integer(i64), intent(in) :: ind, grid(3)
    real(r64), intent(in) :: primitive_vecs(3, 3)
    type(vector_allreps) :: vector_obj

    !vector_obj%grid = grid

    !vector_obj%primitive_vecs = primitive_vecs
    
    vector_obj%muxed_index = ind

    call demux_vector(ind, vector_obj%int, &
         grid, 0_i64)
    
    vector_obj%frac = dble(vector_obj%int)/grid
    
    vector_obj%cart = matmul(primitive_vecs, vector_obj%frac)
  end function create

  subroutine vector_allreps_print(v)
    !! Printer
    
    type(vector_allreps), intent(in) :: v

    print*, v%muxed_index
    print*, v%int
    print*, v%frac
    print*, v%cart
  end subroutine vector_allreps_print

  pure function vector_allreps_add(v1, v2, grid, primitive_vecs) result(v3)
    !! Adder
    
    type(vector_allreps), intent(in) :: v1
    type(vector_allreps), intent(in) :: v2
    integer(i64), intent(in) :: grid(3)
    real(r64), intent(in) :: primitive_vecs(3, 3)
    type(vector_allreps) :: v3

    !This arithmetic is independent of mesh density
    v3%frac = v1%frac .umklapp. v2%frac
    v3%cart = matmul(primitive_vecs, v3%frac)

    !This bit is depedent on the mesh density
    v3%int = nint(v3%frac*grid)
    v3%muxed_index = mux_vector(v3%int, grid, 0_i64)
  end function vector_allreps_add

  pure function vector_allreps_sub(v1, v2, grid, primitive_vecs) result(v3)
    !! Subtracter
    
    type(vector_allreps), intent(in) :: v1
    type(vector_allreps), intent(in) :: v2
    integer(i64), intent(in) :: grid(3)
    real(r64), intent(in) :: primitive_vecs(3, 3)
    type(vector_allreps) :: v3

    !This arithmetic is independent of mesh density
    v3%frac = v1%frac .umklapp. -v2%frac
    v3%cart = matmul(primitive_vecs, v3%frac)

    !This bit is depedent on the mesh density
    v3%int = nint(v3%frac*grid)
    v3%muxed_index = mux_vector(v3%int, grid, 0_i64)
  end function vector_allreps_sub

  pure function vector_allreps_change_grid(vin, grid) result(vout)
    !! Change grid.
    !! Beware, this is only well-defined if the vector is representable
    !! in both original and the new grids.

    type(vector_allreps), intent(in) :: vin
    integer(i64), intent(in) :: grid(3)
    type(vector_allreps) :: vout

    !Copy grid independent sector
    vout%frac = vin%frac
    vout%cart = vin%cart

    !This bit is depedent on the mesh density
    vout%int = nint(vin%frac*grid)
    vout%muxed_index = mux_vector(vout%int, grid, 0_i64)
  end function vector_allreps_change_grid
end module vector_allreps_module

module isotopes_module
  
  use precision, only: r64
  
  implicit none

  private
  public :: isotopes
  
  type isotopes
     integer :: numisotopes = 0
     real(r64), allocatable :: masses(:)
     real(r64), allocatable :: abundances(:)
   contains
     !TODO Add getters for:
     !VCA mass
     !VCA g-factor
     !DIB mass
     !DIB g-factor
  end type isotopes

  interface isotopes
     !! Create new collection of isotopes 
     module procedure :: constructor
  end interface isotopes
       
contains
  
  function constructor(masses, abundances) result(this)
    type(isotopes) :: this

    this%numisotopes = size(masses)
    allocate(this%masses(this%numisotopes), this%abundances(this%numisotopes))
    this%masses = masses
    this%abundances = abundances
  end function constructor
end module isotopes_module

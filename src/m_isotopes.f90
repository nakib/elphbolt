module m_isotopes
  
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
     module procedure :: constructor
  end interface isotopes
       
contains
  
  function constructor(masses, abundances) result(this)
    real(r64), intent(in) :: masses(:), abundances(:)
    type(isotopes) :: this

    this%numisotopes = size(masses)
    allocate(this%masses(this%numisotopes), this%abundances(this%numisotopes))
    this%masses = masses
    this%abundances = abundances
  end function constructor
end module m_isotopes

module isotopes_module
  
  use precision, only: r64
  
  implicit none

  private
  public :: isotopes
  
  type isotopes
     !! Datatype for storing information related to isotopes of a species
     integer :: numisotopes = 0
     !! Number of isotopes
     real(r64), allocatable :: masses(:)
     !! Mass of each isotope in atomic mass units
     real(r64), allocatable :: abundances(:)
     !! Abundance of each isoptope in percent
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
    real(r64), intent(in) :: masses(:)
    !! Array containing the mass of each isotope in atomic mass units
    real(r64), intent(in) :: abundances(:)
    !! Array containing the abundance of each isotope in percent
    type(isotopes) :: this

    this%numisotopes = size(masses)
    allocate(this%masses(this%numisotopes), this%abundances(this%numisotopes))
    this%masses = masses
    this%abundances = abundances
  end function constructor
end module isotopes_module

! Copyright 2023 elphbolt contributors.
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

module autodiff_m
  !! Module containing a simple implemention of the forward
  !! automatic differentiation. At the moment, only the real(64)
  !! kind is supported along with just the basis arithmetic operations.

  use params, only: r64

  implicit none

  private
  public :: autodiff, &
       operator(+), operator(-), operator(*), operator(/)

  type autodiff
     !! Augmented arithmetic of reals.

     real(r64) :: f
     real(r64) :: df
  end type autodiff

  interface operator(+)
     module procedure add
     module procedure add_scalar_right
     module procedure add_scalar_left
  end interface operator(+)

  interface operator(-)
     module procedure sub
     module procedure sub_scalar_right
     module procedure sub_scalar_left
     module procedure neg
  end interface operator(-)

  interface operator(*)
     module procedure mult
     module procedure mult_scalar_right
     module procedure mult_scalar_left
  end interface operator(*)

  interface operator(/)
     module procedure div
     module procedure div_scalar_right
     module procedure div_scalar_left
  end interface operator(/)

contains

  pure function add(x, y) result(z)
    type(autodiff), intent(in) :: x
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = x%f + y%f
    z%df = x%df + y%df
  end function add

  pure function add_scalar_right(x, s) result(z)
    type(autodiff), intent(in) :: x
    real(r64), intent(in) :: s
    type(autodiff) :: z

    z%f = x%f + s
    z%df = x%df
  end function add_scalar_right

  pure function add_scalar_left(s, y) result(z)
    real(r64), intent(in) :: s
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = s + y%f
    z%df = y%df
  end function add_scalar_left

  pure function sub(x, y) result(z)
    type(autodiff), intent(in) :: x
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z = add(x, autodiff(-y%f, -y%df))
  end function sub

  pure function sub_scalar_right(x, s) result(z)
    type(autodiff), intent(in) :: x
    real(r64), intent(in) :: s
    type(autodiff) :: z

    z = add_scalar_right(x, -1.0*s)
  end function sub_scalar_right

  pure function sub_scalar_left(s, y) result(z)
    real(r64), intent(in) :: s
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z = add_scalar_left(-1.0*s, y)
  end function sub_scalar_left

  pure function neg(x) result(y)
    type(autodiff), intent(in) :: x
    type(autodiff) :: y

    y = autodiff(-x%f, -x%df)
  end function neg

  pure function mult(x, y) result(z)
    type(autodiff), intent(in) :: x
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = x%f*y%f
    z%df = x%df*y%f + y%df*x%f
  end function mult

  pure function mult_scalar_right(x, s) result(z)
    type(autodiff), intent(in) :: x
    real(r64), intent(in) :: s
    type(autodiff) :: z

    z%f = s*x%f
    z%df = s*x%df
  end function mult_scalar_right

  pure function mult_scalar_left(s, y) result(z)
    real(r64), intent(in) :: s
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = s*y%f
    z%df = s*y%df
  end function mult_scalar_left

  pure function div(x, y) result(z)
    type(autodiff), intent(in) :: x
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = x%f/y%f
    z%df = (x%df*y%f - y%df*x%f)/y%f**2
  end function div

  pure function div_scalar_right(x, s) result(z)
    type(autodiff), intent(in) :: x
    real(r64), intent(in) :: s
    type(autodiff) :: z

    z%f = x%f/s
    z%df = x%df/s
  end function div_scalar_right

  pure function div_scalar_left(s, y) result(z)
    real(r64), intent(in) :: s
    type(autodiff), intent(in) :: y
    type(autodiff) :: z

    z%f = s/y%f
    z%df = (-s*y%df)/y%f**2
  end function div_scalar_left
end module autodiff_m

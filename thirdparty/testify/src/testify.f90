! Copyright (C) 2023 Nakib Haider Protik <nakib.haider.protik@gmail.com>
! This file is part of testify <https://github.com/nakib/testify>.
!
! testify is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! testify is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with testify. If not, see <http://www.gnu.org/licenses/>.
module testify_m
  
  use iso_fortran_env, only : &
       i8 => int8, i16 => int16, i32 => int32, i64 => int64, &
       r32 => real32, r64 => real64, r128 => real128
  
  implicit none

  private
  public testify
  
  type testify
     character(:), private, allocatable :: name
     logical, private :: status = .true.
     integer, private :: test_count = 0
     integer, private :: pass_count = 0
     
   contains
     procedure, public :: rename, report, get_status
     generic, public :: operator(+) => compose
     generic,  public :: assert => assert_scalar, assert_array
     procedure, private :: assert_scalar, assert_array, assert_elemental, &
          compose
  end type testify

  interface testify
     module procedure testify_name
     module procedure testify_array
  end interface testify
contains

  function testify_name(name) result(test_object)
    character(*), intent(in) :: name
    type(testify) :: test_object

    test_object%name = name
  end function testify_name

  function testify_array(array) result(summed)
    type(testify), intent(in) :: array(:)
    type(testify) :: summed
    integer :: itest

    summed = array(1)
    do itest = 2, size(array)
       summed = summed + array(itest)
    end do
  end function testify_array
  
  pure function compose(self, other) result(composed)
    class(testify), intent(in) :: self
    type(testify), intent(in) :: other
    type(testify) :: composed

    composed%name = self%name // '; ' // other%name
    composed%status = self%status .and. other%status
    composed%test_count = self%test_count + other%test_count
    composed%pass_count = self%pass_count + other%pass_count
  end function compose
  
  subroutine rename(self, name)
    class(testify), intent(inout) :: self
    character(*), intent(in) :: name

    self%name = name
  end subroutine rename

  pure logical function get_status(self)
    class(testify), intent(in) :: self

    get_status = self%status
  end function get_status
    
  subroutine assert_scalar(self, val, ref, tol)
    class(testify), intent(inout) :: self
    class(*), intent(in) :: val, ref
    class(*), intent(in), optional :: tol

    self%test_count = self%test_count + 1

    if(present(tol)) then
       self%status = self%assert_elemental(val, ref, tol) .eqv. .true.
    else
       self%status = self%assert_elemental(val, ref) .eqv. .true.
    end if

    if(self%status) self%pass_count = self%pass_count + 1

    call print_test_result(self%name, self%status)
  end subroutine assert_scalar
  
  subroutine assert_array(self, val, ref, tol)
    class(testify), intent(inout) :: self
    class(*), intent(in) :: val(:), ref(:)
    class(*), intent(in), optional :: tol
    
    self%test_count = self%test_count + 1

    if(present(tol)) then
       self%status = all(self%assert_elemental(val, ref, tol) .eqv. .true.)
    else
       self%status = all(self%assert_elemental(val, ref) .eqv. .true.)
    end if

    if(self%status) self%pass_count = self%pass_count + 1

    call print_test_result(self%name, self%status)
  end subroutine assert_array
  
  elemental logical function assert_elemental(self, val, ref, tol)
    class(testify), intent(in) :: self
    class(*), intent(in) :: val, ref
    class(*), intent(in), optional :: tol

    assert_elemental = .false.
    
    select type(val)
    type is(logical)
       select type(ref)
       type is(logical)
          assert_elemental = val .eqv. ref
       end select

    type is(character(*))
       select type(ref)
       type is(character(*))
          assert_elemental = val == ref
       end select

    type is(integer(i8))
       select type(ref)
       type is(integer(i8))
          assert_elemental = val == ref
       end select

    type is(integer(i16))
       select type(ref)
       type is(integer(i16))
          assert_elemental = val == ref
       end select

    type is(integer(i32))
       select type(ref)
       type is(integer(i32))
          assert_elemental = val == ref
       end select

    type is(integer(i64))
       select type(ref)
       type is(integer(i64))
          assert_elemental = val == ref
       end select

    type is(real(r32))
       select type(ref)
       type is(real(r32))
          if(present(tol)) then
             select type(tol)
             type is(real(r32))
                assert_elemental = abs(val - ref) < tol
             end select
          else
             assert_elemental = val == ref
          end if
       end select
       
    type is(real(r64))
       select type(ref)
       type is(real(r64))
          if(present(tol)) then
             select type(tol)
             type is(real(r64))
                assert_elemental = abs(val - ref) < tol
             end select
          else
             assert_elemental = val == ref
          end if
       end select

    type is(real(r128))
       select type(ref)
       type is(real(r128))
          if(present(tol)) then
             select type(tol)
             type is(real(r128))
                assert_elemental = abs(val - ref) < tol
             end select
          else
             assert_elemental = val == ref
          end if
       end select

    type is(complex(r32))
       select type(ref)
       type is(complex(r32))
          if(present(tol)) then
             select type(tol)
             type is(real(r32))
                assert_elemental = &
                     (abs(val%re - ref%re) < tol) .and. &
                     (abs(val%im - ref%im) < tol)
             end select
          else
             assert_elemental = val == ref
          end if
       end select

    type is(complex(r64))
       select type(ref)
       type is(complex(r64))
          if(present(tol)) then
             select type(tol)
             type is(real(r64))
                assert_elemental = &
                     (abs(val%re - ref%re) < tol) .and. &
                     (abs(val%im - ref%im) < tol)
             end select
          else
             assert_elemental = val == ref
          end if
       end select

    type is(complex(r128))
       select type(ref)
       type is(complex(r128))
          if(present(tol)) then
             select type(tol)
             type is(real(r128))
                assert_elemental = &
                     (abs(val%re - ref%re) < tol) .and. &
                     (abs(val%im - ref%im) < tol)
             end select
          else
             assert_elemental = val == ref
          end if
       end select
    end select
  end function assert_elemental

  subroutine print_test_result(name, passed)
    logical, intent(in) :: passed
    character(*), intent(in) :: name

    character(:), allocatable :: message_head, message_butt

    message_head = name // ' => ' // achar(27)
    message_butt = achar(27) // '[0m'
    if(passed) then
       print *, message_head // '[32m PASSED! :)' // message_butt
    else
       print *, message_head // '[31m FAILED! :(' // message_butt
    end if
  end subroutine print_test_result
  
  subroutine report(self)
    class(testify), intent(in) :: self

    print*, '+--------------------------------------------------+'
    print*, ' Tests carried out: [', self%name , ']'
    print*, ' Total number of tests: ', self%test_count
    print*, ' Number of tests passed: ', self%pass_count
    print*, '+--------------------------------------------------+'
  end subroutine report
end module testify_m

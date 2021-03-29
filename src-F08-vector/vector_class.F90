!!
!! VECTOR_CLASS
!!
!! A base class for defining an abstract mathematical vector that
!! defines an interface for basic BLAS-like vector space operations.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2021  Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! CLASS METHODS
!!
!!  CLASS(VECTOR) :: V
!!
!!  CALL V%CLONE(C [,N]) creates a clone C of V. The allocatable CLASS(VECTOR)
!!  variable C is allocated with the same dynamic type and internal structure
!!  as V, but the values of its vector elements are not defined (see COPY).
!!  C may be a rank-1 array, in which case N must also be specified, and its
!!  value is the size to which C is allocated.
!!
!!    This is a generic interface with deferred specific procedures CLONE1 and
!!    CLONE2 that subclasses must implement.
!!
!!  CALL V%COPY(SRC) copies the vector element values from the CLASS(VECTOR)
!!  variable SRC to V. V and SRC must have the same dynamic type and be
!!  compatibly defined, as is the case if V were a clone of SRC.
!!
!!    This is implemented using the non-virtual interface (NVI) pattern.
!!    Subclasses implement the private deferred procedure COPY_ and are
!!    assured that SRC and V will have the same dynamic type.
!!
!!  CALL V%SETVAL(VAL) sets the vector element values to the real scalar VAL.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  CALL V%SCALE(A) multiplies the vector elements of V by the real scalar A.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  CALL V%UPDATE(A, X [,B [,Y [,C]]]) performs SAXPY-like vector updates:
!!    V <-- A*X + V
!!    V <-- A*X + B*V
!!    V <-- A*X + B*Y + V
!!    V <-- A*X + B*Y + C*V
!!  V, X, and Y must all have the same dynamic type and be compatibly defined.
!!  A, B, and C are real scalars.
!!
!!    This is a generic procedure with specific procedures implemented using
!!    the NVI pattern. Subclasses implement the private deferred procedures
!!    UPDATE1_, UPDATE2_, UPDATE3_, UPDATE4_ corresponding to the different
!!    versions of the update, and are assured that V, X, and Y will all have
!!    the same dynamic type.
!!
!!  V%DOT(Y) is the dot product of V with the CLASS(VECTOR) variable Y. V and
!!  Y must have the same dynamic type and be compatibly defined.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  V%NORM2() is the l2 norm of V; i.e., SQRT(V%DOT(V))
!!
!!    This is a deferred procedure implemented by subclasses.
!!

module vector_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: vector
  contains
    generic :: clone => clone1, clone2
    procedure(clone1), deferred :: clone1 ! private except to subclass
    procedure(clone2), deferred :: clone2 ! private except to subclass
    procedure :: copy
    procedure(setval), deferred :: setval
    procedure(scale), deferred :: scale
    generic :: update => update1, update2, update3, update4
    procedure, private :: update1, update2, update3, update4
    procedure :: dot
    procedure(norm2), deferred :: norm2
    ! remaining are all private except to subclass
    procedure(copy), deferred :: copy_
    procedure(update1), deferred :: update1_
    procedure(update2), deferred :: update2_
    procedure(update3), deferred :: update3_
    procedure(update4), deferred :: update4_
    procedure(dot), deferred :: dot_
  end type

  abstract interface
    subroutine clone1(this, clone)
      import vector
      class(vector), intent(in)  :: this
      class(vector), allocatable, intent(out) :: clone
    end subroutine
    subroutine clone2(this, clone, n)
      import vector
      class(vector), intent(in)  :: this
      class(vector), allocatable, intent(out) :: clone(:)
      integer, intent(in) :: n
    end subroutine
  end interface

  abstract interface
    subroutine setval(this, val)
      import vector, r8
      class(vector), intent(inout) :: this
      real(r8), intent(in) :: val
    end subroutine
  end interface

  abstract interface
    subroutine scale(this, a)
      import vector, r8
      class(vector), intent(inout) :: this
      real(r8), intent(in) :: a
    end subroutine
  end interface

  abstract interface
    function norm2(this)
      import vector, r8
      class(vector), intent(in) :: this
      real(r8) :: norm2
    end function
  end interface

contains

  recursive subroutine copy(dest, src)
    class(vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    if (same_type_as(dest, src)) then
      call dest%copy_(src)
    else
      error stop 'incompatible arguments to VECTOR%COPY'
    end if
  end subroutine

  recursive function dot(x, y)
    class(vector), intent(in) :: x, y
    real(r8) :: dot
    if (same_type_as(x, y)) then
      dot = x%dot_(y)
    else
      error stop 'incompatible arguments to VECTOR%DOT'
    end if
  end function

  !! Conventional SAXPY procedure: y <-- y + a*x
  recursive subroutine update1(this, a, x)
    class(vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    if (a == 0.0_r8) return
    if (same_type_as(this, x)) then
      call this%update1_(a, x)
    else
      error stop 'incompatible arguments to VECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  recursive subroutine update2(this, a, x, b)
    class(vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    if (a == 0.0_r8) then
      call this%scale(b)
    else if (same_type_as(this, x)) then
      call this%update2_(a, x, b)
    else
      error stop 'incompatible arguments to VECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  recursive subroutine update3(this, a, x, b, y)
    class(vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    if (a == 0.0_r8) then
      call update1(this, b, y)
    else if (b == 0.0_r8) then
      call update1(this, a, x)
    else if (same_type_as(this, x) .and. same_type_as(this, y)) then
      call this%update3_(a, x, b, y)
    else
      error stop 'incompatible arguments to VECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  recursive subroutine update4(this, a, x, b, y, c)
    class(vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    if (a == 0.0_r8) then
      call update2(this, b, y, c)
    else if (b == 0.0_r8) then
      call update2(this, a, x, c)
    else if (same_type_as(this, x) .and. same_type_as(this, y)) then
      call this%update4_(a, x, b, y, c)
    else
      error stop 'incompatible arguments to VECTOR%UPDATE'
    end if
  end subroutine

end module vector_class

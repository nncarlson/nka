!!
!! GRID_VECTOR_TYPE
!!
!! An implementation of the VECTOR base class that stores cell-centered data
!! associated with a structured 2D grid. The data is stored naturally as a
!! rank-2 array. The array includes a layer of ghost cells. For convenience
!! all of the methods operate on this ghost data with the exception of the
!! reduction methods (dot product and 2-norm) which necessarily exclude these
!! from their computations.
!!
!! This is used by the nka_example.F90 program.
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

module grid_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  implicit none
  private

  type, extends(vector), public :: grid_vector
    integer :: nx, ny
    real(r8), allocatable :: array(:,:)
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setval
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dot_
    procedure :: norm2 => norm2_
    !! Additional procedures specific to this type
    generic :: init => init_dim, init_mold
    procedure, private :: init_dim, init_mold
  end type

contains

  !! Cell-centered data corresponding to a structured grid with NX zones in
  !! the x direction and NY zones in the y direction. The cell data are indexed
  !! as (i,j) with i=1..NX and j=1..NY. The array includes ghost cell elements
  !! with i = 0 or NX+1, and j=0 or NY+1.

  subroutine init_dim(this, nx, ny)
    class(grid_vector), intent(out) :: this
    integer, intent(in) :: nx, ny
    allocate(this%array(0:nx+1,0:ny+1))
    this%nx = nx
    this%ny = ny
  end subroutine

  subroutine init_mold(this, mold)
    class(grid_vector), intent(out) :: this
    class(grid_vector), intent(in)  :: mold
    allocate(this%array, source=mold%array)
  end subroutine

  subroutine clone1(this, clone)
    class(grid_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(clone, source=this)
  end subroutine

  subroutine clone2(this, clone, n)
    class(grid_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    allocate(clone(n), source=this)
  end subroutine

  subroutine copy_(dest, src)
    class(grid_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (grid_vector)
      dest%array(:,:) = src%array
    end select
  end subroutine

  subroutine setval(this, val)
    class(grid_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    this%array = val
  end subroutine

  subroutine scale(this, a)
    class(grid_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    this%array = a * this%array
  end subroutine

  !! Conventional SAXPY procedure: y <-- a*x + y
  subroutine update1_(this, a, x)
    class(grid_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (grid_vector)
      this%array = a * x%array + this%array
    end select
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  subroutine update2_(this, a, x, b)
    class(grid_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (grid_vector)
      this%array = a * x%array + b * this%array
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  subroutine update3_(this, a, x, b, y)
    class(grid_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (grid_vector)
      select type (y)
      class is (grid_vector)
        this%array = a * x%array + b * y%array + this%array
      end select
    end select
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  subroutine update4_(this, a, x, b, y, c)
    class(grid_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (grid_vector)
      select type (y)
      class is (grid_vector)
        this%array = a * x%array + b * y%array + c * this%array
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(grid_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    integer :: i, j
    select type (y)
    class is (grid_vector)
      dp = 0.0_r8
      do j = 1, x%ny
        do i = 1, x%nx
          dp = dp + x%array(i,j) * y%array(i,j)
        end do
      end do
    end select
  end function

  function norm2_(this)
    class(grid_vector), intent(in) :: this
    real(r8) :: norm2_
    integer :: i, j
    norm2_ = 0.0_r8
    do j = 1, this%ny
      do i = 1, this%nx
        norm2_ = norm2_ + this%array(i,j)**2
      end do
    end do
    norm2_ = sqrt(norm2_)
  end function

end module grid_vector_type

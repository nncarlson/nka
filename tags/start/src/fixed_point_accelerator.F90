!!
!! The FIXED_POINT_ACCELERATOR Module
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 15 Feb 2004; initial F90 version 1996.
!!
!! This module implements a nonlinear subspace acceleration procedure
!! of Carlson and Miller [1] for the the fixed-point iterative solution
!! of nonlinear equations.  A modified Newton iteration is one example
!! of such an iteration.
!!
!! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
!!     weighted moving finite element code I: in one dimension", SIAM J.
!!     Sci. Comput;, 19 (1998), pp. 728-765..
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 1996, 2004 Neil N. Carlson
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
!! PROGRAMING INTERFACE
!!
!! This module provides the derived data type FPA_STATE with private components
!! that encapsulates the entire state of the acceleration procedure, and the
!! following procedures.  All real arguments are of double precision kind.
!!
!!  CALL FPA_CREATE (THIS, F, MAXV, VTOL)
!!
!!    TYPE(FPA_STATE), INTENT(OUT) :: THIS
!!    REAL(KIND=DP),   INTENT(IN) :: F(:)
!!    INTEGER,         INTENT(IN) :: MAXV
!!    REAL(KIND=DP),   INTENT(IN), OPTIONAL :: VTOL
!!
!!    This creates a new state variable THIS capable of using as many as
!!    MVEC vectors in the acceleration procedure, whose size equals the
!!    size of the vector argument F.  The optional argument VTOL specifies
!!    the the vector drop tolerance: a vector is dropped when the sine of
!!    the angle between the vector the the subspace spanned by the preceding
!!    vectors is less than this value.  The default is 0.01.  The only use
!!    of F is to glean its shape; its value is ignored.
!!
!!  CALL FPA_DESTROY (THIS)
!!
!!    TYPE(FPA_STATE), INTENT(OUT) :: THIS
!!  
!!    This deallocates all the array components of the state THIS and
!!    returns it to its default initialization state.
!!
!!  CALL FPA_RESTART (THIS)
!!
!!    TYPE(FPA_STATE), INTENT(OUT) :: THIS
!!
!!    This sets a flag in the state THIS which will cause the next call
!!    to FPA_CORRECTION to begin accumulating a new subspace, flushing
!!    any previous subspace.
!!
!!  CALL FPA_CORRECTION (THIS, F)
!!
!!    TYPE(FPA_STATE), INTENT(INOUT) :: THIS
!!    REAL(KIND=DP),   INTENT(INOUT) :: F(:)
!!
!!    This call takes the function value F, which would be the correction
!!    vector in the fixed point iteration, and overwrites it with the
!!    accelerated correction computed from the current state THIS.  The
!!    state THIS is updated with input function value F and the returned
!!    correction.
!!
!! USAGE
!!
!!  The following simple example shows the usage of this acceleration
!!  procedure.  For more details, see the associated documentation.
!!  Consider a quasi-Newton iteration for solving the nonlinear system
!!  f(x) = 0.  Suppose pc(y) is some preconditioning procedure that applies
!!  some approximation to the inverse of the Jacobian of f(x) to the vector y.
!!  The original fixed point iteration would look something like
!!
!!    x = 0
!!    do <until converged>
!!      v = pc(f(x))
!!      x = x - v
!!    end do
!!
!!  The accelerated iteration would look something like
!!
!!    call fpa_create (fpa, v, mvec=5)
!!    x = 0
!!    do <until converged>
!!      v = pc(f(x))
!!      call fpa_correction (fpa, v)
!!      x = x - v
!!    end do
!!    call fpa_destroy (fpa)
!!
!! The create and destroy can of course be moved outside any nonlinear
!! solution procedure containing this iteration, and a single state variable
!! used for repeated calls to the procedure.  This avoids the repeated
!! allocations and deallocations of arrays associated with the state variable.
!! In this case, it might be advisable to include a call to fpa_restart
!! before the loop so that each iterative solve starts with clean slate.
!!

#include "f90_assert.fpp"

#ifdef SUPPORTS_TR15581
# define DEFINED allocated
#else
# define DEFINED associated
#endif

module fixed_point_accelerator

  implicit none
  private

  public :: fpa_create, fpa_destroy, fpa_restart, fpa_correction
  
  integer, parameter :: dp = kind(1.0d0)  ! double precision kind

#ifdef SUPPORTS_TR15581
  type, public :: fpa_state
    private
    logical :: initialized = .false.
    logical :: empty = .true.
    integer :: mvec = 0             ! maximum number of vectors
    real(kind=dp) :: tol = 0.01_dp  ! vector drop tolerance
    !! Subspace storage.
    real(kind=dp), allocatable :: v(:,:)   ! correction vectors
    real(kind=dp), allocatable :: w(:,:)   ! function difference vectors
    real(kind=dp), allocatable :: h(:,:)   ! matrix of inner products
    !! Linked-list organization of the vector storage.
    integer :: first, last, free
    integer, allocatable :: next(:), prev(:)
  end type fpa_state
#else
  type, public :: fpa_state
    private
    logical :: initialized = .false.
    logical :: empty = .true.
    integer :: mvec = 0             ! maximum number of vectors
    real(kind=dp) :: tol = 0.01_dp  ! vector drop tolerance
    !! Subspace storage.
    real(kind=dp), pointer :: v(:,:) => null()  ! correction vectors
    real(kind=dp), pointer :: w(:,:) => null()  ! function difference vectors
    real(kind=dp), pointer :: h(:,:) => null()  ! matrix of inner products
    !! Linked-list organization of the vector storage.
    integer :: first, last, free
    integer, pointer :: next(:) => null()
    integer, pointer :: prev(:) => null()
  end type fpa_state
#endif SUPPORTS_TR15581

contains

 !!
 !! FPA_CREATE
 !!

  subroutine fpa_create (this, f, maxv, vtol)

    type(fpa_state), intent(out) :: this
    real(kind=dp),   intent(in) :: f(:)
    integer,         intent(in) :: maxv
    real(kind=dp),   intent(in), optional :: vtol

    integer :: n, d1

    ASSERT( maxv > 0 )

    this%mvec = maxv
    n = this%mvec + 1
    d1 = size(f,dim=1)
    allocate (this%v(d1,n), this%w(d1,n))
    allocate (this%h(n,n), this%next(n), this%prev(n))

    if (present(vtol)) then
      ASSERT( vtol > 0.0_dp)
      this%tol = vtol
    end if

    this%initialized = .true.
    this%empty = .true.

  end subroutine fpa_create

 !!
 !! FPA_DESTROY
 !!

  subroutine fpa_destroy (this)

    type(fpa_state), intent(inout) :: this

    type(fpa_state) :: default_state
    
    ASSERT( this%initialized )
    deallocate (this%v, this%w, this%h, this%next, this%prev)
    this = default_state    ! Set default values

  end subroutine fpa_destroy

 !!
 !! FPA_RESTART
 !!

  subroutine fpa_restart (this)

    type(fpa_state), intent(inout) :: this

    ASSERT( this%initialized )

    this%empty = .true.

  end subroutine fpa_restart

 !!
 !! FPA_CORRECTION
 !!

  subroutine fpa_correction (this, f)

    type(fpa_state), intent(inout) :: this
    real(kind=dp),   intent(inout) :: f(:)

    ! local variables.
    integer :: i, j, k, new, nvec
    real(kind=dp) :: s, hkk, hkj, cj, c(this%mvec+1)

    ASSERT( this%initialized )
    ASSERT( size(f) == size(this%v,dim=1) )

    if (this%empty) then

     !!!
     !!! INITIAL VECTOR

      this%v(:,1) = f      ! Save the (unaccelerated) correction.
      this%w(:,1) = f      ! Save f to compute the difference on the next call.

      !! Initialize the vector linked list.
      this%first   = 1
      this%last    = 1
      this%next(1) = 0
      this%prev(1) = 0

      !! Initialize the free storage linked list.
      this%free = 2
      do k = 2, size(this%next)-1
        this%next(k) = k + 1
      end do
      this%next(size(this%next)) = 0

      this%empty = .false.

    else

     !!!
     !!! NEXT FUNCTION DIFFERENCE W

      this%w(:,this%first) = this%w(:,this%first) - f
      s = 1.0_dp / sqrt(dot_product(this%w(:,this%first), this%w(:,this%first)))

      !! Normalize w_1 and apply same factor to v_1.
      this%v(:,this%first) = s * this%v(:,this%first)
      this%w(:,this%first) = s * this%w(:,this%first)

      !! Update H.
      k = this%next(this%first)
      do while (k /= 0)
        this%h(this%first,k) = dot_product(this%w(:,this%first), this%w(:,k))
        k = this%next(k)
      end do

     !!!
     !!! CHOLESKI FACTORIZATION OF H

      this%h(this%first,this%first) = 1.0_dp
      k = this%next(this%first)
      nvec = 1

      do while (k /= 0)
        nvec = nvec + 1
        if (nvec > this%mvec) then  ! Maintain at most MVEC vectors:
          !! Drop the last vector and update the free storage list.
          ASSERT( this%last == k )
          this%next(this%last) = this%free
          this%free = k
          this%last = this%prev(k)
          this%next(this%last) = 0
          exit
        end if

        hkk = 1.0_dp           ! Single stage of Choleski factorization.
        j = this%first         ! Original matrix kept in lower triangle (unit diagonal).
        do while (j /= k)      ! Upper triangle holds the factorization.
          hkj = this%h(j,k)
          i = this%first
          do while (i /= j)
            hkj = hkj - this%h(k,i) * this%h(j,i)
            i = this%next(i)
          end do
          hkj = hkj / this%h(j,j)
          hkk = hkk - hkj**2
          this%h(k,j) = hkj
          j = this%next(j)
        end do

        if (hkk > this%tol**2) then
          this%h(k,k) = sqrt(hkk)
        else  ! The current w nearly lies in the span of the previous vectors.

          !! Drop this vector
          ASSERT( this%prev(k) /= 0 )
          this%next(this%prev(k)) = this%next(k)
          if (this%next(k) == 0) then
            this%last = this%prev(k)
          else
            this%prev(this%next(k)) = this%prev(k)
          end if

          this%next(k) = this%free    ! update the free storage list,
          this%free = k

          k = this%prev(k)            ! and back-up.
          nvec = nvec - 1

        endif
        k = this%next(k)
      end do

     !!!
     !!! PROJECT F ONTO THE SPAN OF THE W VECTORS.

      !! Forward substitution
      j = this%first
      do while (j /= 0)
        cj = dot_product(f, this%w(:,j))
        i = this%first
        do while (i /= j)
          cj = cj - this%h(j,i) * c(i)
          i = this%next(i)
        end do
        c(j) = cj / this%h(j,j)
        j = this%next(j)
      end do

      !! Backward substitution
      j = this%last
      do while (j /= 0)
        cj = c(j)
        i = this%last
        do while (i /= j)
          cj = cj - this%h(i,j) * c(i)
          i = this%prev(i)
        end do
        c(j) = cj / this%h(j,j)
        j = this%prev(j)
      end do

     !!!
     !!! ACCELERATED CORRECTION

      ASSERT( this%free /= 0 )
      new = this%free                    ! Locate storage for the new vectors.
      this%free = this%next(this%free)

      this%w(:,new) = f                  ! Save the original f for the next call.

      k = this%first                     ! Compute the accelerated correction,
      do while (k /= 0)
        f = f - c(k) * this%w(:,k) + c(k) * this%v(:,k)
        k = this%next(k)
      end do

      this%v(:,new) = f                  !   and save it for the next call.

      this%prev(new) = 0                 ! Prepend the vectors to the list.
      this%next(new) = this%first
      this%prev(this%first) = new
      this%first = new

    end if

  end subroutine fpa_correction

end module fixed_point_accelerator

!!
!! NLK_ACCEL
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! This module implements the nonlinear Krylov accelerator introduced in [1]
!! for inexact Newton's (IN) method, in which the correction equation of
!! Newton's method is only approximately solved because the Jacobian matrix
!! is approximated and/or the linear system is not solved exactly.  Placed
!! in the iteration loop, this black-box accelerator listens to the sequence
!! of inexact corrections and replaces them with accelerated corrections;
!! the resulting method is a type of accelerated inexact Newton (AIN) method.
!! Note that an IN iteration is merely a standard fixed point iteration for
!! a preconditioned system, and so this accelerator is more generally
!! applicable to fixed point iterations.
!!
!! This is a Fortran 2003 adaptation of the original Fortran 95 version that
!! implements the methods as type bound procedures.
!!
!! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
!!     weighted moving finite element code I: in one dimension", SIAM J.
!!     Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2010  Neil N. Carlson
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
!! This module defines the type NLK_ACCEL with private data components, which
!! encapsulates the state of the acceleration procedure, and the following
!! public type bound procedures.
!!
!!  CREATE (VLEN, MVEC) creates an object capable of using MVEC vectors in the
!!    acceleration procedure.  The length of the vectors is specified by VLEN.
!!
!!  DESTROY() deallocates all storage associated with the object.
!!
!!  SET_VEC_TOL (VTOL) sets the vector drop tolerance.  A vector is dropped
!!    from the acceleration subspace when the sine of the angle between the
!!    vector and the subspace spanned by the preceding vectors is less than
!!    this value.  If not set, the default value 0.01 is used.
!!
!!  SET_DOT_PROD (DOT_PROD) sets the procedure used to compute vector dot
!!    products.  DOT_PROD is a procedure pointer with the same interface
!!    as the intrinsic function DOT_PRODUCT with real rank-1 array arguments.
!!    If not set, the default is to use the instrinsic DOT_PRODUCT.  One case
!!    where the user is required to defined the dot product is in a parallel
!!    implementation of an iterative nonlinear solver where the vector
!!    components are distributed across processors.  Here a global dot
!!    product procedure that performs the required communication internally
!!    needs to be supplied.
!!
!!  CORRECTION (F) takes the function value F, which would be the correction
!!    vector in a fixed point iteration, and overwrites it with the accelerated
!!    correction computed from the acceleration subspace stored in the object.
!!    This subspace is updated prior to computing the correction using F and
!!    the previous function value and correction that were cached on the
!!    preceding call to CORRECTION, if any.
!!
!!  RELAX() -- A call to CORRECTION caches the passed function value and
!!    returned accelerated correction in the object for use on the following
!!    call to CORRECTION where they are used to update the acceleration
!!    subspace.  RELAX deletes these pending vectors from the object so that
!!    the next call to CORRECTION will not update the acceleration subspace
!!    prior to computing an accelerated correction.
!!
!!    This can be used to carry over the subspace from one nonlinear solve
!!    to another.  CORRECTION expects that the passed function value is
!!    connected to the preceding correction (if it exists), but this is not
!!    normally true for the first call in a subsequent nonlinear solve, and
!!    would result in the subspace being updated with bogus information.  A
!!    call to RELAX at the end of a nonlinear solve prevents this from
!!    occuring.
!!
!!  RESTART() flushes the acceleration subspace from the object, returning
!!    it to its state after the call to CREATE; the next call to CORRECTION
!!    will start the process of accumulating a new acceleration subspace.
!!    Typical usage would be to call RESTART at the start of each nonlinear
!!    solve in a series of solves, allowing the object to be reused and
!!    avoiding the overhead of repeated creation and destruction.
!!
!!  NUM_VEC() returns the number of vectors in the acceleration subspace.
!!
!!  MAX_VEC() returns the max number of vectors in the acceleration subspace.
!!
!!  VEC_LEN() returns the length of the vectors.
!!
!!  VEC_TOL() returns the vector drop tolerance.
!!
!!  REAL_KIND() returns the real kind parameter value of the object's
!!    real-valued data components.
!!
!!  DEFINED() returns the value true if the object is well-defined; otherwise
!!    it returns the value false.  Defined means that the data components of
!!    the object are properly and consistently defined.  Due to the amount of
!!    effort this function goes through to examine the object, it is is
!!    primarly intended to be used in debugging situations.
!!
!! USAGE
!!
!!  The following simple example shows the usage of this acceleration
!!  procedure.  For more details, see the associated documentation.
!!  Consider an inexact Newton iteration for solving the nonlinear system
!!  f(x) = 0.  Suppose pc(y) is some preconditioning procedure that applies
!!  some approximation to the inverse of the Jacobian of f(x) to the vector y.
!!  The original inexact-Newton iteration (equivalent to the fixed point
!!  iteration for pc(f(x)) = 0) would look something like
!!
!!    x = 0
!!    do <until converged>
!!      v = pc(f(x))
!!      x = x - v
!!    end do
!!
!!  The accelerated inexact Newton (AIN) iteration would look something like
!!
!!    call state%create (size(v), mvec=5)
!!    x = 0
!!    do <until converged>
!!      v = pc(f(x))
!!      call state%correction (v)
!!      x = x - v
!!    end do
!!    call state%destroy (state)
!!
!! The create and destroy can of course be moved outside any nonlinear
!! solution procedure containing this iteration, and a single state variable
!! used for repeated calls to the procedure.  This avoids the repeated
!! allocations and deallocations of arrays associated with the state variable.
!! In this case, one should either include a call to RESTART before the
!! loop so that each iterative solve starts with clean slate, or include a
!! call to NKA_RELAX after the loop so that first call to CORRECTION in
!! the next iterative solve doesn't update the acceleration subspace with
!! bogus information.
!!

#include "f90_assert.fpp"

module nlk_accel_type

  implicit none
  private

  integer, parameter :: r8 = selected_real_kind(15) ! 8-byte IEEE float

  type, public :: nlk_accel
    private
    logical :: subspace = .false.
    logical :: pending  = .false.
    integer :: vlen = 0         ! vector length
    integer :: mvec = 0         ! maximum number of vectors
    real(r8) :: vtol = 0.01_r8  ! vector drop tolerance
    procedure(dp), pointer, nopass :: dp => null()
    !! Subspace storage.
    real(r8), allocatable :: v(:,:)   ! correction vectors
    real(r8), allocatable :: w(:,:)   ! function difference vectors
    real(r8), allocatable :: h(:,:)   ! matrix of inner products
    !! Linked-list organization of the vector storage.
    integer :: first, last, free
    integer, allocatable :: next(:), prev(:)
  contains
    procedure :: create
    procedure :: destroy
    procedure :: set_vec_tol
    procedure :: set_dot_prod
    procedure :: vec_len
    procedure :: num_vec
    procedure :: max_vec
    procedure :: vec_tol
    procedure :: correction
    procedure :: relax
    procedure :: restart
    procedure :: defined
    procedure :: real_kind
  end type nlk_accel

  !abstract interface
  !  real(r8) function dp (x, y)
  !    import r8
  !    real(r8), intent(in) :: x(:), y(:)
  !  end function dp
  !end interface

contains

 !!
 !! NKA_CREATE
 !!

  subroutine create (state, vlen, mvec)

    class(nlk_accel), intent(out) :: state
    integer,         intent(in) :: vlen
    integer,         intent(in) :: mvec

    integer :: n

    ASSERT(mvec > 0)
    ASSERT(vlen >= 0)

    state%dp => dp

    state%vlen = vlen
    state%mvec = mvec
    n = mvec + 1
    allocate(state%v(vlen,n), state%w(vlen,n))
    allocate(state%h(n,n), state%next(n), state%prev(n))

    call restart (state)

    ASSERT(defined(state))

  end subroutine create

 !!
 !! NKA_DESTROY
 !!

  subroutine destroy (state)

    class(nlk_accel), intent(inout) :: state

    !type(nlk_accel) :: default_state

    if (allocated(state%v)) deallocate(state%v)
    if (allocated(state%w)) deallocate(state%w)
    if (allocated(state%h)) deallocate(state%h)
    if (allocated(state%next)) deallocate(state%next)
    if (allocated(state%prev)) deallocate(state%prev)

    !state = default_state    ! Set default values

  end subroutine destroy

  integer function num_vec (state)
    class(nlk_accel), intent(in) :: state
    integer :: k
    num_vec = 0
    k = state%first
    do while (k /= 0)
      num_vec = num_vec + 1
      k = state%next(k)
    end do
    if (state%pending) num_vec = num_vec - 1
  end function num_vec

  integer function max_vec (state)
    class(nlk_accel), intent(in) :: state
    max_vec = state%mvec
  end function max_vec

  integer function vec_len (state)
    class(nlk_accel), intent(in) :: state
    vec_len = state%vlen
  end function vec_len

  real(r8) function vec_tol (state)
    class(nlk_accel), intent(in) :: state
    vec_tol = state%vtol
  end function vec_tol

  subroutine set_vec_tol (state, vtol)
    class(nlk_accel), intent(inout) :: state
    real(r8), intent(in) :: vtol
    ASSERT(vtol > 0.0_r8)
    state%vtol = vtol
  end subroutine set_vec_tol

  subroutine set_dot_prod (state, dot_prod)
    class(nlk_accel), intent(inout) :: state
    procedure(dp), pointer :: dot_prod
    ASSERT(associated(dot_prod))
    state%dp => dot_prod
  end subroutine set_dot_prod

 !!
 !! NKA_CORRECTION
 !!

  subroutine correction (state, f)

    class(nlk_accel), intent(inout) :: state
    real(r8),   intent(inout) :: f(:)

    ! local variables.
    integer :: i, j, k, new, nvec
    real(r8) :: s, hkk, hkj, cj, c(state%mvec+1)

    ASSERT(defined(state))
    ASSERT(size(f) == size(state%v,dim=1))

   !!!
   !!! UPDATE THE ACCELERATION SUBSPACE

    if (state%pending) then

      !! Next function difference w_1.
      state%w(:,state%first) = state%w(:,state%first) - f
      s = sqrt(state%dp(state%w(:,state%first), state%w(:,state%first)))

      !! If the function difference is 0, we can't update the subspace with
      !! this data; so we toss it out and continue.  In this case it is likely
      !! that the outer iterative solution procedure has gone badly awry
      !! (unless the function value is itself 0), and we merely want to do
      !! something reasonable here and hope that situation is detected on the
      !! outside.
      if (s == 0.0_r8) call relax (state)

    end if

    if (state%pending) then

      !! Normalize w_1 and apply same factor to v_1.
      state%v(:,state%first) = state%v(:,state%first) / s
      state%w(:,state%first) = state%w(:,state%first) / s

      !! Update H.
      k = state%next(state%first)
      do while (k /= 0)
        state%h(state%first,k) = state%dp(state%w(:,state%first), state%w(:,k))
        k = state%next(k)
      end do

     !!!
     !!! CHOLESKI FACTORIZATION OF H

      state%h(state%first,state%first) = 1.0_r8
      k = state%next(state%first)
      nvec = 1

      do while (k /= 0)
        nvec = nvec + 1
        if (nvec > state%mvec) then  ! Maintain at most MVEC vectors:
          !! Drop the last vector and update the free storage list.
          ASSERT(state%last == k)
          state%next(state%last) = state%free
          state%free = k
          state%last = state%prev(k)
          state%next(state%last) = 0
          exit
        end if

        hkk = 1.0_r8           ! Single stage of Choleski factorization.
        j = state%first         ! Original matrix kept in lower triangle (unit diagonal).
        do while (j /= k)      ! Upper triangle holds the factorization.
          hkj = state%h(j,k)
          i = state%first
          do while (i /= j)
            hkj = hkj - state%h(k,i) * state%h(j,i)
            i = state%next(i)
          end do
          hkj = hkj / state%h(j,j)
          hkk = hkk - hkj**2
          state%h(k,j) = hkj
          j = state%next(j)
        end do

        if (hkk > state%vtol**2) then
          state%h(k,k) = sqrt(hkk)
        else  ! The current w nearly lies in the span of the previous vectors.

          !! Drop this vector
          ASSERT(state%prev(k) /= 0)
          state%next(state%prev(k)) = state%next(k)
          if (state%next(k) == 0) then
            state%last = state%prev(k)
          else
            state%prev(state%next(k)) = state%prev(k)
          end if

          state%next(k) = state%free    ! update the free storage list,
          state%free = k

          k = state%prev(k)            ! and back-up.
          nvec = nvec - 1

        end if
        k = state%next(k)
      end do

      ASSERT(state%first /= 0)
      state%subspace = .true.
      state%pending  = .false.

    end if

    !! Locate storage for the new vectors.
    ASSERT(state%free /= 0)
    new = state%free
    state%free = state%next(state%free)

    !! Save the original f for the next call.
    state%w(:,new) = f

   !!!
   !!! ACCELERATED CORRECTION

    if (state%subspace) then

      !! Project f onto the span of the w vectors: forward substitution
      j = state%first
      do while (j /= 0)
        cj = state%dp(f, state%w(:,j))
        i = state%first
        do while (i /= j)
          cj = cj - state%h(j,i) * c(i)
          i = state%next(i)
        end do
        c(j) = cj / state%h(j,j)
        j = state%next(j)
      end do

      !! Project f onto the span of the w vectors: backward substitution
      j = state%last
      do while (j /= 0)
        cj = c(j)
        i = state%last
        do while (i /= j)
          cj = cj - state%h(i,j) * c(i)
          i = state%prev(i)
        end do
        c(j) = cj / state%h(j,j)
        j = state%prev(j)
      end do

      !! The accelerated correction
      k = state%first
      do while (k /= 0)
        f = f - c(k) * state%w(:,k) + c(k) * state%v(:,k)
        k = state%next(k)
      end do

    end if

    !! Save the correction for the next call.
    state%v(:,new) = f

    !! Prepend the new vectors to the list.
    state%prev(new) = 0
    state%next(new) = state%first
    if (state%first == 0) then
      state%last = new
    else
      state%prev(state%first) = new
    end if
    state%first = new

    !! The original f and accelerated correction are cached for the next call.
    state%pending = .true.

  end subroutine correction

  real(r8) function dp (x, y)
    real(r8), intent(in) :: x(:), y(:)
    dp = dot_product(x, y)
  end function dp

 !!
 !! NKA_RESTART
 !!

  subroutine restart (state)

    class(nlk_accel), intent(inout) :: state

    integer :: k

    state%subspace = .false.
    state%pending  = .false.

    !! No vectors are stored.
    state%first = 0
    state%last  = 0

    !! Initialize the free storage linked list.
    state%free  = 1
    do k = 1, size(state%next)-1
      state%next(k) = k + 1
    end do
    state%next(size(state%next)) = 0

  end subroutine restart

 !!
 !! NKA_RELAX
 !!

  subroutine relax (state)

    class(nlk_accel), intent(inout) :: state

    integer :: new

    if (state%pending) then

      ASSERT(state%first /= 0)

      !! Drop the pending vectors.
      new = state%first
      state%first = state%next(state%first)
      if (state%first == 0) then
        state%last = 0
      else
        state%prev(state%first) = 0
      end if

      !! Update the free storage list.
      state%next(new) = state%free
      state%free = new

      state%pending = .false.

    end if

  end subroutine relax

 !!
 !! NKA_DEFINED
 !!

  logical function defined (state)

    class(nlk_accel), intent(in) :: state

    integer :: n
    logical, allocatable :: tag(:)

    CHECKLIST: do
      defined = .false.
      if (state%mvec < 1) exit
      if (.not.allocated(state%v)) exit
      if (.not.allocated(state%w)) exit
      if (any(shape(state%v) /= shape(state%w))) exit
      if (size(state%v,dim=1) /= state%vlen) exit
      if (size(state%v,dim=2) /= state%mvec+1) exit
      if (.not.allocated(state%h)) exit
      if (size(state%h,dim=1) /= state%mvec+1) exit
      if (size(state%h,dim=2) /= state%mvec+1) exit
      if (.not.allocated(state%next)) exit
      if (size(state%next) /= state%mvec+1) exit
      if (.not.allocated(state%prev)) exit
      if (size(state%prev) /= state%mvec+1) exit

      if (state%vtol <= 0.0_r8) exit

      n = size(state%next)
      if (any(state%next < 0) .or. any(state%next > n)) exit
      if (state%first < 0 .or. state%first > n) exit
      if (state%free  < 0 .or. state%free  > n) exit

      !! Tag array: each location is either in the free list or vector list.
      allocate(tag(size(state%next)))
      tag = .false.

      !! Check the vector list for consistency.
      if (state%first == 0) then
        if (state%last /= 0) exit
      else
        n = state%first
        if (state%prev(n) /= 0) exit
        tag(n) = .true.
        do while (state%next(n) /= 0)
          if (state%prev(state%next(n)) /= n) exit CHECKLIST
          n = state%next(n)
          if (tag(n)) exit CHECKLIST
          tag(n) = .true.
        end do
        if (state%last /= n) exit
      end if

      !! Check the free list.
      n = state%free
      do while (n /= 0)
        if (tag(n)) exit CHECKLIST
        tag(n) = .true.
        n = state%next(n)
      end do

      !! All locations accounted for?
      if (.not.all(tag)) exit

      defined = .true.
      exit
    end do CHECKLIST

    if (allocated(tag)) deallocate(tag)

  end function defined

  integer function real_kind (state)
    class(nlk_accel), intent(in) :: state
    real_kind = kind(state%h)
  end function real_kind

end module nlk_accel_type

!!
!! The FIXED_POINT_ACCELERATOR Module
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 15 Feb 2004; initial F90 version 1996.
!!
!! This module implements the nonlinear subspace acceleration method of Carlson
!! and Miller [1] for the fixed-point iterative solution of nonlinear equations.
!! An inexact Newton iteration, in which the Newton correction equation is only
!! approximately solved (because the Jacobian matrix is approximated and/or the
!! linear system is not solved exactly), can be interpreted as a fixed-point
!! iteration for a preconditioned equation.  Using this acceleration method in
!! such an iterative method results in a type of accelerated inexact Newton (AIN)
!! scheme.
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
!! following procedures that operate on variables of that type.  All real
!! arguments are of kind FPA_RK (double precision).
!!
!!  CALL FPA_CREATE (STATE, F, MAXV, VTOL)
!!
!!    TYPE(FPA_STATE), INTENT(OUT) :: STATE
!!    REAL(FPA_RK),    INTENT(IN) :: F(:)
!!    INTEGER,         INTENT(IN) :: MAXV
!!    REAL(FPA_RK),    INTENT(IN), OPTIONAL :: VTOL
!!
!!    This creates a new state variable STATE capable of using as many as
!!    MVEC vectors in the acceleration procedure, whose size equals the
!!    size of the vector argument F.  The optional argument VTOL specifies
!!    the the vector drop tolerance: a vector is dropped when the sine of
!!    the angle between the vector the the subspace spanned by the preceding
!!    vectors is less than this value.  The default is 0.01.  The only use
!!    of F is to glean its shape; its value is ignored.
!!
!!  CALL FPA_DESTROY (STATE)
!!
!!    TYPE(FPA_STATE), INTENT(INOUT) :: STATE
!!
!!    This deallocates all the array components of the state STATE and
!!    returns it to its default initialization state.
!!
!!  CALL FPA_CORRECTION (STATE, F, DP)
!!
!!    TYPE(FPA_STATE), INTENT(INOUT) :: STATE
!!    REAL(FPA_RK),    INTENT(INOUT) :: F(:)
!!    OPTIONAL :: DP
!!
!!    This call takes the function value F, which would be the correction
!!    vector in a fixed point iteration, and overwrites it with the
!!    accelerated correction computed from the subspace stored in STATE.
!!    This acceleration subspace is updated prior to computing the correction
!!    using F and previous function value and correction that were cached on
!!    the preceding call to FPA_CORRECTION (if any).  The input function
!!    value F and the returned correction are then cached in STATE for the
!!    next call.
!!
!!    DP is an optional procedure argument having the same interface as the
!!    intrinsic function DOT_PRODUCT.  If DP is present, it is used to compute
!!    vector dot products instead of DOT_PRODUCT.  For example, in a parallel
!!    implementation of an iterative nonlinear solve where the vector components
!!    are distributed across processors, a global dot product procedure, which
!!    performs the required communication internally, should be passed as DP
!!
!!  CALL FPA_RELAX (STATE)
!!
!!    TYPE(FPA_STATE), INTENT(INOUT) :: STATE
!!
!!    A call to FPA_CORRECTION caches the passed function value and returned
!!    accelerated correction in the state variable for use on the following
!!    call where they are used to update the acceleration subspace.  FPA_RELAX
!!    deletes these pending vectors from STATE, so that the next call to
!!    FPA_CORRECTION will not update the acceleration subspace prior to
!!    computing an accelerated correction.
!!
!!    This can be used to carry over the subspace from one nonlinear solve
!!    to another.  FPA_CORRECTION expects that the passed function value is
!!    connected to the preceding correction (if it exists), but this is not
!!    normally true for the first call in a subsequent nonlinear solve, and
!!    would result in the subspace being updated with bogus information.  A
!!    call to FPA_RELAX at the end of a nonlinear solve prevents this from
!!    occuring.  Be very cautious when using this strategy; it may not be
!!    appropriate.
!!
!!  CALL FPA_RESTART (STATE)
!!
!!    TYPE(FPA_STATE), INTENT(OUT) :: STATE
!!
!!    This call flushes the acceleration subspace from STATE, returning STATE
!!    to its state as returned by FPA_CREATE; the next call to FPA_CORRECTION
!!    will start the process of accumulating a new subspace.  Typical usage
!!    would be to call FPA_RESTART at the start of each in a series of
!!    nonlinear solves, allowing a single state to be reused and avoiding the
!!    overhead of repeated creation and destruction of state variables.
!!
!!  FPA_DEFINED(STATE)
!!
!!    TYPE(FPA_STATE), INTENT(IN) :: STATE
!!
!!    This logical function returns the value true if STATE is defined;
!!    otherwise it returns the value false.  Defined means that the components
!!    of STATE are properly defined, and in particular that the linked-list
!!    storage of vectors is coherently defined.  Due to the amount of effort
!!    this function goes through to examine the structure, it is primarily
!!    intended to be used in debugging situations.
!!
!! USAGE
!!
!!  The following simple example shows the usage of this acceleration
!!  procedure.  For more details, see the associated documentation.
!!  Consider an inexact-Newton iteration for solving the nonlinear system
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
!!  The accelerated inexact-Newton (AIN) iteration would look something like
!!
!!    call fpa_create (state, v, mvec=5)
!!    x = 0
!!    do <until converged>
!!      v = pc(f(x))
!!      call fpa_correction (state, v)
!!      x = x - v
!!    end do
!!    call fpa_destroy (state)
!!
!! The create and destroy can of course be moved outside any nonlinear
!! solution procedure containing this iteration, and a single state variable
!! used for repeated calls to the procedure.  This avoids the repeated
!! allocations and deallocations of arrays associated with the state variable.
!! In this case, one should either include a call to FPA_RESTART before the
!! loop so that each iterative solve starts with clean slate, or include a
!! call to FPA_RELAX after the loop so that first call to FPA_CORRECTION in
!! the next iterative solve doesn't update the acceleration subspace with
!! bogus information.
!!

#include "f90_assert.fpp"

#ifdef SUPPORTS_TR15581
# define _DEFINED_ allocated
#else
# define _DEFINED_ associated
#endif

module fixed_point_accelerator

  use fpa_kinds
  implicit none
  private

  public :: fpa_create, fpa_destroy, fpa_restart, fpa_correction, fpa_relax, fpa_defined, fpa_rk

  integer, parameter :: rk = fpa_rk

#ifdef SUPPORTS_TR15581
  type, public :: fpa_state
    private
    logical :: subspace = .false.
    logical :: pending  = .false.
    integer :: mvec = 0             ! maximum number of vectors
    real(kind=rk) :: tol = 0.01_rk  ! vector drop tolerance
    !! Subspace storage.
    real(kind=rk), allocatable :: v(:,:)   ! correction vectors
    real(kind=rk), allocatable :: w(:,:)   ! function difference vectors
    real(kind=rk), allocatable :: h(:,:)   ! matrix of inner products
    !! Linked-list organization of the vector storage.
    integer :: first, last, free
    integer, allocatable :: next(:), prev(:)
  end type fpa_state
#else
  type, public :: fpa_state
    private
    logical :: subspace = .false.
    logical :: pending  = .false.
    integer :: mvec = 0             ! maximum number of vectors
    real(kind=rk) :: tol = 0.01_rk  ! vector drop tolerance
    !! Subspace storage.
    real(kind=rk), pointer :: v(:,:) => null()  ! correction vectors
    real(kind=rk), pointer :: w(:,:) => null()  ! function difference vectors
    real(kind=rk), pointer :: h(:,:) => null()  ! matrix of inner products
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

  subroutine fpa_create (state, f, maxv, vtol)

    type(fpa_state), intent(out) :: state
    real(kind=rk),   intent(in) :: f(:)
    integer,         intent(in) :: maxv
    real(kind=rk),   intent(in), optional :: vtol

    integer :: n, d1

    ASSERT( maxv > 0 )

    state%mvec = maxv
    n = state%mvec + 1
    d1 = size(f,dim=1)
    allocate (state%v(d1,n), state%w(d1,n))
    allocate (state%h(n,n), state%next(n), state%prev(n))

    if (present(vtol)) then
      ASSERT( vtol > 0.0_rk )
      state%tol = vtol
    end if

    call fpa_restart (state)

    ASSERT( fpa_defined(state) )

  end subroutine fpa_create

 !!
 !! FPA_DESTROY
 !!

  subroutine fpa_destroy (state)

    type(fpa_state), intent(inout) :: state

    type(fpa_state) :: default_state

    if (_DEFINED_(state%v)) deallocate(state%v)
    if (_DEFINED_(state%w)) deallocate(state%w)
    if (_DEFINED_(state%h)) deallocate(state%h)
    if (_DEFINED_(state%next)) deallocate(state%next)
    if (_DEFINED_(state%prev)) deallocate(state%prev)

    state = default_state    ! Set default values

  end subroutine fpa_destroy

 !!
 !! FPA_CORRECTION
 !!

  subroutine fpa_correction (state, f, dp)

    type(fpa_state), intent(inout) :: state
    real(kind=rk),   intent(inout) :: f(:)

    !! Optional dot product procedure to use instead of the intrinsic DOT_PRODUCT.
    interface
      pure function dp (x, y)
        use fpa_kinds
        real(fpa_rk), intent(in) :: x(:), y(:)
        real(fpa_rk) :: dp
      end function dp
    end interface
    optional :: dp

    ! local variables.
    integer :: i, j, k, new, nvec
    real(kind=rk) :: s, hkk, hkj, cj, c(state%mvec+1)

    ASSERT( fpa_defined(state) )
    ASSERT( size(f) == size(state%v,dim=1) )

   !!!
   !!! UPDATE THE ACCELERATION SUBSPACE

    if (state%pending) then

      !! Next function difference w_1.
      state%w(:,state%first) = state%w(:,state%first) - f
      if (present(dp)) then
        s = sqrt(dp(state%w(:,state%first), state%w(:,state%first)))
      else
        s = sqrt(dot_product(state%w(:,state%first), state%w(:,state%first)))
      end if

      !! If the function difference is 0, we can't update the subspace with
      !! this data; so we toss it out and continue.  In this case it is likely
      !! that the outer iterative solution procedure has gone badly awry
      !! (unless the function value is itself 0), and we merely want to do
      !! something reasonable here and hope that situation is detected on the
      !! outside.
      if (s == 0.0_rk) call fpa_relax (state)

    end if

    if (state%pending) then

      !! Normalize w_1 and apply same factor to v_1.
      state%v(:,state%first) = state%v(:,state%first) / s
      state%w(:,state%first) = state%w(:,state%first) / s

      !! Update H.
      k = state%next(state%first)
      do while (k /= 0)
        if (present(dp)) then
          state%h(state%first,k) = dp(state%w(:,state%first), state%w(:,k))
        else
          state%h(state%first,k) = dot_product(state%w(:,state%first), state%w(:,k))
        end if
        k = state%next(k)
      end do

     !!!
     !!! CHOLESKI FACTORIZATION OF H

      state%h(state%first,state%first) = 1.0_rk
      k = state%next(state%first)
      nvec = 1

      do while (k /= 0)
        nvec = nvec + 1
        if (nvec > state%mvec) then  ! Maintain at most MVEC vectors:
          !! Drop the last vector and update the free storage list.
          ASSERT( state%last == k )
          state%next(state%last) = state%free
          state%free = k
          state%last = state%prev(k)
          state%next(state%last) = 0
          exit
        end if

        hkk = 1.0_rk           ! Single stage of Choleski factorization.
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

        if (hkk > state%tol**2) then
          state%h(k,k) = sqrt(hkk)
        else  ! The current w nearly lies in the span of the previous vectors.

          !! Drop this vector
          ASSERT( state%prev(k) /= 0 )
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

      ASSERT( state%first /= 0 )
      state%subspace = .true.

    end if

    !! Locate storage for the new vectors.
    ASSERT( state%free /= 0 )
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
        if (present(dp)) then
          cj = dp(f, state%w(:,j))
        else
          cj = dot_product(f, state%w(:,j))
        endif
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

  end subroutine fpa_correction

 !!
 !! FPA_RESTART
 !!

  subroutine fpa_restart (state)

    type(fpa_state), intent(inout) :: state

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

  end subroutine fpa_restart

 !!
 !! FPA_RELAX
 !!

  subroutine fpa_relax (state)

    type(fpa_state), intent(inout) :: state

    integer :: new

    if (state%pending) then

      ASSERT( state%first /= 0 )

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

  end subroutine fpa_relax

 !!
 !! FPA_DEFINED
 !!

  logical function fpa_defined (state)

    type(fpa_state), intent(in) :: state

    integer :: n
    logical, allocatable :: tag(:)

    CHECKLIST: do
      fpa_defined = .false.
      if (state%mvec < 1) exit
      if (.not._DEFINED_(state%v)) exit
      if (.not._DEFINED_(state%w)) exit
      if (any(shape(state%v) /= shape(state%w))) exit
      if (size(state%v,dim=2) /= state%mvec+1) exit
      if (.not._DEFINED_(state%h)) exit
      if (size(state%h,dim=1) /= state%mvec+1) exit
      if (size(state%h,dim=2) /= state%mvec+1) exit
      if (.not._DEFINED_(state%next)) exit
      if (size(state%next) /= state%mvec+1) exit
      if (.not._DEFINED_(state%prev)) exit
      if (size(state%prev) /= state%mvec+1) exit

      if (state%tol <= 0.0_rk) exit

      n = size(state%next)
      if (any(state%next < 0) .or. any(state%next > n)) exit
      if (any(state%prev < 0) .or. any(state%prev > n)) exit
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

      fpa_defined = .true.
      exit
    end do CHECKLIST

    if (allocated(tag)) deallocate(tag)

  end function fpa_defined

end module fixed_point_accelerator

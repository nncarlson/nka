!!
!!  NKA_EXAMPLE
!!
!!  Example/test program for the NKA_TYPE module.
!!
!!  Neil N. Carlson <neil.n.carlson@gmail.com> 15 Feb 2004
!!  Last revised 4 Jul 2009
!!
!!  This program exercises the nonlinear Krylov acceleration procedure
!!  implemented by the NONLINEAR_KRYLOV_ACCELERATOR module by applying it to
!!  the solution of the nonlinear elliptic problem
!!
!!      -Div[(a + u)Grad[u]] = q,  u:[0,1]^2 -> R
!!
!!  with zero boundary values.  We use a mimetic discretization over a
!!  completely regular rectangular mesh.  The scalar field u is discretized
!!  in the cell-based space.  The discrete equations are of the form r(u) = 0,
!!  where r(u) = A(u) u - q.  A(u) is a symmetric positive definite matrix
!!  that depends on the value of the scalar field u.  To solve r(u) = 0 we
!!  use a fixed point iteration applied to a preconditioned r(u).  The
!!  preconditioner is one or more passes of SSOR to approximately solve
!!  A(u)^{-1} r(u) (here u is fixed).  We use a=0.02 and a uniform unit
!!  source field q.
!!
!!  We first solve the equation with an accelerated iteration, and then,
!!  for comparison, solve again without acceleration.  In each case, the
!!  l2 norm of r(u), the reduction factor in the norm, and the net convergence
!!  rate are printed for each iteration.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2004, 2009  Neil N. Carlson
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

program nka_example

  use nka_type
  implicit none

  integer,  parameter :: dp = kind(1.0d0)
  integer,  parameter :: N = 50
  real(dp), parameter :: A = 0.02_dp

  type :: system
    integer :: nx, ny
    real(dp) :: a, hx, hy
    real(dp), pointer :: ax(:,:) => null(), ay(:,:) => null(), ac(:,:) => null(), q(:,:) => null()
    type(nka) :: nka
  end type system

  type(system) :: sys
  real(dp) :: u(N*N), omega
  integer :: nsweep

  !! Initialize the system
  sys%a = A
  sys%nx = N
  sys%ny = N
  sys%hx = 1.0_dp / sys%nx
  sys%hy = 1.0_dp / sys%ny
  allocate(sys%ax(sys%nx+1,sys%ny), sys%ay(sys%nx,sys%ny+1), sys%ac(sys%nx,sys%ny))
  allocate(sys%q(sys%nx,sys%ny))
  call nka_init (sys%nka, size(u), mvec=5)
  sys%q = sys%hx*sys%hy

  nsweep = 2
  omega = 1.4_dp

  call solve (sys, nsweep, omega, use_nka=.true., sol=u)
  call solve (sys, nsweep, omega, use_nka=.false., sol=u)

  call nka_delete (sys%nka)

contains

  subroutine solve (sys, nsweep, omega, use_nka, sol)

    type(system), intent(inout) :: sys
    integer,  intent(in)  :: nsweep
    real(dp), intent(in)  :: omega
    logical,  intent(in)  :: use_nka
    real(dp), intent(out) :: sol(:)

    integer :: itr
    real(dp) :: rnorm, rnorm0, red, rate
    real(dp), target :: upad(0:sys%nx+1,0:sys%ny+1)
    real(dp), pointer :: u(:,:)
    real(dp) :: r(sys%nx*sys%ny)
    integer, parameter :: MAXITR = 999
    real(dp), parameter :: TOL = 1.0d-6

    if (use_nka) then
      write(unit=*,fmt='(/,a,/)') 'ACCELERATED SOLVE'
    else
      write(unit=*,fmt='(/,a,/)') 'UNACCELERATED SOLVE'
    end if
    write(unit=*,fmt='(a4,a14,a13,a8)') 'Iter', 'Residual Norm', 'Reduction', 'Rate'

    upad = 0.0_dp  ! zero initial solution and zero value for ghost cells
    u => upad(1:sys%nx,1:sys%ny)  ! the true part of the solution vector

    call residual (sys, upad, r)
    rnorm0 = l2norm(r)
    write(*,'(i3,a,es14.6)') 0, ':', rnorm0

    do itr = 1, MAXITR
      call ssor_pc (sys, nsweep, omega, r)
      if (use_nka) call nka_accel_update (sys%nka, r)
      u = u - reshape(r,shape(u))

      call residual (sys, upad, r)
      rnorm = l2norm(r)
      red = rnorm / rnorm0
      rate = red**(1.0_dp/itr)
      write(unit=*,fmt='(i3,a,es14.6,es13.3,f8.3)') itr, ':', rnorm, red, rate

      if (rnorm < TOL * rnorm0) exit
    end do

    sol = reshape(u, shape(sol))

  end subroutine solve

  subroutine residual (sys, u, r)

    type(system), intent(inout) :: sys
    real(dp), intent(in) :: u(0:,0:)  ! solution array padded with ghost cells
    real(dp), intent(out) :: r(sys%nx,sys%ny)

    integer :: j, k

    call update_system (sys, u)

    do k = 1, sys%ny
      do j = 1, sys%nx
        r(j,k) = sys%ac(j,k) * u(j,k) - sys%ax(j,k) * u(j-1,k) - sys%ax(j+1,k) * u(j+1,k) &
                                      - sys%ay(j,k) * u(j,k-1) - sys%ay(j,k+1) * u(j,k+1) &
                                      - sys%q(j,k)
      end do
    end do

  end subroutine residual

  subroutine update_system (sys, upad)

    type(system), intent(inout) :: sys
    real(dp), intent(in) :: upad(0:,0:) ! solution vector padded with ghost cells

    integer :: j, k
    real(dp) :: t, rx, ry

    rx = sys%hx / sys%hy
    ry = sys%hy / sys%hx

    sys%ax = 0.0_dp
    sys%ay = 0.0_dp
    do k = 1, sys%ny
      do j = 1, sys%nx
        t = 1.0_dp / (sys%a + upad(j,k))
        sys%ax(j,k)   = sys%ax(j,k)   + rx*t
        sys%ax(j+1,k) = sys%ax(j+1,k) + rx*t
        sys%ay(j,k)   = sys%ay(j,k)   + ry*t
        sys%ay(j,k+1) = sys%ay(j,k+1) + ry*t
      end do
    end do

    sys%ax = 2.0_dp / sys%ax
    sys%ay = 2.0_dp / sys%ay

    do k = 1, sys%ny
      do j = 1, sys%nx
        sys%ac(j,k) = sys%ax(j,k) + sys%ax(j+1,k) + sys%ay(j,k) + sys%ay(j,k+1)
      end do
    end do

  end subroutine update_system

  subroutine ssor_pc (sys, nsweep, omega, r)

    type(system), intent(in)    :: sys
    integer,      intent(in)    :: nsweep
    real(dp),     intent(in)    :: omega
    real(dp),     intent(inout) :: r(sys%nx,sys%ny)

    integer :: i, j, k
    real(dp) :: z(0:sys%nx+1,0:sys%ny+1)

    z = 0.0_dp

    do i = 1, nsweep
      !! Forward sweep
      do k = 1, sys%ny
        do j = 1, sys%nx
          z(j,k) = (1 - omega)*z(j,k) + omega*(r(j,k) &
                     + sys%ax(j,k)*z(j-1,k) + sys%ax(j+1,k)*z(j+1,k)  &
                     + sys%ay(j,k)*z(j,k-1) + sys%ay(j,k+1)*z(j,k+1))/sys%ac(j,k)
        end do
      end do
      !! Backward sweep
      do k = sys%ny, 1, -1
        do j = sys%nx, 1, -1
          z(j,k) = (1 - omega)*z(j,k) + omega*(r(j,k) &
                     + sys%ax(j,k)*z(j-1,k) + sys%ax(j+1,k)*z(j+1,k)  &
                     + sys%ay(j,k)*z(j,k-1) + sys%ay(j,k+1)*z(j,k+1))/sys%ac(j,k)
        end do
      end do
    end do

    !! Copy result to return array.
    r = z(1:sys%nx,1:sys%ny)

  end subroutine ssor_pc

  real(dp) function l2norm (x)
    real(dp), intent(in) :: x(:)
    l2norm = sqrt(sum(x**2))
  end function l2norm

end program nka_example

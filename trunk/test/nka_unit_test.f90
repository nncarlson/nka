!!
!!  NKA_UNIT_TEST
!!
!!  Unit test program for the NONLINEAR_KRYLOV_ACCELERATOR module.
!!
!!  Neil N. Carlson <nnc@newmexico.com> 15 Feb 2004
!!  Last revised 14 Oct 2006
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

program nka_unit_test

use nonlinear_krylov_accelerator
implicit none

integer, parameter :: dp = kind(1.0d0)
integer, parameter :: N = 50
real(kind=dp), parameter :: A = 0.02_dp

type :: system
  integer :: nx, ny
  type(nka_state) :: fpa
  real(kind=dp) :: a
  real(kind=dp), pointer :: ax(:,:) => null(), ay(:,:) => null(), ac(:,:) => null(), q(:,:) => null()
end type system

type(system) :: sys
real(kind=dp) :: u(N*N), h

!! Initialize the system
sys%nx = N
sys%ny = N
sys%a = A
allocate(sys%ax(sys%nx+1,sys%ny), sys%ay(sys%nx,sys%ny+1), sys%ac(sys%nx,sys%ny), sys%q(sys%nx,sys%ny))
call nka_create (sys%fpa, size(u), maxv=5)
h = 1.0_dp / N
sys%q = h**2

call solve (sys, u, nsweep=2, omega=1.4_dp, use_fpa=.true.)
call solve (sys, u, nsweep=2, omega=1.4_dp, use_fpa=.false.)

call nka_destroy(sys%fpa)

contains

  subroutine solve (sys, sol, nsweep, omega, use_fpa)

    type(system), intent(inout) :: sys
    real(kind=dp), intent(out) :: sol(:)
    integer, intent(in) :: nsweep
    real(kind=dp), intent(in) :: omega
    logical, intent(in) :: use_fpa

    integer :: itr
    real(kind=dp) :: rnorm, rnorm0, red, rate, rflat(sys%nx*sys%ny)
    real(kind=dp), dimension(0:sys%nx+1,0:sys%ny+1), target :: upad, rpad
    real(kind=dp), dimension(:,:), pointer :: u, r

    if (use_fpa) then
      write(unit=*,fmt='(/,a,/)') 'ACCELERATED SOLVE'
    else
      write(unit=*,fmt='(/,a,/)') 'UNACCELERATED SOLVE'
    end if
    write(unit=*,fmt='(a4,a14,a13,a8)') 'Iter', 'Residual Norm', 'Reduction', 'Rate'

    u => upad(1:sys%nx,1:sys%ny)
    r => rpad(1:sys%nx,1:sys%ny)

    upad = 0.0_dp  ! We start with a zero initial solution, and for the padded vectors,
    rpad = 0.0_dp  ! we need to zero the values in the surrounding ghost cells.

    r = residual(sys, upad)
    rnorm0 = norm(r)
    write(unit=*,fmt='(i3,a,es14.6)') 0, ':', rnorm0

    do itr = 1, 999
      rpad = ssor(sys, nsweep, omega, rpad)
      if (use_fpa) then
        rflat = reshape(r, shape(rflat))
        call nka_correction (sys%fpa, rflat)
        r = reshape(rflat, shape(r))
      end if
      u = u - r

      r = residual(sys, upad)
      rnorm = norm(r)
      red = rnorm / rnorm0
      rate = red**(1.0_dp/itr)
      write(unit=*,fmt='(i3,a,es14.6,es13.3,f8.3)') itr, ':', rnorm, red, rate

      if (rnorm < 1.0e-6 * rnorm0) exit
    end do

    sol = reshape(u, shape(sol))

  end subroutine solve
  
  function residual (sys, u) result (r)
  
    type(system),  intent(inout) :: sys
    real(kind=dp), intent(in)    :: u(0:,0:)
    real(kind=dp) :: r(sys%nx,sys%ny)
    
    integer :: j, k
    
    call update_system (sys, u)
    
    do k = 1, sys%ny
      do j = 1, sys%nx
        r(j,k) = sys%ac(j,k) * u(j,k) - sys%ax(j,k) * u(j-1,k) - sys%ax(j+1,k) * u(j+1,k) &
                                      - sys%ay(j,k) * u(j,k-1) - sys%ay(j,k+1) * u(j,k+1) &
                                      - sys%q(j,k)
      end do
    end do
    
  end function residual
  
  subroutine update_system (sys, u)
  
    type(system),  intent(inout) :: sys
    real(kind=dp), intent(in)    :: u(0:,0:)
  
    integer :: j, k
    real(kind=dp) :: t
    
    sys%ax = 0.0_dp
    sys%ay = 0.0_dp
    do k = 1, sys%ny
      do j = 1, sys%nx
        t = 1.0_dp / (sys%a + u(j,k))
        sys%ax(j,k)   = sys%ax(j,k)   + t
        sys%ax(j+1,k) = sys%ax(j+1,k) + t
        sys%ay(j,k)   = sys%ay(j,k)   + t
        sys%ay(j,k+1) = sys%ay(j,k+1) + t
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
  
  function ssor (sys, nsweep, omega, r) result (z)
  
    type(system),  intent(in) :: sys
    integer,       intent(in) :: nsweep
    real(kind=dp), intent(in) :: omega
    real(kind=dp), intent(in) :: r(0:,0:)
    real(kind=dp) :: z(0:sys%nx+1,0:sys%ny+1)
    
    integer :: i, j, k
    
    z = r
    
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
    
  end function ssor
    
  function norm (x) result (normx)
    real(kind=dp), intent(in) :: x(:,:)
    real(kind=dp) :: normx
    normx = sqrt(sum(x**2))
  end function norm

end program nka_unit_test

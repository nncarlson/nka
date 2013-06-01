!!
!!  NKA_EXAMPLE
!!
!!  An example/test program for the NKA_TYPE module.
!!
!!  Neil N. Carlson <neil.n.carlson@gmail.com> 15 Feb 2004
!!  Last revised 31 May 2013
!!
!!  This program exercises the nonlinear Krylov acceleration procedure
!!  implemented by the NKA_TYPE module by applying it to the solution of
!!  the nonlinear elliptic problem
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
!!  A(u)^{-1} r(u) (here u is fixed).  We use a uniform unit source field q.
!!
!!  Several options can be specified on the command line (use --help to get
!!  a full list and their default values):
!!
!!    -a A         The equation parameter A.
!!    -n N         Number of grid cells in each coordinate direction.
!!    --sweeps S   Use S SSOR sweeps as the preconditioner.
!!    --omega W    SSOR over-relaxation parameter W.
!!    --nka-vec M  Use an NKA acceleration subspace of dimension at most M;
!!                 0 turns off NKA acceleration.
!!
!!  The program writes the l2 norm of r(u), the reduction factor in the norm,
!!  and the net convergence rate for each iteration.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2004, 2009, 2013  Neil N. Carlson
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SYSTEM_TYPE -- Defines the discrete nonlinear system for the PDE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module system_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: system
    private
    integer :: nx, ny
    real(r8) :: a, hx, hy
    real(r8), allocatable :: ax(:,:), ay(:,:), ac(:,:), q(:,:)
  contains
    procedure :: init => system_init
    procedure :: residual => system_residual
    procedure :: pc_ssor => system_pc_ssor
    procedure :: size => system_size
    procedure :: xsize => system_xsize
    procedure :: ysize => system_ysize
  end type system

contains

  subroutine system_init (this, a, nx, ny)
    class(system), intent(out) :: this
    real(r8), intent(in) :: a
    integer,  intent(in) :: nx, ny
    ASSERT(a > 0.0_r8)
    ASSERT(nx >= 3)
    ASSERT(ny >= 3)
    this%a  = a
    this%nx = nx
    this%ny = ny
    this%hx = 1.0_r8 / nx
    this%hy = 1.0_r8 / ny
    allocate(this%ax(nx+1,ny), this%ay(nx,ny+1), this%ac(nx,ny))
    allocate(this%q(nx,ny))
    this%q = 1.0_r8
  end subroutine system_init

  pure integer function system_size (this)
    class(system), intent(in) :: this
    system_size = this%nx * this%ny
  end function system_size

  pure integer function system_xsize (this)
    class(system), intent(in) :: this
    system_xsize = this%nx
  end function system_xsize

  pure integer function system_ysize (this)
    class(system), intent(in) :: this
    system_ysize = this%ny
  end function system_ysize

  subroutine system_residual (this, uext, r)
    class(system), intent(inout) :: this
    real(r8), intent(in)  :: uext(0:,0:)  ! solution array extended with boundary data
    real(r8), intent(out) :: r(:,:)
    integer :: j, k
    ASSERT(ubound(uext,1) == this%nx+1)
    ASSERT(ubound(uext,2) == this%ny+1)
    ASSERT(size(r,1) == this%nx)
    ASSERT(size(r,2) == this%ny)
    call update_system (this, uext)
    do k = 1, this%ny
      do j = 1, this%nx
        r(j,k) = this%ac(j,k)*uext(j,k) - this%ax(j,k)*uext(j-1,k) - this%ax(j+1,k)*uext(j+1,k) &
                                        - this%ay(j,k)*uext(j,k-1) - this%ay(j,k+1)*uext(j,k+1) &
                                        - this%q(j,k)
      end do
    end do
  end subroutine system_residual

  subroutine update_system (this, uext)
    class(system), intent(inout) :: this
    real(r8), intent(in) :: uext(0:,0:) ! solution array extended with boundary data
    integer :: j, k
    real(r8) :: t
    this%ax = 0.0_r8
    this%ay = 0.0_r8
    do k = 1, this%ny
      do j = 1, this%nx
        t = 1.0_r8 / (this%a + uext(j,k))
        this%ax(j,k)   = this%ax(j,k)   + (t*this%hx**2)
        this%ax(j+1,k) = this%ax(j+1,k) + (t*this%hx**2)
        this%ay(j,k)   = this%ay(j,k)   + (t*this%hy**2)
        this%ay(j,k+1) = this%ay(j,k+1) + (t*this%hy**2)
      end do
    end do
    this%ax = 2.0_r8 / this%ax
    this%ay = 2.0_r8 / this%ay
    do k = 1, this%ny
      do j = 1, this%nx
        this%ac(j,k) = this%ax(j,k) + this%ax(j+1,k) + this%ay(j,k) + this%ay(j,k+1)
      end do
    end do
  end subroutine update_system

  subroutine system_pc_ssor (this, nsweep, omega, r)
    class(system), intent(in)    :: this
    integer,       intent(in)    :: nsweep
    real(r8),      intent(in)    :: omega
    real(r8),      intent(inout) :: r(:,:)
    integer :: i, j, k
    real(r8) :: z(0:this%nx+1,0:this%ny+1)
    ASSERT(size(r,1) == this%nx)
    ASSERT(size(r,2) == this%ny)
    ASSERT(nsweep >= 1)
    ASSERT(omega > 0.0_r8)
    z = 0.0_r8
    do i = 1, nsweep
      !! Forward sweep
      do k = 1, this%ny
        do j = 1, this%nx
          z(j,k) = (1 - omega)*z(j,k) + omega*(r(j,k) &
                     + this%ax(j,k)*z(j-1,k) + this%ax(j+1,k)*z(j+1,k)  &
                     + this%ay(j,k)*z(j,k-1) + this%ay(j,k+1)*z(j,k+1))/this%ac(j,k)
        end do
      end do
      !! Backward sweep
      do k = this%ny, 1, -1
        do j = this%nx, 1, -1
          z(j,k) = (1 - omega)*z(j,k) + omega*(r(j,k) &
                     + this%ax(j,k)*z(j-1,k) + this%ax(j+1,k)*z(j+1,k)  &
                     + this%ay(j,k)*z(j,k-1) + this%ay(j,k+1)*z(j,k+1))/this%ac(j,k)
        end do
      end do
    end do
    !! Copy result to return array.
    r = z(1:this%nx,1:this%ny)
  end subroutine system_pc_ssor

end module system_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SOLVER_TYPE -- Defines the iterative solver for the discrete nonlinear system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use system_type
  use nka_type
  implicit none
  private

  type, public :: solver
    private
    type(system), pointer :: sys => null()
    integer :: nsweep
    real(r8) :: omega
    type(nka), allocatable :: accel
  contains
    procedure :: init => solver_init
    procedure :: solve => solver_solve
  end type solver

contains

  subroutine solver_init (this, sys, nsweep, omega, mvec)
    class(solver), intent(out) :: this
    type(system), intent(in), target :: sys
    integer,  intent(in) :: nsweep  ! number of SSOR preconditioner sweeps
    real(r8), intent(in) :: omega   ! SSOR preconditioner over-relaxation factor
    integer,  intent(in) :: mvec    ! NKA subspace size
    ASSERT(nsweep > 0)
    ASSERT(omega > 0.0_r8)
    ASSERT(mvec >= 0)
    this%sys => sys
    this%nsweep = nsweep
    this%omega  = omega
    if (mvec > 0) then
      allocate(this%accel)
      call this%accel%init (sys%size(), mvec=mvec)
    end if
  end subroutine solver_init

  subroutine solver_solve (this, uext)
    class(solver), intent(inout) :: this
    real(r8), intent(inout), target :: uext(0:,0:)  ! solution and boundary values
    integer :: itr, nx, ny
    real(r8) :: rnorm, rnorm0, red, rate
    real(r8), pointer :: u(:,:), r(:,:)
    real(r8), target :: r1(this%sys%size())
    integer, parameter :: MAXITR = 999
    real(r8), parameter :: TOL = 1.0d-6
    ASSERT(ubound(uext,1) == this%sys%xsize()+1)
    ASSERT(ubound(uext,2) == this%sys%ysize()+1)
    write(*,fmt='(a4,a14,a13,a8)') 'Iter', 'Residual Norm', 'Reduction', 'Rate'
    nx = this%sys%xsize()
    ny = this%sys%ysize()
    r(1:nx,1:ny) => r1    ! r1 is a rank-1 view of the residual array
    u => uext(1:nx,1:ny)  ! trim off the boundary data from the solution array
    call this%sys%residual (uext, r)
    rnorm0 = l2norm(r1)
    write(*,'(i3,a,es14.6)') 0, ':', rnorm0
    do itr = 1, MAXITR
      call this%sys%pc_ssor (this%nsweep, this%omega, r)
      if (allocated(this%accel)) call this%accel%accel_update (r1)
      u = u - r
      call this%sys%residual (uext, r)
      rnorm = l2norm(r1)
      red = rnorm / rnorm0
      rate = red**(1.0_r8/itr)
      write(*,fmt='(i3,a,es14.6,es13.3,f8.3)') itr, ':', rnorm, red, rate
      if (rnorm < TOL * rnorm0) exit
    end do
  end subroutine solver_solve

  real(r8) function l2norm (x)
    real(r8), intent(in) :: x(:)
    l2norm = sqrt(sum(x**2))
  end function l2norm

end module solver_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COMMAND_LINE -- Defines a parser for the command line arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module command_line

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: parse_command_line

  !! Command line parameters with default values.
  integer, public :: nx=50, ny=50, nsweep=2, mvec=0
  real(r8), public :: a=0.02_r8, omega=1.4_r8

  character(:), allocatable :: prog

contains

  subroutine parse_command_line

    character(64) :: arg
    integer :: n, num_arg, ios

    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.)
    prog = trim(arg(n+1:))

    n = 0
    num_arg = command_argument_count()

    do while (n < num_arg)
      n = n + 1
      call get_command_argument (n, arg)
      select case (arg)
      case ('--help')
        call usage_summary
      case ('-a')
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) a
        if (ios /= 0) call usage_halt ('invalid value for -a: ' // trim(arg))
        if (a <= 0.0_r8) call usage_halt ('invalid value for -a: must be >0')
      case ('-n')
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) nx
        if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))
        if (nx < 3) call usage_halt ('invalid value for -n: must be >2')
        ny = nx
      case ('--sweeps')
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) nsweep
        if (ios /= 0) call usage_halt ('invalid value for --sweeps: ' // trim(arg))
        if (nsweep <= 0) call usage_halt ('invalid value for --sweeps: must be >0')
      case ('--omega')
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) omega
        if (ios /= 0) call usage_halt ('invalid value for --omega: ' // trim(arg))
        if (omega <= 0.0_r8) call usage_halt ('invalid value for --omega: must be >0')
      case ('--nka-vec')
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(arg,*,iostat=ios) mvec
        if (ios /= 0) call usage_halt ('invalid value for --nka-vec: ' // trim(arg))
        if (mvec < 0) call usage_halt ('invalid value for --nka-vec: must be >=0')
      case default
        call usage_halt ('invalid option: ' // trim(arg))
      end select
    end do

  end subroutine parse_command_line

  subroutine usage_summary
    write(*,fmt='(a,/)') 'Usage: ' // prog // ' [options]'
    write(*,fmt='(a)')   'Solve the nonlinear elliptic equation -Div[(U+A)Grad[U]] = 1 on [0,1]^2'
    write(*,fmt='(a,/)') 'with zero boundary values, using a preconditioned fixed point iteration.'
    write(*,fmt='(a)') 'Options:'
    write(*,fmt='(a)') ' -a A         The equation parameter A (>0, default 0.02).'
    write(*,fmt='(a)') ' -n N         Number of grid cells in each coordinate direction (default 50).'
    write(*,fmt='(a)') ' --sweeps S   Use S SSOR sweeps as the preconditioner (>0, default 2).'
    write(*,fmt='(a)') ' --omega W    SSOR over-relaxation parameter W (>0, default 1.4).'
    write(*,fmt='(a)') ' --nka-vec M  Use an NKA acceleration subspace of dimension at most M;'
    write(*,fmt='(a)') '              0 turns off NKA acceleration (default).'
    stop
  end subroutine usage_summary

  subroutine usage_halt (message)
    character(*), intent(in) :: message
    write(*,fmt='(a)') prog // ': ' // message
    write(*,fmt='(a)') 'Use `' // prog // ' --help'' to get usage information.'
    stop
  end subroutine usage_halt

end module command_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program nka_example

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use system_type
  use solver_type
  use command_line
  implicit none

  call run

contains

  !! Wrapping things in a call forces automatic finalization of all objects,
  !! making it easier to find memory leaks.

  subroutine run

    type(solver) :: sol
    type(system), target :: sys
    real(r8), allocatable :: u(:,:)

    call parse_command_line

    !! Define the discrete problem.
    call sys%init (a=a, nx=nx, ny=ny)

    !! Define the solver and then solve.
    call sol%init (sys, nsweep=nsweep, omega=omega, mvec=mvec)
    allocate(u(0:sys%xsize()+1,0:sys%ysize()+1))
    u = 0.0_r8  ! initial solution guess plus boundary data
    call sol%solve (u)

  end subroutine run

end program nka_example

!!  F90_ASSERT -- C-style assertions for Fortran.
!!
!!    Neil N. Carlson <nnc@newmexico.com>
!!
!!  Usage: At the top of the source file, include the preprocessor file
!!  f90_assert.fpp which defines the Assert() preprocessor macro:
!!
!!    #include "f90_assert.fpp"
!!
!!  If the macro NDEBUG is defined (-D NDEBUG) when the file is passed
!!  through the preprocessor, lines of the form
!!
!!    Assert( <scalar logical expression> )
!!
!!  will be expanded to Fortran code which tests whether the logical
!!  expression is true, and if not, calls the following routine which
!!  will print the file name and line number and then halt execution.
!!  If the macro NDEBUG is not defined, then the Assert() is expanded
!!  to a Fortran comment line.
!!
!!  This is intentionally NOT a module procedure.
!!
!!  NB: Use with Fortran-aware preprocessors like fpp is robust.  One
!!  can use the C preprocessor cpp, but if the expanded macro extends
!!  the line past 132 characters, a compiler error will probably result.
!!
!!  NB: This version prints to unit 0 which is typically unix stderr.
!!

subroutine f90_assert (file, line)

  character(len=*), intent(in) :: file
  integer,          intent(in) :: line

  write(unit=0,fmt='(a,i4.4)') "Assertion failed at " // file // ":", line
  stop

end subroutine f90_assert

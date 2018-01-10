module module_compiler_dependant
!-----------------------------------------------------------------------
! Module to define entities that depend on gfortran compiler
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 12/01/2012
!
! PUBLIC CONSTANTS (they will be compiler-independant in future versions of Fortran):
!   real64, double: kind selector for real types requiring double precision.
!   input_unit:  standard input unit number
!   output_unit: standard output unit number
!   error_unit:  standard error unit number
!   iostat_end:  value returned by iostat= to indicate a end-of-file condition
!   iostat_eor:  value returned by iostat= to indicate a end-of-record condition
!
! PUBLIC PROCEDURES:
!   dir_exists: checks the existence of a directory
!   execute_command_line: call command shell
!-----------------------------------------------------------------------
implicit none

!Constants
integer, parameter, public :: real64 = selected_real_kind(15, 307) !already in iso_Fortran_env (gfortran 4.5.1)
integer, parameter, public :: double = selected_real_kind(15, 307)
integer, parameter, public :: input_unit  =  5 !already defined in module iso_Fortran_env (gfortran 4.3.3)
integer, parameter, public :: output_unit =  6 !already defined in module iso_Fortran_env (gfortran 4.3.3)
integer, parameter, public :: error_unit  =  0 !already defined in module iso_Fortran_env (gfortran 4.3.3)
integer, parameter, public :: iostat_end  = -1 !already defined in module iso_Fortran_env (gfortran 4.3.3)
integer, parameter, public :: iostat_eor  = -2 !already defined in module iso_Fortran_env (gfortran 4.3.3)

contains

!-----------------------------------------------------------------------
! dir_exists: checks the existence of a directory
!-----------------------------------------------------------------------
function dir_exists(dirname) result(res)
character(*) :: dirname
logical :: res
integer :: ios

inquire(file = dirname, exist = res, iostat = ios)
if (ios /= 0) res = .false.
end function

!-----------------------------------------------------------------------
! execute_command_line: call command shell
! REMARK: its header agrees with the Fortran 2008 intrinsic procedure
!-----------------------------------------------------------------------
subroutine execute_command_line(command, wait, exitstat, cmdstat, cmdmsg)
character(*),           intent(in)    :: command
logical,      optional, intent(in)    :: wait
integer,      optional, intent(inout) :: exitstat
integer,      optional, intent(out)   :: cmdstat
character(*), optional, intent(inout) :: cmdmsg
integer :: res, ierrno

call system(command, res)
if (present(cmdstat)) cmdstat = res
if (present(exitstat)) exitstat = ierrno() 
end subroutine

end module

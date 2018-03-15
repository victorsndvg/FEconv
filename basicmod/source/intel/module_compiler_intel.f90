module module_compiler_dependant_bmod
  !! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to define entities that depend on Intel compiler
!!
!! @note Several constants are declared in intrinsic module `iso_Fortran_env`:
!!   `input_unit`, `output_unit`, `error_unit`, `iostat_end`, `iostat_eor` (supported in ifort since version 11.1) and
!!   `real64` (supported in ifort since version 12.0).
!! For compatibility issues, to call the aforementioned constants from `module_compiler_dependant` is still permitted.
!! __However, we recommend to use such constants directly from `iso_Fortran_env` whenever it is possible.__
!!
!! @note `execute_command_line` is an intrinsic Fortran 2008 procedure (supported in ifort since version 15.0).
!! For compatibility issues, the aforementioned procedure can be called from `basicmod`.
!! __However, we recommend to use the intrinsic procedure `execute_command_line` whenever it is possible.__
!!
!! @warning If you use a version of ifort older than 12.0, you must comment here the use of module `iso_Fortran_env` and
!! uncomment the manual declaration of the aforementioned constants.
!!
!! @warning If you use the GFortran compiler, you must compile your code including `module_compiler_gfortran.f90` instead.
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
!   modification time: get modification time of a file
!-----------------------------------------------------------------------
use ifport, only: system, ierrno, GETFILEINFOQQ, FILE$INFO, FILE$FIRST, FILE$ERROR, GETLASTERRORQQ, ERR$NOMEM, ERR$NOENT, &
                  UNPACKTIMEQQ
use iso_Fortran_env, only: real64, input_unit, output_unit, error_unit, iostat_end, iostat_eor
implicit none
!Constants
integer, parameter, public :: double = selected_real_kind(15, 307) !! Alternative definition of `real64`.

!Manual declaration of constants (uncomment for versions of ifort older than 12.0)
!integer, parameter, public :: real64 = selected_real_kind(15, 307) ! Kind selector for real types requiring double precision.
!integer, parameter, public :: double = selected_real_kind(15, 307) ! Alternative definition of real64.
!integer, parameter, public :: input_unit  =  5  ! Standard input unit number.
!integer, parameter, public :: output_unit =  6 ! Standard output unit number.
!integer, parameter, public :: error_unit  =  0 ! Standard error unit number.
!integer, parameter, public :: iostat_end  = -1 ! Value returned by iostat= to indicate a end-of-file condition.
!integer, parameter, public :: iostat_eor  = -2 ! Value returned by iostat= to indicate a end-of-record condition.

contains

!-----------------------------------------------------------------------
! dir_exists
!-----------------------------------------------------------------------
function dir_exists(dirname, ios) result(res)
!! Checks the existence of a directory.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dir_exists`  
!! `implicit none`  
!! `integer :: ios`  
!! `if (.not. dir_exists('./code/err.txt')) print*, ios`  
!! `end program`  
!!
!! @note This procedure is compiler-dependant because GFortran and ifort use different arguments in `inquire` to check
!! the existence of a directory.
character(*), intent(in)  :: dirname !! Directory name.
integer,      intent(out) :: ios     !! I/O error number.
logical                   :: res     !! `.true.` if directory exists; `.false.` otherwise.

inquire(directory = dirname, exist = res, iostat = ios)
end function

!-----------------------------------------------------------------------
! execute_command_line
!-----------------------------------------------------------------------
subroutine execute_command_line(command, wait, exitstat, cmdstat, cmdmsg)
!! Calls command shell.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: execute_command_line`  
!! `implicit none`  
!! `integer :: exits, cmds`  
!! `call execute_command_line('program.exe', exitstat=exits, cmdstat=cmds)`  
!! `if (cmds /= 0) write(error_unit,'(a)') cmds`  
!! `print*, exits`  
!! `end program`  
!!
!! @note Since `execute_command_line` is a Fortran 2008 intrinsic procedure, we only recommend to use the subroutine defined here
!! with old versions of the compiler that do not support Fortran 2008.
character(*),           intent(in)    :: command  !! Holds the command line to be executed.
logical,      optional, intent(in)    :: wait     !! Disabled, for compatibility purposes only.
integer,      optional, intent(inout) :: exitstat !! Returns the value of ierrno.
integer,      optional, intent(out)   :: cmdstat  !! Returns a non-zero if some problem appears; zero otherwise.
character(*), optional, intent(inout) :: cmdmsg   !! Disabled, for compatibility purposes only.
integer :: res

res = system(command)
if (present(cmdstat)) cmdstat = res
if (present(exitstat)) exitstat = ierrno()
end subroutine

!-----------------------------------------------------------------------
! modification time
!-----------------------------------------------------------------------
function modification_time(filename, errmsg) result(values)
!! Get the modification time of a file.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: modification_time, maxpath, error`  
!! `implicit none`  
!! `character(maxpath) :: errmsg`  
!! `integer :: values(8)`  
!! `values = modification_time('file.txt', errmsg)`  
!! `if (values(1) <= 0) then`  
!! &nbsp;&nbsp;`call error(trim(errmsg))`  
!! `else`  
!! &nbsp;&nbsp;`print*, values`  
!! `end if`  
!! `end program`    
!!
!! @warning `GETFILEINFOQQ` requires, as a second argument, an `INTEGER(4)` on IA-32 architecture and an `INTEGER(8)`
!! on Intel(R) 64 and IA-64 architectures. If you pretend to compile on IA-32 architecture, you must manually change the integer
!! kind in the declaration of `handle`.
character(*),           intent(in)  :: filename   !! Filename.
character(*), optional, intent(out) :: errmsg     !! Error message.
integer                             :: values(8)  !! Return `values(1) = -1` if error; otherwise, `values` with modification time
                                                  !! in the `date_and_time()` format.
integer(4) :: length
integer(8) :: handle !On Intel(R) 64 architecture and IA-64 architecture
!integer(4) :: handle !On IA-32 architecture
integer(2) :: iyr, imon, iday, ihr, imin, isec
type(FILE$INFO) :: buffer

handle = FILE$FIRST
length = GETFILEINFOQQ(filename, buffer, handle)

if (handle == FILE$ERROR) then
  select case(GETLASTERRORQQ())
  case(ERR$NOMEM)
    if (present(errmsg)) errmsg = '(module_compiler_intel/modification time) Out of memory'
  case(ERR$NOENT)
    if (present(errmsg)) errmsg = '(module_compiler_intel/modification time) The file or path specified was not found: '//&
    &trim(filename)
  case default
    if (present(errmsg)) errmsg = '(module_compiler_intel/modification time) Unable to find file or path name: '//trim(filename)
  end select
  values(1) = -1
elseif (length == 0) then
  if (present(errmsg)) errmsg = '(module_compiler_intel/modification time) No matching files were found: '//trim(filename)
  values(1) = 0
else
  call UNPACKTIMEQQ(buffer%LASTWRITE, iyr, imon, iday, ihr, imin, isec)
  !convert to the data_and_time() format
  values(1) = iyr  !4-digit year
  values(2) = imon !month of the year
  values(3) = iday !day of the month
  values(4) = -1   !time difference with respect to UTC in minutes (NOT SET)
  values(5) = ihr  !hour of the day (0 to 23) - local time
  values(6) = imin !minutes of the hour (0 to 59) - local time
  values(7) = isec !seconds of the minute (0 to 59) - local time
  values(8) = -1   !milliseconds of the second (0 to 999) - local (NOT SET)
end if
end function
end module

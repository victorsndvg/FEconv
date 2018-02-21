module module_report_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to manage messages and errors.  
! 
! PRIVATE CONSTANTS:
!   DEBUG_: index for level and output arrays
!   INFO_:  index for level and output arrays
!   WARN_:  index for level and output arrays
!   ERROR_: index for level and output arrays
! 
! PUBLIC PROCEDURES:
!   report_option: sets the report options
!   report: reports a message for a given level
!   debug: reports a debug message
!   info: reports an info
!   warning: reports a warning
!   error: reports an error
! 
! PRIVATE PROCEDURES:
!   get_unit: seeks a non-connected unit number
!   add_sep: adds a slash (if necessary) at the end of the string
!   lcase: converts a string to lower case
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: error_unit, output_unit
use module_os_dependant_bmod, only: maxpath, slash

implicit none

!Constants
integer, private :: DEBUG_ = 1 !! DEBUG   index for level and output arrays.
integer, private :: INFO_  = 2 !! INFO    index for level and output arrays.
integer, private :: WARN_  = 3 !! WARNING index for level and output arrays.
integer, private :: ERROR_ = 4 !! ERROR   index for level and output arrays.

!Variables
character(7), dimension(4), private :: level   = ['DEBUG  ', 'INFO   ', 'WARNING', 'ERROR  '] !! Existing report levels.
character(4), dimension(4), private :: output  = ['none',    'none',    'std ',    'std ']    !! Actual report output.
character(maxpath), private :: rfile  = 'OUTPUT.txt' !! Actual report file.
character(maxpath), private :: folder = '.'          !! Actual report folder.

!Private procedures
private :: get_unit, add_sep, lcase

contains

!-----------------------------------------------------------------------
! report_option: sets the report options
!-----------------------------------------------------------------------
subroutine report_option(opt, sval)
!! Sets the report options.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: report_option`  
!! `implicit none`  
!! `call report_option('folder', './out')`  
!! `call report_option('file', './errors.txt')`  
!! `call report_option('output', 'std')`  
!! `call report_option('error', 'file')`  
!! `! Error messages are sent to a file, the rest of them are sent to the standard output`  
!! `end program`  
!!
!! @note Option `folder` allows you to change the folder where the outfput file is saved. By default, `.`  
!!
!! @note Option `file` allows you to change the output file. By default, `OUTPUT.txt`  
!!
!! @note Options `debug`, `info`, `warning` and `error` allows you to change the type of output for the given level, 
!! whereas option `output` allows you to change the type of output for all the levels. In all cases, the valid values are:  
!! * `none`, no output.  
!! * `file`, output is only redirected to a file.  
!! * `std`, output is only redirected to the standard output (`stderr` for errors, `stdout` otherwise).  
!! * `all`, output is redirected to a file and to the standard output.  
!! By default, output is `none` for `debug` and `info`, and `std` for `warning` and `error`.  
!!
!! @warning An option that changes the type of output for a level overwrites the previous settings for this level.  
character(*), intent(in) :: opt  !! Option; it can be `file`, `folder`, `output`, `debug`, `info`, `warning` or `error`.
character(*), intent(in) :: sval !! Value for the given option. See valid values in the notes of the [[report_option(subroutine)]]
                                 !! documentation.

select case(lcase(opt))
case('file')
  rfile = sval
case('folder')
  folder = add_sep(sval)
case('output', 'debug', 'info', 'warning', 'error')
  if (sval /= 'none' .and. sval /= 'std' .and. sval /= 'file' .and. sval /= 'all') &
    call error('(report/set_option) this option requires the 2nd argument be none|std|file|all.')
  select case(lcase(opt))
  case('output');  output         = sval
  case('debug');   output(DEBUG_) = sval
  case('info');    output(INFO_)  = sval
  case('warning'); output(WARN_)  = sval
  case('error');   output(ERROR_) = sval
  end select
case default
  call error('(report/set_option) "'//trim(opt)//'" is not one of the valid options: file, folder, output, debug, info, &
  &warning, error.')
end select
end subroutine

!-----------------------------------------------------------------------
! report: reports a message for a given level
!-----------------------------------------------------------------------
subroutine report(lev, msg)
!! Reports a message for a given level.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: report`  
!! `implicit none`  
!! `call report('debug',   'This is a debug message.')`  
!! `call report('info',    'This is an informative message.')`  
!! `call report('warning', 'This is a warning message.')`  
!! `call report('error',   'This is an error message; it also stops the program.')`  
!! `end program`  
!!
!! @note Valid levels are `"debug"`, `"info"`, `"warning"` and `"error"`.  
!!
!! @note The standard output means `stderr` for errors and `stdout` otherwise.  
!!
!! @warning A message for the level `"error"` is always followed by a stop with exit code 1.  
character(*), intent(in) :: lev !! Level.
character(*), intent(in) :: msg !! Message.
character(8)  :: date
character(10) :: time
character(20) :: dt
integer :: iu, il

select case(lcase(lev))
case('debug');   il = DEBUG_
case('info');    il = INFO_
case('warning'); il = WARN_
case('error');   il = ERROR_
case default
  write(error_unit,'(a)') 'ERROR: (report/report) "'//trim(lev)//'" is not one of the valid levels: debug, info, warning, error.'
  stop 1
end select
!file write
if (output(il) == 'file' .or. output(il) == 'all') then
  iu = get_unit()
  if (iu > 0) then
    open(iu, file = trim(add_sep(folder))//adjustl(rfile), position = "append")
    call date_and_time(date, time)
    dt = date(7:8)//'/'//date(5:6)//'/'//date(1:4)//'  '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(iu,*) dt//'  '//trim(level(il))//': '//trim(msg)
    close(iu)
  end if
end if
!standard write
if (output(il) == 'std' .or. output(il) == 'all') then
  if (il == ERROR_) then
    write(error_unit,'(a)')  trim(level(il))//': '//trim(msg)
    stop 1
  else
    write(output_unit,'(a)') trim(level(il))//': '//trim(msg)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! debug: reports a debug message
!-----------------------------------------------------------------------
subroutine debug(msg)
!! Reports a debug message.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: debug`  
!! `implicit none`  
!! `call debug('This is a debug message.')`  
!! `end program`  
!!
!! @note `call debug(msg)` is equivalent to `call report('debug', msg)`.  
character(*), intent(in) :: msg !! Message.

call report('DEBUG', msg)
end subroutine

!-----------------------------------------------------------------------
! info: reports an info
!-----------------------------------------------------------------------
subroutine info(msg)
!! Reports an informative message.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: info`  
!! `implicit none`  
!! `call info('This is an informative message.')`  
!! `end program`  
!!
!! @note `call info(msg)` is equivalent to `call report('info', msg)`.  
character(*), intent(in) :: msg !! Message.

call report('INFO', msg)
end subroutine

!-----------------------------------------------------------------------
! warning: reports a warning
!-----------------------------------------------------------------------
subroutine warning(msg)
!! Reports a warning message.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: warning`  
!! `implicit none`  
!! `call warning('This is a warning message.')`  
!! `end program`  
!!
!! @note `call warning(msg)` is equivalent to `call report('warning', msg)`.  
character(*), intent(in) :: msg !! Message.

call report('WARNING', msg)
end subroutine

!-----------------------------------------------------------------------
! error: report an error
!-----------------------------------------------------------------------
subroutine error(msg)
!! Reports an error and stops the program with exit code 1.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: error`  
!! `implicit none`  
!! `call error('This is an error message; it also stops the program.')`  
!! `end program`  
!!
!! @note `call error(msg)` is equivalent to `call report('error', msg)`.  
character(*), intent(in) :: msg !! Message.

call report('ERROR', msg)
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! add_sep: add a slash (if necessary) at the end of the string
!-----------------------------------------------------------------------
function add_sep(cad) result(res)
 !! Adds a slash (if necessary) at the end of the string.  
character(*) :: cad !! Input string
character(len(cad) + 1) :: res !! Output string

res = cad
if (index(cad,slash(),back=.true.).ne.len_trim(cad)) res = trim(cad)//slash()
end function

!-----------------------------------------------------------------------
! get_unit: get a non-connected unit number
! Remark: since module_files uses this module, we cannot call get_unit from there
!-----------------------------------------------------------------------
function get_unit() result(next)
 !! Gets a non-connected unit number.  
 !!
 !! @note Since `module_files` uses this module, we cannot call its `get_unit` from there  
implicit none
integer :: next !! Non-connected unit number
integer, parameter :: min = 10, max = 999
logical :: open_

do next = min, max
  inquire(unit = next, opened = open_)
     if (.not. open_) return
end do
next = -1
end function

!-----------------------------------------------------------------------
! lcase: converts a string to lower case
! Remark: since module_convers uses this module, we cannot call lcase from there
!-----------------------------------------------------------------------
function lcase(str) result(res)
 !! Converts a string to lower case.
 !!
 !! @note Since `module_convers` uses this module, we cannot call its `lcase` from there.
character(*), intent(in) :: str !! Input string
character(len(str))      :: res !! Output string in lower case
integer :: diff, i

res = str
diff = ichar('A') - ichar('a')
do i = 1, len_trim(str)
  if (str(i:i) < 'A' .or. str(i:i) > 'Z') cycle
  res(i:i) = char(ichar(str(i:i)) - diff)
end do
end function

end module

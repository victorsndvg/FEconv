module module_report
!-----------------------------------------------------------------------
! Module for messages and errors management
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 07/03/2012
!
! PUBLIC PROCEDURES:
!   report_option: sets the report options (stop, file, folder and level)
!   error: sends a report and stops the program (depending on break)
!   info: sends a report
!   error_found: checks whether or not an error was found
!
! PRIVATE PROCEDURES:
!   get_unit: seeks a non-connected unit number
!   add_sep: adds a slash (if necessary) at the end of the string
!-----------------------------------------------------------------------
use module_compiler_dependant, only: error_unit, output_unit
use module_os_dependant, only: maxpath, set_os, slash
implicit none

!Variables
logical,            private :: break  = .true.       !whether or not to stop the program in error()
character(maxpath), private :: rfile  = 'OUTPUT.txt' !actual report file
character(maxpath), private :: folder = '.'          !actual report folder
character(maxpath), private :: level  = 'stdout'     !actual report level ('none'|'stdout'|'file'|'all')
logical,            private :: errfound = .false.    !whether or not an error was found

!Private procedures
private :: get_unit, add_sep

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! report_option: sets the report options
!-----------------------------------------------------------------------
subroutine report_option(opt, s, l)
character(*), intent(in) :: opt
character(*), optional, intent(in) :: s
logical,      optional, intent(in) :: l

select case(opt)
case('stop')
  if (.not. present(l)) call error('(report/set_option) option "stop", value must be logical.')
  break = l
case('file')
  if (.not. present(s)) call error('(report/set_option) option "file", value must be character.')
  rfile = s
case('folder')
  if (.not. present(s)) call error('(report/set_option) option "folder", value must be character.')
  folder = add_sep(s)
case('level')
  if (.not. present(s)) call error('(report/set_option) option "level", value must be character.')
  level = s
  if (level == 'file' .or. level == 'all') call set_os()
case default
  call error('(report/set_option) "'//trim(opt)//'" is not a report option: stop, file, folder, level.')
end select

end subroutine
  
!***********************************************************************
! OUTPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! error: report an error
!-----------------------------------------------------------------------
subroutine error(err)
character(*), intent(in) :: err
character(8)  :: date
character(10) :: time
character(20) :: dt
integer :: iu

errfound = .true.
!file write
if (level == 'file' .or. level == 'all') then
  iu = get_unit()
  if (iu > 0) then
    open(iu, file = trim(add_sep(folder))//adjustl(rfile), position = "append")
    call date_and_time(date, time)
    dt = date(7:8)//'/'//date(5:6)//'/'//date(1:4)//'  '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(iu,*) dt//'  '//'ERROR: '//trim(err)
    close(iu)
  end if
end if
!standard write
if (level == 'stdout' .or. level == 'all') then
  write(error_unit,'(a)') 'ERROR: '//trim(err)
end if
!program stops
 if (break) stop 1
end subroutine

!-----------------------------------------------------------------------
! info: report an info
!-----------------------------------------------------------------------
subroutine info(inf)
character(*), intent(in) :: inf
character(8)  :: date
character(10) :: time
character(20) :: dt
integer :: iu

!file write
if (level == 'file' .or. level == 'all') then
  iu = get_unit()
  if (iu > 0) then
    open(iu, file = trim(add_sep(folder))//adjustl(rfile), position = "append")
    call date_and_time(date, time)
    dt = date(7:8)//'/'//date(5:6)//'/'//date(1:4)//'  '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(iu,*) dt//'  '//'INFO: '//trim(inf)
    close(iu)
  end if
end if
!standard write
if (level == 'stdout' .or. level == 'all') then
  write(output_unit,'(a)') 'INFO: '//trim(inf)
end if
end subroutine

!-----------------------------------------------------------------------
! error_found: checks whether or not an error was found
!-----------------------------------------------------------------------
function error_found() result(res)
logical :: res
res = errfound
end function

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! add_sep: add a slash (if necessary) at the end of the string
!-----------------------------------------------------------------------
function add_sep(cad) result(res)
character(*) :: cad
character(len(cad) + 1) :: res

res = cad
if (index(cad,slash(),back=.true.).ne.len_trim(cad)) res = trim(cad)//slash()
end function

!-----------------------------------------------------------------------
! get_unit: get a non-connected unit number
! Remark: since module_files uses this module, we cannot call get_unit from there
!-----------------------------------------------------------------------
function get_unit() result(next)
implicit none
integer :: next
integer, parameter :: min = 10, max = 999
logical :: open_
do next = min, max
  inquire(unit = next, opened = open_)
     if (.not. open_) return
end do
next = -1
end function

end module

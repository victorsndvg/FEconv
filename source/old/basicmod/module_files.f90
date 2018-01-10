module module_files
!-----------------------------------------------------------------------
! Module for files management
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 20/01/2012
!
! PUBLIC PROCEDURES:
!   get_unit: seeks a non connected unit number
!   file_exists: checks the existence of a file
!   get_name: removes the path from a filename
!-----------------------------------------------------------------------
use module_os_dependant, only: slash
use module_report, only: error, info
use module_convers, only: string
implicit none

contains

!--------------------------------------------------------------------
! get_unit: gets a non-connected unit number
!--------------------------------------------------------------------
function get_unit() result(next)
integer :: next, ios
integer, parameter :: min = 10, max = 999
logical :: open_

do next = min, max
  inquire(unit = next, opened = open_, iostat = ios)
     if (ios /= 0) call error('get_unit/inquire, #'//trim(string(ios)))
     if (.not. open_) return
end do
call error('get_unit, unused unit number not found')
end function

!--------------------------------------------------------------------
! file_exists: checks the existence of a file
!--------------------------------------------------------------------
function file_exists(filename) result(res)
character(*) :: filename
integer :: iu, ios
logical :: res, exists, already_open

res =.true.
inquire(file = filename, exist = exists, opened = already_open, iostat = ios)
if (ios /= 0) call error('file_exists/inquire, #'//trim(string(ios)))
if (.not. exists) then
  res = .false.; return
elseif (.not. already_open) then
  iu = get_unit()
  open(iu, file = filename, iostat = ios)
  if (ios /= 0) then
    res = .false.; return
  else
     close(iu)
  end if
end if
end function

!-----------------------------------------------------------------------
! get_name: removes the path from a filename
! REMARK: depends on set_os() via slash()
!-----------------------------------------------------------------------
function get_name(str) result(res)
character(*), intent(in) :: str
character(len(str)) :: res
integer :: p

p = index(str, slash(), back = .true.)
res = str(p+1:len_trim(str))
end function

end module

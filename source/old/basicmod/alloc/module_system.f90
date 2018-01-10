module module_system
!-----------------------------------------------------------------------
! Module for system commands
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 20/01/2012
!
! PUBLIC PROCEDURES:
!   ls: list the files of a folder
!   mkdir: create a directory
!   rm: delete a file
!-----------------------------------------------------------------------
use module_compiler_dependant, only: execute_command_line, dir_exists
use module_os_dependant, only: get_os, slash, maxpath
use module_report, only: info
use module_convers, only: string
use module_alloc_char_r1, only: dealloc, set, reduce
use module_files, only: get_unit
implicit none

contains

!-----------------------------------------------------------------------
! ls: list the files of 'folder'
!-----------------------------------------------------------------------
function ls(folder, list) result(ios)
character(*), intent(in)  :: folder
character(*), allocatable :: list(:)
character(maxpath) :: fichero, fname
integer :: ios, pos, iu, iu2, es
logical :: ex, op

call dealloc(list)
!create list of files
if (trim(get_os()) == 'windows') then
  call execute_command_line('dir /b '//trim(folder)//'> ls_files_auxiliary_file.txt', cmdstat=ios, exitstat = es)
else
  call execute_command_line('ls -1 '//trim(folder)//'> ls_files_auxiliary_file.txt', cmdstat=ios, exitstat = es)
endif
if (ios /= 0) then; call info('ls_files/system, error #'//trim(string(es))); return; endif
!read file
iu = get_unit()
open (unit=iu, file='ls_files_auxiliary_file.txt', form='formatted', iostat=ios, status='old', position='append')
if (ios /= 0) then; call info('ls_files/open, #'//trim(string(ios))); return; endif
!avoid ios = -1 in the last line
write(iu, *) ' '
rewind(iu)
pos = 0
do
  read (unit=iu, fmt='(a)', iostat=ios) fichero
  if (ios /= 0) exit
  !check filename
  fname = trim(folder)//slash()//trim(fichero)
  inquire(file=fname, exist=ex, opened=op, iostat=ios)
  if (ios /= 0 .or. .not. ex) cycle
  if (.not. op) then
    iu2 = get_unit()
    open (unit=iu2, file=fname, iostat=ios, status='old')
    if (ios /= 0) then; close(unit=iu2, iostat=ios); cycle; end if
    close(unit=iu2)
  end if
  !add it to the list
  call set(list, fichero, pos+1, fit=.false.); pos = pos+1
enddo
close(iu)
ios = rm('ls_files_auxiliary_file.txt')
call reduce(list, pos)
end function

!-----------------------------------------------------------------------
! mkdir: create a directory
!-----------------------------------------------------------------------
function mkdir(folder) result(res)
character(*), intent(in) :: folder
integer :: res

call execute_command_line('mkdir '//trim(folder), cmdstat=res)
if (.not. dir_exists(folder)) res = -1
end function

!-----------------------------------------------------------------------
! rm: delete a file
!-----------------------------------------------------------------------
function rm(str) result(res)
character(*), intent(in) :: str
integer :: res

if (trim(get_os()) == 'windows') then
  call execute_command_line('del '//trim(str), cmdstat=res)
else
  call execute_command_line('rm  '//trim(str), cmdstat=res)
endif
if (res/=0) call info('del_file, file cannot be deleted: '//trim(str))
end function

end module

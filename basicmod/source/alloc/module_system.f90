module module_system_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for system call commands.  
!!
!! @warning Procedures [[ls(function)]] and [[rm(function)]] use the OS-dependant function 
!! `[[module_os_dependant_bmod(module):get_os(function)]]` to select the terminal command. Please see documentation of 
!! procedure `[[module_os_dependant_bmod(module):set_os(subroutine)]]` to be aware of the OS-detection method.  
!!
!
! PUBLIC PROCEDURES:
!   ls: list the files of a folder
!   mkdir: create a directory
!   rm: delete a file
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: execute_command_line, dir_exists
use module_os_dependant_bmod, only: get_os, slash, maxpath
use module_report_bmod, only: info
use module_convers_bmod, only: string
use module_alloc_char_r1_bmod, only: dealloc, set, reduce
use module_files_bmod, only: get_unit
implicit none

contains

!-----------------------------------------------------------------------
! ls
!-----------------------------------------------------------------------
function ls(folder, list) result(ios)
!! List the files included in a folder.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: maxpath, trim, error, ls`  
!! `implicit none`  
!! `character(maxpath), allocatable :: list(:)`  
!! `integer :: i`  
!! `if (ls('.', list) /= 0) call error('Unsuccessful ls')`  
!! `print*, (trim(list(i)),' ', i=1,size(list, 1))`  
!! `end program`  
!!
!! @note This procedure creates the temporary file *ls_files_auxiliary_file.txt*; the file is deleted at the end with 
!! `[[rm(function)]]`.  
!!
!! @note This procedure uses the compiler-dependant subroutine 
!! `[[module_compiler_dependant_bmod(module):execute_command_line(subroutine)]]` to make a system call.  
!!
!! @warning This procedure uses the OS-dependant function `[[module_os_dependant_bmod(module):get_os(function)]]` 
!! to select the terminal command. Please see documentation of procedure 
!! `[[module_os_dependant_bmod(module):set_os(subroutine)]]` to be aware of the OS-detection method.  
!!
!! @warning `ls` only lists files, not subfolders.  
character(*), intent(in)  :: folder   !! Folder.
character(*), allocatable :: list(:)  !! Returned file list.
integer :: ios                        !! Zero if successful; non-zero otherwise.
character(maxpath) :: fichero, fname
integer :: pos, iu, iu2, es
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
! mkdir
!-----------------------------------------------------------------------
function mkdir(folder) result(res)
!! Create a folder.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: error, mkdir`  
!! `implicit none`  
!! `if (mkdir('test_folder') /= 0) call error('Unsuccessful mkdir')`  
!! `end program`  
!!
!! @note This procedure uses the compiler-dependant procedures 
!! `[[module_compiler_dependant_bmod(module):execute_command_line(subroutine)]]` to make a system call and 
!! `[[module_compiler_dependant_bmod(module):dir_exists(function)]]` to check whether or not the folder was created.  
character(*), intent(in) :: folder !! Folder.
integer                  :: res    !! Zero if successful; non-zero otherwise.
integer :: ios

call execute_command_line('mkdir '//trim(folder), cmdstat=res)
if (.not. dir_exists(folder, ios)) res = -1
end function

!-----------------------------------------------------------------------
! rm
!-----------------------------------------------------------------------
function rm(str) result(res)
!! Delete a file.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: error, rm`  
!! `implicit none`  
!! `if (rm('test_folder') /= 0) call error('Unsuccessful rm')`  
!! `end program`  
!!
!! @note This procedure uses the compiler-dependant subroutine 
!! `[[module_compiler_dependant_bmod(module):execute_command_line(subroutine)]]` to make a system call.  
!!
!! @warning This procedure uses the OS-dependant function `[[module_os_dependant_bmod(module):get_os(function)]]` 
!! to select the terminal command. Please see documentation of procedure 
!! `[[module_os_dependant_bmod(module):set_os(subroutine)]]` to be aware of the OS-detection method.  
character(*), intent(in) :: str !! File or folder specification.
integer                  :: res !! Zero if successful; non-zero otherwise.

if (trim(get_os()) == 'windows') then
  call execute_command_line('del '//trim(str), cmdstat=res)
else
  call execute_command_line('rm  '//trim(str), cmdstat=res)
endif
if (res/=0) call info('(module_system/rm) file cannot be deleted: '//trim(str))
end function

end module

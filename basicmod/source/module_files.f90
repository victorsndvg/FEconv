module module_files_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! This module eases the work with files.
!!
!! @note `get_unit` performs the same work than the Fortran 2008 intrinsic function `newunit`; `get_unit` is included 
!! here for compatibility issues, but __we recommend to use the intrinsic procedure `newunit` whenever it is possible.__  
! 
! PUBLIC PROCEDURES:
!   get_unit: seeks a non connected unit number
!   file_exists: checks the existence of a file
!   get_name: removes the path from a filename
!   .IsNewerThan.: checks the modification times of two files.
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: modification_time
use module_os_dependant_bmod, only: slash
use module_report_bmod, only: error, info
use module_convers_bmod, only: string
implicit none

!Interfaces
interface operator (.IsNewerThan.);  module procedure is_newer_than_prv; end interface

private :: is_newer_than_prv

contains

!--------------------------------------------------------------------
! get_unit
!--------------------------------------------------------------------
function get_unit() result(next)
!! Gets a non-connected unit number; it searches unit numbers from 10 to 999 using `inquire` and return an `inquire` error, 
!! if exists.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: get_unit`  
!! `implicit none`  
!! `integer :: u`  
!! `u=get_unit()`  
!! `open(u,'file.txt')`  
!! `end program`  
!!
!! @note `get_unit` performs the same work than the Fortran 2008 intrinsic function `newunit`; `get_unit` is included 
!! here for compatibility issues, but __we recommend to use the intrinsic procedure `newunit` whenever it is possible.__
integer :: next !! Non-connected unit number.
integer :: ios
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
!! Checks the existence of a file.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: file_exists`  
!! `implicit none`  
!! `print*, file_exists('file.txt')`  
!! `end program`  
!!
!! @note `file_exists` checks the existence of the file using `inquire`; it will return any error produced by `inquire`. 
character(*) :: filename !! Filename to verify their existence.
logical      :: res      !! `.true.` if the file exists; `.false.` otherwise.
integer :: iu, ios
logical :: exists, already_open

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
! get_name
!-----------------------------------------------------------------------
function get_name(str) result(res)
!! Removes the path from a filename.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: get_name`  
!! `implicit none`  
!! `print*, get_name('../file.txt')`  
!! `end program`  
!!
!! @note `get_name` depends on `slash` to know the folder separation symbol; please check Notes 
!! of procedure `[[module_os_dependant_bmod(module):set_os(subroutine)]]` to be aware of the OS-detection method.  
character(*), intent(in) :: str !! Filename.
character(len(str))      :: res !! Name of the file after removing the path.
integer :: p

p = index(str, slash(), back = .true.)
res = str(p+1:len_trim(str))
end function

!--------------------------------------------------------------------
! .IsNewerThan.
!--------------------------------------------------------------------
function is_newer_than_prv(f1, f2) result(res)
!! Checks the modification times of two files.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: .isNewerThan.`  
!! `implicit none`  
!! `print*, 'file1.txt' .IsNewerThan. 'file2.txt'`  
!! `end program`  
character(*), intent(in) :: f1 !! First filename.
character(*), intent(in) :: f2 !! Second filename.
logical :: res
integer :: v1(8), v2(8)

if (.not. file_exists(f1) .or. .not. file_exists(f2)) then
  res = .false.
  return
end if
v1 = modification_time(f1)
v2 = modification_time(f2)

res = .true.
if (v1(1) >  v2(1)) return
if (v1(1) == v2(1) .and. v1(2) >  v2(2)) return
if (v1(1) == v2(1) .and. v1(2) == v2(2) .and. v1(3) >  v2(3)) return
if (v1(1) == v2(1) .and. v1(2) == v2(2) .and. v1(3) == v2(3) .and. v1(5) >  v2(5)) return
if (v1(1) == v2(1) .and. v1(2) == v2(2) .and. v1(3) == v2(3) .and. v1(5) == v2(5) .and. v1(6) >  v2(6)) return
if (v1(1) == v2(1) .and. v1(2) == v2(2) .and. v1(3) == v2(3) .and. v1(5) == v2(5) .and. v1(6) == v2(6) .and. v1(7) >  v2(7)) return
res = .false.
end function

end module

module module_set_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for set management.  
!!
!! @note In GFortran 4.8 a generic interface cannot contain both functions and subroutines. 
!! This is why pairs like `[[module_set_bmod(module):unique(function)]]` and `[[module_set_bmod(module):sunique(subroutine)]]`
!! have different public names.  
!!
!! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
!! Thus, since `[[module_set_bmod(module):unique(function)]]` returns an array, this code is correct:  
!! `integer, allocatable :: x(:)`  
!! `allocate(x, source=unique(...))`  
!! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
!! `integer, allocatable :: x(:)`  
!! `x = unique(...)`  
!! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
!! `[[module_set_bmod(module):sunique(subroutine)]]`.
!
! PUBLIC PROCEDURES:
!   unique: returns the values of the vector with no repetitions
!   sunique: subrutine for unique
!   setdiff: returns the values in one vector that are not in the other
!   ssetdiff: subrutine for setdiff
!-----------------------------------------------------------------------
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error
implicit none

contains

!-----------------------------------------------------------------------
! sunique
!-----------------------------------------------------------------------
subroutine sunique(v, res)
!! Returns the input array without repeated terms (subroutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: sunique`  
!! `integer :: v(6)=[1, 2, 2, 3, 3, 4]`  
!! `integer, allocatable :: res(:)`  
!! `call sunique(v, res)`  
!! `print*, res`  
!! `end program`  
integer, intent(in)               :: v(:)   !! Input array.
integer, allocatable, intent(out) :: res(:) !! Array without repeated terms (not sorted).
integer, dimension(size(v,1)) :: tmp
integer :: ios, n, i
character(maxpath) :: err_msg

n = 0
do i = 1, size(v,1) !loop to find the elements
  if (all(tmp(1:n) /= v(i))) then
    n = n + 1
    tmp(n) = v(i)
  endif
enddo
allocate(res(n), stat=ios, errmsg=err_msg)
if (ios /= 0) call error('(module_set/sunique), unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end subroutine

!-----------------------------------------------------------------------
! unique
!-----------------------------------------------------------------------
function unique(v) result(res)
!! Returns the input array without repeated terms (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: unique`  
!! `integer :: v(6)=[1, 2, 2, 3, 3, 4]`  
!! `print*, unique(v)`  
!! `end program`  
!!
!! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
!! Thus, since [[module_set_bmod(module):unique(function)]]` returns an array, this code is correct:  
!! `integer, allocatable :: x(:)`  
!! `allocate(x, source=unique(...))`  
!! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
!! `integer, allocatable :: x(:)`  
!! `x = unique(...)`  
!! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
!! [[module_set_bmod(module):sunique(subroutine)]]`.  
integer, intent(in)  :: v(:)   !! Input array.
integer, allocatable :: res(:) !! Array without repeated terms (not sorted).
integer, dimension(size(v,1)) :: tmp
integer :: ios, n, i
character(maxpath) :: err_msg

call sunique(v, res)
end function

!-----------------------------------------------------------------------
! ssetdiff
!-----------------------------------------------------------------------
subroutine ssetdiff(v, w, res)
!! Returns the complement of `w` in `v` (subroutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: ssetdiff`  
!! `integer :: v(6) = [1, 2, 3, 4, 5, 6]`  
!! `integer :: w(4) = [1, 4, 2, 6]`  
!! `integer, allocatable :: res(:)`  
!! `call ssetdiff(v, w, res)`  
!! `print*, res`  
!! `end program`  
integer, intent(in)               :: v(:)   !! Universal set.
integer, intent(in)               :: w(:)   !! Set.
integer, allocatable, intent(out) :: res(:) !! Complement of `w` in `v` (not sorted, no repeated terms).
integer, dimension(size(v,1)) :: tmp
integer :: ios, n, i
character(maxpath) :: err_msg

n = 0
do i = 1, size(v,1) !loop to find the elements
  if (all(w /= v(i))) then
    if (any(tmp(1:n) == v(i))) cycle
    n = n + 1
    tmp(n) = v(i)
  endif
enddo
allocate(res(n), stat=ios, errmsg=err_msg)
if (ios /= 0) call error('module_set/ssetdiff), unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end subroutine

!-----------------------------------------------------------------------
! setdiff
!-----------------------------------------------------------------------
function setdiff(v, w) result(res)
!! Returns the complement of `w` in `v` (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: setdiff`  
!! `integer :: v(6) = [1, 2, 3, 4, 5, 6]`  
!! `integer :: w(4) = [1, 4, 2, 6]`  
!! `print*, setdiff(v, w)`  
!! `end program`  
!!
!! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
!! Thus, since [[module_set_bmod(module):setdiff(function)]]` returns an array, this code is correct:  
!! `integer, allocatable :: x(:)`  
!! `allocate(x, source=setdiff(...))`  
!! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
!! `integer, allocatable :: x(:)`  
!! `x = setdiff(...)`  
!! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
!! [[module_set_bmod(module):ssetdiff(subroutine)]]`.  
integer, intent(in)  :: v(:)   !! Universal set.
integer, intent(in)  :: w(:)   !! Set.
integer, allocatable :: res(:) !! Complement of `w` in `v` (not sorted, no repeated terms).
integer, dimension(size(v,1)) :: tmp
integer :: ios, n, i
character(maxpath) :: err_msg

call ssetdiff(v, w, res)
end function

end module

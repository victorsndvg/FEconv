module module_set
!-----------------------------------------------------------------------
! Module for set management
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 20/05/2013
!
! PUBLIC PROCEDURES:
!   unique: returns the values of the vector with no repetitions
!   sunique: subrutine for unique
!   setdiff: returns the values in one vector that are not in the other
!   ssetdiff: subrutine for setdiff
!
! REMARKS:
!   In gfortran 4.8 a generic interface cannot contain both functions and 
!   and subroutines. Thus,  unique and sunique must be different procedures
!   In gfortran 4.8, but not in ifort 13.0,  this statement is valid:
!
!     integer, allocatable :: x(:)
!     x = unique(...) 
! 
!   In ifort 13.0, but not in gfortran 4.8, this statement is valid: 
!
!     allocate(x, source=unique(...))
!
!   For compatibility between ifort 13.0 and gfortran 4.8, use sunique 
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error
implicit none

contains

!-----------------------------------------------------------------------
! unique (function): returns the same values as in v with no repetitions.
! REMARK: res will not be sorted
!-----------------------------------------------------------------------
function unique(v) result(res)
integer, intent(in) :: v(:)
integer, allocatable :: res(:)
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
if (ios /= 0) call error('unique, unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end function

!-----------------------------------------------------------------------
! sunique (subroutine): returns the same values as in v with no repetitions.
! REMARK: res will not be sorted
!-----------------------------------------------------------------------
subroutine sunique(v, res)
integer, intent(in) :: v(:)
integer, allocatable, intent(out) :: res(:)
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
if (ios /= 0) call error('unique, unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end subroutine

!-----------------------------------------------------------------------
! setdiff (function): returns the values in v that are not in w
! REMARK: res will not be sorted but its elements will be unique
!-----------------------------------------------------------------------
function setdiff(v, w) result(res)
integer, intent(in) :: v(:), w(:)
integer, allocatable :: res(:)
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
if (ios /= 0) call error('unique, unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end function

!-----------------------------------------------------------------------
! ssetdiff (subroutine): returns the values in v that are not in w
! REMARK: res will not be sorted but its elements will be unique
!-----------------------------------------------------------------------
subroutine ssetdiff(v, w, res)
integer, intent(in) :: v(:), w(:)
integer, allocatable, intent(out) :: res(:)
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
if (ios /= 0) call error('unique, unable to allocate output variable: '//trim(err_msg))
res = tmp(1:n)
end subroutine

end module

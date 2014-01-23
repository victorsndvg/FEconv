module module_alloc_char_r1
!-----------------------------------------------------------------------
! Module for memory allocation of character rank 1 allocatable arrays
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 02/03/2012
!
! PUBLIC PROCEDURES:
!   dealloc: dealloc memory
!   alloc: alloc memory
!   extend: extend the extension of an array
!   set: set a scalar or a vector in the array
!   insert: insert a scalar in the array
!   reduce: reduce the array
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error, info
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_vector_prv
private :: insert_prv
private :: search_multiple

!Interface
interface dealloc; module procedure    dealloc_prv; end interface
interface   alloc; module procedure      alloc_prv; end interface
interface  extend; module procedure     extend_prv; end interface
interface  reduce; module procedure     reduce_prv; end interface
interface     set; module procedure set_scalar_prv; end interface
interface     set; module procedure set_vector_prv; end interface
interface  insert; module procedure     insert_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory 
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
character(*), allocatable :: v(:)
integer :: res
character(maxpath) :: cad
  
if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_char_r1/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory 
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d)
character(*), allocatable :: v(:)
integer, intent(in) :: d
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d) then; v = ' '; return; end if
  call dealloc(v)
end if
allocate(v(d), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_char_r1/alloc) unable to allocate variable: '//trim(cad))
v = ' '
end subroutine

!-----------------------------------------------------------------------
! extend: extend the array to contain position (d)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d, fit)
character(*), allocatable     :: v(:)
integer, intent(in)           :: d !new dimension given by the user
logical, intent(in), optional :: fit
character(len(v)), allocatable :: temp(:)
integer :: res, s, ns 
character(maxpath) :: cad

if (.not. allocated(v)) then
  !DIMENSIONS
  if (present(fit)) then     
    if (fit) then; ns = d                        !we must fit to dimension given as argument
    else; ns = search_multiple(DEFAULT_ALLOC, d) !a multiple of DEFAULT_ALLOC must be taken as new dimension
    end if
  else; ns = d                                   !fit is not present, the same as if it where .true.
  end if
  !ALLOCATION
  allocate(v(ns), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_alloc_char_r1/extend) unable to allocate variable v: '//trim(cad))
  v = ' '
else !v is already allocated
  s = size(v,1)
  if (d > s) then !reallocation is mandatory
    !DIMENSIONS
    if (present(fit)) then 
      if (fit) then; ns = d            !we must fit to dimension given as argument, if necessary
      else; ns = search_multiple(s, d) !a multiple of the current size must be taken as new dimension
      end if
    else; ns = d                       !fit is not present, the same as if it where .true.
    end if
    !REALLOCATION
    allocate(temp(ns), stat = res, errmsg = cad)
    if (res /= 0) call error('(module_alloc_int_r2/extend) unable to allocate variable temp: '//trim(cad))
    temp(1:s)    = v
    temp(s+1:ns) = ' '
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d)
character(*), allocatable :: v(:)
integer, intent(in)       :: d
character(len(v)), allocatable :: temp(:)

!if (d == 0) then
!  call info('(module_alloc_char_r1/reduce) Given dimension is zero, variable will be deallocated')
!  call dealloc(v)
!else
if (.not. allocated(v)) then 
  call info('(module_alloc_char_r1/reduce) Variable not allocated'); return
end if
if (size(v,1) == d) return !v has the right size
if (size(v,1) <  d) then   !d is too large
  call info('(module_alloc_char_r1/reduce) Given dimension is too large to reduce'); return
end if
call alloc(temp, d)
temp = v(1:d)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set: set a scalar in the array
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d, fit)
character(*), allocatable     :: v(:)
character(*), intent(in)      :: val
integer, intent(in)           :: d
logical, intent(in), optional :: fit

call extend(v, d, fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! set: set a vector in the array
!-----------------------------------------------------------------------
subroutine set_vector_prv(v, val, d, fit)
character(*), allocatable     :: v(:)
character(*), intent(in)      :: val(:)
integer, intent(in)           :: d(:)
logical, intent(in), optional :: fit

call extend(v, maxval(d), fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! insert: insert a scalar in the array
!-----------------------------------------------------------------------
subroutine insert_prv(v, val, d, fit)
character(*), allocatable     :: v(:)
character(*), intent(in)      :: val
integer, intent(in)           :: d
logical, intent(in), optional :: fit

call extend(v, max(size(v,1)+1, d), fit)
v(d+1:size(v,1)) = v(d:size(v,1)-1)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! search_multiple: search the smallest value of 2 to the power of a that is bigger than b
! 2^n*a > b  <=>  n > log2(b/a)
!-----------------------------------------------------------------------
integer function search_multiple(a,b)
integer, intent(in) :: a, b

if (b > a) then 
  search_multiple = int(2**real(ceiling(log(real(b)/a)/log(2.)))*a)
else 
  search_multiple = a
end if
end function 

end module

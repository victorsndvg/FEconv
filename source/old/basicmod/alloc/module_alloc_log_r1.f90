module module_alloc_log_r1
!-----------------------------------------------------------------------
! Module for memory allocation of logical rank 1 allocatable arrays
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 02/10/2012
!
! PUBLIC PROCEDURES:
!   dealloc: dealloc memory
!   alloc: alloc memory
!   extend: extend the extension of an array
!   set: set a scalar or a vector in the array
!   add: add a scalar or a vector in the array
!   insert: insert a scalar in the array
!   reduce: reduce the array
!   find_first: funtion to find the first occurrence of a value in an array
!   find: funtion to find all the occurrences of a value in an array
!   sfind: subrutine for find
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_alloc_int_r1, only: alloc
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_vector_prv
private :: add_scalar_prv, add_vector_prv
private :: insert_prv
private :: find_first_prv, find_sca_prv, sfind_sca_prv
private :: search_multiple

!Interface
interface       dealloc; module procedure       dealloc_prv; end interface
interface         alloc; module procedure         alloc_prv; end interface
interface        extend; module procedure        extend_prv; end interface
interface        reduce; module procedure        reduce_prv; end interface
interface           set; module procedure    set_scalar_prv; end interface
interface           set; module procedure    set_vector_prv; end interface
interface           add; module procedure    add_scalar_prv; end interface
interface           add; module procedure    add_vector_prv; end interface
interface        insert; module procedure        insert_prv; end interface
interface    find_first; module procedure    find_first_prv; end interface
interface          find; module procedure      find_sca_prv; end interface
interface         sfind; module procedure     sfind_sca_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory 
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
logical, allocatable :: v(:)
integer :: res
character(maxpath) :: cad
  
if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_log_r1/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory 
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d)
logical, allocatable :: v(:)
integer, intent(in) :: d
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d) then; v = .false.; return; end if
  call dealloc(v)
end if
allocate(v(d), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_log_r1/alloc) unable to allocate variable: '//trim(cad))
v = .false.
end subroutine

!-----------------------------------------------------------------------
! extend: extend the array to contain position (d)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d, fit)
logical, allocatable          :: v(:), temp(:)
integer, intent(in)           :: d !new dimension given by the user
logical, intent(in), optional :: fit
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
  if (res /= 0) call error('(module_alloc_log_r1/extend) unable to allocate variable v: '//trim(cad))
  v = .false.
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
    temp(s+1:ns) = .false.
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d)
logical, allocatable :: v(:), temp(:)
integer, intent(in)  :: d

!if (d == 0) then
!  call info('(module_alloc_log_r1/reduce) Given dimension is zero, variable will be deallocated')
!  call dealloc(v)
!else
if (.not. allocated(v)) then 
  call info('(module_alloc_log_r1/reduce) Variable not allocated'); return
end if
if (size(v,1) == d) return !v has the right size
if (size(v,1) <  d) then   !d is too large
  call info('(module_alloc_log_r1/reduce) Given dimension is too large to reduce'); return
end if
call alloc(temp, d)
temp = v(1:d)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set: set a scalar in the array
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d, fit)
logical, allocatable :: v(:)
logical, intent(in)  :: val
integer, intent(in)  :: d
logical, intent(in), optional :: fit

call extend(v, d, fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! set: set a vector in the array
!-----------------------------------------------------------------------
subroutine set_vector_prv(v, val, d, fit)
logical, allocatable :: v(:)
logical, intent(in)  :: val(:)
integer, intent(in)  :: d(:)
logical, intent(in), optional :: fit

call extend(v, maxval(d), fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! add: add a scalar in the array
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d, fit)
logical, allocatable :: v(:)
logical, intent(in)  :: val
integer, intent(in)  :: d
logical, intent(in), optional :: fit

call extend(v, d, fit)
v(d) = v(d) .or. val
end subroutine

!-----------------------------------------------------------------------
! add: add a vector in the array
!-----------------------------------------------------------------------
subroutine add_vector_prv(v, val, d, fit)
logical, allocatable :: v(:)
logical, intent(in)  :: val(:)
integer, intent(in)  :: d(:)
logical, intent(in), optional :: fit

call extend(v, maxval(d), fit)
v(d) = v(d) .or. val
end subroutine

!-----------------------------------------------------------------------
! insert: insert a scalar in the array
!-----------------------------------------------------------------------
subroutine insert_prv(v, val, d, used, fit)
logical, allocatable :: v(:)
logical, intent(in)  :: val
integer, intent(in)  :: d
integer, intent(in), optional :: used
logical, intent(in), optional :: fit
integer :: s

if (present(used)) then; s = max(used+1, d)
else; s = max(size(v,1)+1, d)
end if

call extend(v, s, fit)
v(d+1:size(v,1)) = v(d:size(v,1)-1)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! find_first_prv: find the first occurrence of a scalar value in an array
!-----------------------------------------------------------------------
function find_first_prv(v, val) result(res)
logical,           intent(in) :: v(:) ! array
logical, optional, intent(in) :: val  ! value to search
integer                       :: res  ! result
integer :: i
logical :: lval

lval = .true.
if (present(val)) lval = val
res = 0
do i = 1, size(v,1)
  if (v(i) .eqv. lval) then
    res = i
    return
  end if  
end do
end function

!-----------------------------------------------------------------------
! find: find all the occurrences of scalar value in an array
!-----------------------------------------------------------------------
function find_sca_prv(v, val) result(res)
logical,             intent(in) :: v(:)   ! array
logical, optional,   intent(in) :: val    ! value to search
integer, allocatable            :: res(:) ! results
integer :: i, n, p
logical :: lval

lval = .true.
if (present(val)) lval = val
! allocate res
n = size(pack(v,v .eqv. lval),1)
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  if (v(i) .eqv. lval) then
    res(p) = i
    p = p+1
  end if  
  if (p > n) return
end do
end function

!-----------------------------------------------------------------------
! sfind (subroutine): find all the occurrences of scalar value in an array
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, res)
logical,              intent(in)  :: v(:)   ! array
logical,              intent(in)  :: val    ! value to search
integer, allocatable, intent(out) :: res(:) ! results
integer :: i, n, p

! allocate res
n = size(pack(v,v .eqv. val),1)
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  if (v(i) .eqv. val) then
    res(p) = i
    p = p+1
  end if  
  if (p > n) return
end do
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

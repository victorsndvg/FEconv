module module_alloc_real64_r1
!-----------------------------------------------------------------------
! Module for memory allocation of real64 rank 1 allocatable arrays
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 30/05/2013
!
! PUBLIC PROCEDURES:
!   dealloc: dealloc memory
!   alloc: alloc memory
!   extend: extend the extension of an array
!   set: set a scalar or a vector in the array
!   add: add a scalar or a vector in the array
!   insert: insert a scalar in the array
!   reduce: reduce the array
!   bsearch: search a value in a sorted vector
!   insert_sorted: insert a value in a sorted vector
!   sort: sort the array using the binary search
!   ssort: subrutine for sort
!   find_first: funtion to find the first occurrence of a value in an array
!   find: funtion to find all the occurrences of a value in an array
!     Private functions under 'find' are:
!     find_sca: find a scalar value in an array
!     find_vec: find an array of values in an array
!   sfind: subrutine for find
!
! REMARKS:
!   In gfortran 4.8 a generic interface cannot contain both functions and 
!   and subroutines. Thus,  find and sfind must be different procedures
!   In gfortran 4.8, but not in ifort 13.0,  this statement is valid:
!
!     real(real64), allocatable :: x(:)
!     x = find(...) 
! 
!   In ifort 13.0, but not in gfortran 4.8, this statement is valid: 
!
!     allocate(x, source=find(...))
!
!   For compatibility between ifort 13.0 and gfortran 4.8, use sfind
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
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
private :: insert_prv, bsearch_prv, insert_sorted_prv, sort_prv, ssort_prv
private :: find_first_prv, find_sca_prv, find_vec_prv, sfind_sca_prv, sfind_vec_prv
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
interface       bsearch; module procedure       bsearch_prv; end interface
interface insert_sorted; module procedure insert_sorted_prv; end interface
interface          sort; module procedure          sort_prv; end interface
interface         ssort; module procedure         ssort_prv; end interface
interface    find_first; module procedure    find_first_prv; end interface
interface          find; module procedure      find_sca_prv; end interface
interface          find; module procedure      find_vec_prv; end interface
interface         sfind; module procedure     sfind_sca_prv; end interface
interface         sfind; module procedure     sfind_vec_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory 
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
real(real64), allocatable :: v(:)
integer :: res
character(maxpath) :: cad
  
if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r1/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory 
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d)
real(real64), allocatable :: v(:)
integer, intent(in) :: d
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d) then; v = 0; return; end if
  call dealloc(v)
end if
allocate(v(d), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r1/alloc) unable to allocate variable: '//trim(cad))
v = 0
end subroutine

!-----------------------------------------------------------------------
! extend: extend the array to contain position (d)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d, fit)
real(real64), allocatable          :: v(:), temp(:)
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
  if (res /= 0) call error('(module_alloc_real64_r1/extend) unable to allocate variable v: '//trim(cad))
  v = 0
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
    if (res /= 0) call error('(module_alloc_real64_r2/extend) unable to allocate variable temp: '//trim(cad))
    temp(1:s)    = v
    temp(s+1:ns) = 0
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d)
real(real64), allocatable :: v(:), temp(:)
integer, intent(in)  :: d

!if (d == 0) then
!  call info('(module_alloc_real64_r1/reduce) Given dimension is zero, variable will be deallocated')
!  call dealloc(v)
!else
if (.not. allocated(v)) then 
  call info('(module_alloc_real64_r1/reduce) Variable not allocated'); return
end if
if (size(v,1) == d) return !v has the right size
if (size(v,1) <  d) then   !d is too large
  call info('(module_alloc_real64_r1/reduce) Given dimension is too large to reduce'); return
end if
call alloc(temp, d)
temp = v(1:d)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set: set a scalar in the array
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d, fit)
real(real64), allocatable :: v(:)
real(real64), intent(in)  :: val
integer, intent(in)  :: d
logical, intent(in), optional :: fit

call extend(v, d, fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! set: set a vector in the array
!-----------------------------------------------------------------------
subroutine set_vector_prv(v, val, d, fit)
real(real64), allocatable :: v(:)
real(real64), intent(in)  :: val(:)
integer, intent(in)  :: d(:)
logical, intent(in), optional :: fit

call extend(v, maxval(d), fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! add: add a scalar in the array
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d, fit)
real(real64), allocatable :: v(:)
real(real64), intent(in)  :: val
integer, intent(in)  :: d
logical, intent(in), optional :: fit

call extend(v, d, fit)
v(d) = v(d) + val
end subroutine

!-----------------------------------------------------------------------
! add: add a vector in the array
!-----------------------------------------------------------------------
subroutine add_vector_prv(v, val, d, fit)
real(real64), allocatable :: v(:)
real(real64), intent(in)  :: val(:)
integer, intent(in)  :: d(:)
logical, intent(in), optional :: fit

call extend(v, maxval(d), fit)
v(d) = v(d) + val
end subroutine

!-----------------------------------------------------------------------
! insert: insert a scalar in the array
!-----------------------------------------------------------------------
subroutine insert_prv(v, val, d, used, fit)
real(real64), allocatable :: v(:)
real(real64), intent(in)  :: val
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
! bsearch: If 'u' is in 'x' then its index is returned.
! Otherwise − ind is returned, where ind is the index where 'u'
! should be inserted in 'x'
!-----------------------------------------------------------------------
function bsearch_prv(x, u, n) result(pos)
real(real64), intent(in) :: x(:) ! array, where to search
real(real64), intent(in) :: u    ! value to search
integer, intent(in) :: n    ! number of components in the array
integer :: pos, i, j, k

pos = -1
if (size(x,1) <= 0 .or. n <= 0) return
if (u < x(1))  return

i=1; j=n
do
  k=(i+j)/2
  if (u < x(k)) then
    j=k
  else
    i=k
  end if
  if (i+1 >= j) exit
end do
if (u== x(i)) then
  pos = i
elseif ((i+1) <= n) then
  if (u == x(i+1)) then
    pos = i+1
  elseif (u < x(i+1)) then
    pos = -(i+1)
  else
    pos = -(i+2)
  end if  
else
  pos = -(i+1)
end if
end function

!-----------------------------------------------------------------------
! insert_sorted: insert a value 'val' in a sorted vector 'v'
!-----------------------------------------------------------------------
subroutine insert_sorted_prv(v, val, used, fit)
real(real64), allocatable             :: v(:) ! vector
real(real64),           intent(in)    :: val  ! value to insert
integer, optional,      intent(inout) :: used ! total number of elements
logical, optional,      intent(in)    :: fit  ! whether to fit or not
integer :: pos, n

!not allocated
if (.not. allocated(v)) then
  call set(v, val, 1, fit)
  if (present(used)) used = 1
  return
end if

!search among the existing elements
if (present(used)) then; n = used
else;                    n = size(v,1)
end if
pos = bsearch(v(1:n), val, n)
!insert and return
if (pos < 0) then
  call insert(v, val, -pos, n, fit) !if n = size(v), reallocation is made
  if (present(used)) used = n+1
  return
end if
end subroutine

!-----------------------------------------------------------------------
! sort (function): sort the array using the binary search
!
! WARNING: sort does not check whether v is allocated 
! (the reason is to operate also with non allocatable arrays)
!-----------------------------------------------------------------------
function sort_prv(v) result(res)
real(real64), dimension(:), intent(in) :: v
real(real64), allocatable :: res(:)
integer :: i, p

if (size(v,1) <= 0) return ! (gfortran 4.8) size is not reliable when array is not allocated
call alloc(res,size(v,1))
res = v
do i = 2, size(v,1)
  p = bsearch(res(1:i), v(i), i)
  p = abs(p)
  if (p < i) then
    res(p+1:i) = res(p:i-1)
    res(p) = v(i)
  end if
end do
end function

!-----------------------------------------------------------------------
! ssort (subroutine): sort the array using the binary search
!
! WARNING: the array must be allocatable in order to ensure a correct size with gfortran
!-----------------------------------------------------------------------
subroutine ssort_prv(v)
real(real64), allocatable, intent(inout) :: v(:)
real(real64) :: vi
integer :: i, p

if (.not. allocated(v)) return
if (size(v,1) <= 1) return
do i = 2, size(v,1)
  vi = v(i)
  p = bsearch(v(1:i), vi, i)
  p = abs(p)
  if (p < i) then
    v(p+1:i) = v(p:i-1)
    v(p) = vi
  end if
end do
end subroutine

!-----------------------------------------------------------------------
! find_first_prv: find the first occurrence of a scalar value in an array
!-----------------------------------------------------------------------
function find_first_prv(v, val) result(res)
real(real64),             intent(in) :: v(:) ! array
real(real64),             intent(in) :: val  ! value to search
integer                         :: res  ! result
integer :: i

res = 0
do i = 1, size(v,1)
  if (v(i) == val) then
    res = i
    return
  end if  
end do
end function

!-----------------------------------------------------------------------
! find: find all the occurrences of scalar value in an array
!-----------------------------------------------------------------------
function find_sca_prv(v, val) result(res)
real(real64),             intent(in) :: v(:)   ! array
real(real64),             intent(in) :: val    ! value to search
integer, allocatable            :: res(:) ! results
integer :: i, n, p

! allocate res
n = size(pack(v,v==val),1)
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  if (v(i) == val) then
    res(p) = i
    p = p+1
  end if  
  if (p > n) return
end do
end function

!-----------------------------------------------------------------------
! find: find all the occurrences of an array of values in an array
!-----------------------------------------------------------------------
function find_vec_prv(v, val) result(res)
real(real64),              intent(in)  :: v(:)   ! array
real(real64),              intent(in)  :: val(:) ! values to search
integer, allocatable              :: res(:) ! results
integer :: i, j, n, p

! allocate res
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v,v==val(j)),1)
end do  
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(val,1)
    if (v(i) == val(j)) then
      res(p) = i
      p = p+1
    end if  
  end do  
  if (p > n) return
end do
end function

!-----------------------------------------------------------------------
! sfind (subroutine): find all the occurrences of scalar value in an array
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, res)
real(real64),              intent(in)  :: v(:)   ! array
real(real64),              intent(in)  :: val    ! value to search
integer, allocatable, intent(out) :: res(:) ! results
integer :: i, n, p

! allocate res
n = size(pack(v,v==val),1)
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  if (v(i) == val) then
    res(p) = i
    p = p+1
  end if  
  if (p > n) return
end do
end subroutine

!-----------------------------------------------------------------------
! sfind (subroutine): find all the occurrences of an array of values in an array
!-----------------------------------------------------------------------
subroutine sfind_vec_prv(v, val, res)
real(real64),              intent(in)  :: v(:)   ! array
real(real64),              intent(in)  :: val(:) ! values to search
integer, allocatable, intent(out) :: res(:) ! results
integer :: i, j, n, p

! allocate res
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v,v==val(j)),1)
end do  
call alloc(res, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(val,1)
    if (v(i) == val(j)) then
      res(p) = i
      p = p+1
    end if  
  end do  
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

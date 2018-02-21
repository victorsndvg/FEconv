module module_alloc_real64_r1_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for memory allocation of real64 rank 1 allocatable arrays.  
!!
!! @note Procedures `set`, `add` and `insert` allocate or extend the array when it is necessary.  
!!
!! @note In GFortran 4.8 a generic interface cannot contain both functions and subroutines. 
!! This is why pairs like `[[module_alloc_real64_r1_bmod(module):find(interface)]]` and 
!! `[[module_alloc_real64_r1_bmod(module):sfind(interface)]]` have different public names.  
!!
!! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
!! Thus, since `[[module_alloc_real64_r1_bmod(module):find(interface)]]` returns an array, this code is correct:  
!! `integer, allocatable :: x(:)`  
!! `allocate(x, source=find(...))`  
!! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
!! `integer, allocatable :: x(:)`  
!! `x = find(...)`  
!! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
!! `[[module_alloc_real64_r1_bmod(module):sfind(interface)]]`. 
!! This warning also aplies to `sort`.  
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
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error, info
use module_alloc_common_bmod, only: DEFAULT_ALLOC, csize, search_multiple
use module_alloc_int_r1_bmod, only: alloc
implicit none

!Constants
private :: DEFAULT_ALLOC

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_vector_prv
private :: add_scalar_prv, add_vector_prv
private :: insert_prv, bsearch_prv, insert_sorted_prv, sort_prv, ssort_prv
private :: find_first_prv, find_sca_prv, find_vec_prv, sfind_sca_prv, sfind_vec_prv
private :: search_multiple

!Interfaces
interface dealloc
  !! Deallocates an array.  
  !!
  !! @note The advantages of using this procedure intead the intrinsic procedure `deallocate` are:  
  !!  - Testing the previous deallocation of `v` is automatically done.  
  !!  - Deallocation is called with error catching arguments and every error is shown.  
  module procedure dealloc_prv
end interface

interface alloc
  !! Allocates an array.  
  !!
  !! @note If `v` is allocated and its extension is `d`, allocation is not carried out; 
  !! otherwise, `v` is deallocated before allocation. 
  !! If there is not any error, `v` is set to `0`.  
  module procedure alloc_prv
end interface

interface extend
  !! Extends the array to contain position `d`.  
  !!
  !! @note If `v` is not allocated and  
  !! -`fit` is present and `.false.`, `v` will be allocated with extension `2^n DEFAULT_ALLOC`, where `n` is 
  !! the smallest exponent such that `2^n DEFAULT_ALLOC > d` (`DEFAULT_ALLOC` is a private parameter with value 1000).  
  !! - Otherwise, `v` is allocated with extension `d`.  
  !!
  !! @note If `v` is already allocated and  
  !! - `d` is larger than `size(v)`, `fit` is present and `.false.`, `v` will be reallocated with extension `2^n size(v)`, 
  !! where `n` is the smallest exponent such that `2^n size(v) > d`.  
  !! - Otherwise, `v` is reallocated with extension `d`.  
  !!
  !! @note Argument `fit` is `.true.` by default.  
  !!
  !! @note New positions created in the extension are set to 0.  
  module procedure extend_prv
end interface

interface reduce
  !! Reduces the extension of an array.  
  !!
  !! @note If `v` is not allocated or `d` is larger that its size, the procedure sends an info and returns. 
  !! If `v` is allocated and its extension is `d`, the procedure returns. 
  !! Otherwise, the procedure reduces the extension of `v` to `d`.  
  module procedure reduce_prv
end interface

interface set
  !! Sets a scalar o vector value in an array.  
  !!
  !! @note In the scalar case the function `extends` is initially called: `call extend(v,d,fit)`. Then `val` 
  !! is assigned to `v(d)`. 
  !! In the vector case the function `extends` is initially called: `call extend(v, maxval(d), fit)`. Then 
  !! `val` is assigned to `v(d)`.  
  !!
  !! @warning When `val` and `d` are arrays, they must have the same shape. This check is not done due to performance reasons.  
  module procedure set_scalar_prv
  module procedure set_vector_prv
end interface

interface add
  !! Adds a scalar or vector value to an array.  
  !!
  !! @note In the scalar case the function `extend` is initially called: `call extend(v, d, fit)`. Then `val` is added to `v(d)`. 
  !! In the vector case the function `extend` is initially called: `call extend(v, maxval(d), fit)`. Then `val` is added  
  !! to `v(d)`.
  !!
  !! @warning When `val` and `d` are arrays, they must have the same shape. This check is not done due to performance reasons.  
  module procedure add_scalar_prv
  module procedure add_vector_prv
end interface

interface insert
  !! Inserts a scalar value in an array.  
  !!
  !! @note If `used` is present, `v` is extended until `max(used+1,d)`; otherwise, it is extended until 
  !! `max(size(v,1)+1, d)`. 
  !! Then `v(d:size(v)-1)` are shifted one position to the right and `val` is assigned to `v(d)`.  
  !!
  !! @warning Note that when `used` is not present, __the array is effectively extended in every call__.  
  module procedure insert_prv
end interface

interface bsearch
  !! Returns the position of a value in an array or (minus) the position where the value should be inserted in.  
  !!
  !! @note If the value is in the array, then its position is returned; otherwise, -p is returned, where p is the position 
  !! where the value should be inserted, maintaining the array sorted.  
  !!
  !! @warning The array must be previously sorted in order to use `bsearch`.  
  module procedure bsearch_prv
end interface

interface insert_sorted
  !! Inserts a scalar value in a sorted array.  
  !!
  !! @note The value is searched in the array using `bsearch`. If the array is not allocated, the value is inserted in the first 
  !! position using `set`. 
  !! If `used` is present, the procedure only works with `v(1:used)`; otherwise, the full array is considered.  
  module procedure insert_sorted_prv
end interface

interface sort
  !! Sorts an array using binary search (function version).  
  !!
  !! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since 
  !! version 13.0). Thus, this code is correct:  
  !! `real(real64), allocatable :: x(:)`  
  !! `allocate(x, source=sort(...))`  
  !! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
  !! `real(real64), allocatable :: x(:)`  
  !! `x = sort(...)`  
  !! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version `ssort`.
  module procedure sort_prv
end interface

interface ssort
  !! Sorts the allocatable array using the binary search  (subrutine version).
  module procedure ssort_prv
end interface

interface find_first
  !! Finds the first occurrence of a scalar value in an array.
  module procedure find_first_prv
end interface

interface sfind
  !! Finds all the occurrences of a scalar value or an array of values in an array (subrutine version).
  module procedure sfind_sca_prv
  module procedure sfind_vec_prv
end interface

interface find
  !! Finds all the occurrences of a scalar value or an array of values in an array (function version).  
  !!
  !! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0).
  !! Thus, if `find` returns an array, this code is correct:  
  !! `integer, allocatable :: x(:)`  
  !! `allocate(x, source=find(...))`  
  !! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
  !! `integer, allocatable :: x(:)`  
  !! `x = find(...)`  
  !! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
  !! `sfind`.  
  module procedure find_sca_prv
  module procedure find_vec_prv
end interface

contains

!-----------------------------------------------------------------------
! dealloc_prv (dealloc)
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
!! Deallocates an real64 array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dealloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call dealloc(v)`  
!! `end program`  
real(real64), allocatable :: v(:) !! Array to deallocate.
integer :: res
character(maxpath) :: cad

if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r1/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc_prv (alloc)
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d)
!! Allocates an real64 array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: alloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call alloc(v, 20)`  
!! `end program`  
real(real64), allocatable :: v(:) !! Array to allocate.
integer, intent(in)       :: d    !! Array extension.
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d) then; v = 0._real64; return; end if
  call dealloc(v)
end if
allocate(v(d), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r1/alloc) unable to allocate variable: '//trim(cad))
v = 0._real64
end subroutine

!-----------------------------------------------------------------------
! extend_prv (extend)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d, fit)
!! Extends the array to contain position `d`.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: extend`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call extend(v, 40)`  
!! `end program`  
real(real64), allocatable     :: v(:)    !! Array.
integer, intent(in)           :: d       !! Position to be included in the extension.
logical, intent(in), optional :: fit     !! Whether to fit the array extension to the maximum position; by default, `.true.`
real(real64), allocatable :: temp(:)
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
  v = 0._real64
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
    if (res /= 0) call error('(module_alloc_real64_r1/extend) unable to allocate variable temp: '//trim(cad))
    temp(1:s)    = v
    temp(s+1:ns) = 0._real64
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce_prv (reduce)
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d)
!! Reduces the extension of an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, reduce`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call alloc(v, 40)`  
!! `call reduce(v, 20)`  
!! `end program`  
real(real64), allocatable :: v(:)   !! Array.
integer, intent(in)       :: d      !! Array extension.
real(real64), allocatable :: temp(:)

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
! set_scalar_prv (set)
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d, fit)
!! Sets a scalar in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call set(v, 40._real64, 10, fit=.false.)`  
!! `end program`  
real(real64), allocatable     :: v(:) !! Array.
real(real64),      intent(in) :: val  !! Value.
integer,           intent(in) :: d    !! Position.
logical, optional, intent(in) :: fit  !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, d, fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! set_vector_prv (set)
!-----------------------------------------------------------------------
subroutine set_vector_prv(v, val, d, fit)
!! Sets a vector in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call set(v, real([40, 60], real64), [10,13], fit=.false.)`  
!! `end program`  
real(real64), allocatable     :: v(:)   !! Array.
real(real64), intent(in)      :: val(:) !! Value.
integer, intent(in)           :: d(:)   !! Position.
logical, intent(in), optional :: fit    !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, maxval(d), fit)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! add_scalar_prv (add)
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d, fit)
!! Adds a scalar value in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call add(v, 40._real64, 10, fit=.false.)`  
!! `end program`  
real(real64), allocatable     :: v(:) !! Array.
real(real64), intent(in)      :: val  !! Value.
integer, intent(in)           :: d    !! Position.
logical, intent(in), optional :: fit  !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, d, fit)
v(d) = v(d) + val
end subroutine

!-----------------------------------------------------------------------
! add_vector_prv (add)
!-----------------------------------------------------------------------
subroutine add_vector_prv(v, val, d, fit)
!! Adds a vector in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call add(v, real([40, 60], real64), [10,13], fit=.false.)`  
!! `end program`  
real(real64), allocatable     :: v(:)   !! Array.
real(real64), intent(in)      :: val(:) !! Value.
integer, intent(in)           :: d(:)   !! Position.
logical, intent(in), optional :: fit    !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, maxval(d), fit)
v(d) = v(d) + val
end subroutine

!-----------------------------------------------------------------------
! insert_prv (insert)
!-----------------------------------------------------------------------
subroutine insert_prv(v, val, d, used, fit)
!! Inserts a scalar in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: insert`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call insert(v, 40._real64, 10, used=7, fit=.false.)`  
!! `end program`  
real(real64), allocatable     :: v(:) !! Array.
real(real64), intent(in)      :: val  !! Value.
integer, intent(in)           :: d    !! Position.
integer, intent(in), optional :: used !! Number of rows actually used in vector.
logical, intent(in), optional :: fit  !! Whether to fit the array extension to the maximum position; by default, `.true.`
integer :: s

if (present(used)) then; s = max(used+1, d)
else; s = max(size(v,1)+1, d)
end if

call extend(v, s, fit)
v(d+1:size(v,1)) = v(d:size(v,1)-1)
v(d) = val
end subroutine

!-----------------------------------------------------------------------
! bsearch_prv (bsearch)
!-----------------------------------------------------------------------
function bsearch_prv(x, u, n) result(pos)
!! Returns the position of a value in an array or (minus) the position where the value should be inserted in.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: bsearch`  
!! `implicit none`  
!! `real(real64) :: pos, v(5) = real([1, 3, 4, 8, 10], real64)`  
!! `pos = bsearch(v, 2._real64, 5)`  
!! `end program`  
real(real64), intent(in) :: x(:) !! Array.
real(real64), intent(in) :: u    !! Value.
integer, intent(in)      :: n    !! Number of positions actually used in the array.
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
if (abs(u-x(i)) <= epsilon(1._real64)) then
  pos = i
elseif ((i+1) <= n) then
  if (abs(u-x(i+1)) <= epsilon(1._real64)) then
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
! insert_sorted_prv (insert_sorted)
!-----------------------------------------------------------------------
subroutine insert_sorted_prv(v, val, used, fit)
!! Inserts a scalar value in a sorted array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: insert_sorted`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call insert_sorted(v, 40._real64, used=7, fit=.false.)`  
!! `end program`  
real(real64), allocatable        :: v(:) !! Array.
real(real64),      intent(in)    :: val  !! Value.
integer, optional, intent(inout) :: used !! Number of rows actually used in the array.
logical, optional, intent(in)    :: fit  !! Whether to fit the array extension to the maximum position; by default, `.true.`
integer :: pos
integer :: n

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
! sort_prv (sort)
!
! Note: sort cannot be written in terms of ssort because the latter changes the argument
!-----------------------------------------------------------------------
function sort_prv(v) result(res)
!! Sorts the array using the binary search (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: sort`  
!! `implicit none`  
!! `real(real64) :: v(5) = real([2, 3, 4, 1, 9], real64)`  
!! `real(real64), allocatable :: res(:)`  
!! `res = sort(v)                 ! it works in GFortran 4.8, but not in ifort   13.0`  
!! `allocate(res, source=sort(v)) ! it works in ifort   13.0, but not in GFortran 4.8`  
!! `end program`  
real(real64), intent(in)  :: v(:)   !! Array.
real(real64), allocatable :: res(:) !! Sorted array.
integer :: i, p

if (csize(dbl1=v, d=1) <= 0) return
call alloc(res, size(v,1))
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
! ssort_prv (ssort)
!-----------------------------------------------------------------------
subroutine ssort_prv(v)
!! Sorts the allocatable array using the binary search (subrutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, ssort`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call alloc(v, 3)`  
!! `v = real([1, 3, 2], real64)`  
!! `call ssort(v)`  
!! `end program`  
real(real64), intent(inout) :: v(:) !! Array.
real(real64) :: vi
integer :: i, p

if (csize(dbl1=v, d=1) <= 1) return
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
! find_first_prv (find_first)
!-----------------------------------------------------------------------
function find_first_prv(v, val) result(res)
!! Finds the first occurrence of a scalar value in an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, find_first`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:)`  
!! `call alloc(v, 6)`  
!! `v = real([1, 2, 4, 5, 4, 6], real64)`  
!! `print*, find_first(v, 4._real64)`  
!! `end program`  
real(real64), intent(in) :: v(:) !! Array.
real(real64), intent(in) :: val  !! Value to search.
integer                  :: res  !! Position of the first occurrence; return 0 if nothing was found.
integer :: i

res = 0
if (csize(dbl1=v, d=1) <= 0) return
do i = 1, size(v,1)
  if (abs(val-v(i)) <= epsilon(1._real64)) then
    res = i
    return
  end if
end do
end function

!-----------------------------------------------------------------------
! sfind_sca_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, res)
!! Finds all the occurrences of scalar value in an array (subrutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, sfind`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:), pos(:)`  
!! `call alloc(v, 6)`  
!! `v = real([1, 2, 4, 5, 4, 6], real64)`  
!! `call sfind(v, 4, pos)`  
!! `end program`  
real(real64),         intent(in)  :: v(:)   !! Array.
real(real64),         intent(in)  :: val    !! Value to search.
integer, allocatable, intent(out) :: res(:) !! Array of positions; returns a zero size array if nothing was found.
integer :: i, n, p

! allocate res
if (csize(dbl1=v, d=1) <= 0) return
n = size(pack(v, abs(v-val) <= epsilon(1._real64)), 1)
call alloc(res, n)
if (n <= 0) return
! find positions
p = 1
do i = 1, size(v,1)
  if (abs(v(i)-val) <= epsilon(1._real64)) then
    res(p) = i
    p = p+1
  end if
  if (p > n) return
end do
end subroutine

!-----------------------------------------------------------------------
! sfind_vec_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_vec_prv(v, val, res)
!! Finds all the occurrences of an array of values in an array (subrutine version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, sfind`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:), pos(:)`  
!! `call alloc(v, 6)`  
!! `v = real([1, 6, 4, 5, 4, 6], real64)`  
!! `call sfind(v, [4,6], pos)`  
!! `end program`  
real(real64),         intent(in)  :: v(:)   !! Array.
real(real64),         intent(in)  :: val(:) !! Values to search.
integer, allocatable, intent(out) :: res(:) !! Array of positions; a zero-sized array if nothing was found.
integer :: i, j, n, p

! allocate res
if (csize(dbl1=v, d=1) <= 0) return
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v,abs(v-val(j)) <= epsilon(1._real64)),1)
end do
call alloc(res, n)
if (n <= 0) return
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(val,1)
    if (abs(v(i)-val(j)) <= epsilon(1._real64)) then
      res(p) = i
      p = p+1
      exit
    end if
  end do
  if (p > n) return
end do
end subroutine

!-----------------------------------------------------------------------
! find_sca_prv (find)
!-----------------------------------------------------------------------
function find_sca_prv(v, val) result(res)
!! Finds all the occurrences of scalar value in an array (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, find`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:), pos(:)`  
!! `call alloc(v, 6)`  
!! `v = real([1, 6, 4, 5, 4, 6], real64)`  
!! `allocate(pos, source=find(v, 4))`  
!! `! the latter is not valid in Gfortan 4.9 and below; write the following code or use sfind instead:`  
!! `! pos = find(v, 4)`  
!! `end program`  
real(real64), intent(in) :: v(:)   !! Array.
real(real64), intent(in) :: val    !! Value to search.
integer, allocatable     :: res(:) !! Array of positions; a zero-sized array if nothing was found.

call sfind(v, val, res)
end function

!-----------------------------------------------------------------------
! find_vec_prv (find)
!-----------------------------------------------------------------------
function find_vec_prv(v, val) result(res)
!! Finds all the occurrences of an array of values in an array (function version).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: find`  
!! `implicit none`  
!! `integer :: v(6) = real([1, 6, 4, 5, 4, 6], real64)`  
!! `real(real64), allocatable :: pos(:)`  
!! `allocate(pos, source=find(v, [4,6]))`  
!! `! the latter is not valid in Gfortan 4.9 and below; write the following code or use sfind instead:`  
!! `! pos = find(v, [4,6])`  
!! `end program`  
real(real64), intent(in)  :: v(:)   !! Array.
real(real64), intent(in)  :: val(:) !! Values to search.
integer, allocatable      :: res(:) !! Array of positions; a zero-sized array if nothing was found.

call sfind(v, val, res)
end function
end module

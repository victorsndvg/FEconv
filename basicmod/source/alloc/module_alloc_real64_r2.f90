module module_alloc_real64_r2_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for memory allocation of real64 rank 2 allocatable arrays.  
!!
!! @note Procedures `set`, `add` and `insert` allocate or extend the array when it is necessary.  
!!
!! @note In GFortran 4.8 a generic interface cannot contain both functions and subroutines. 
!! This is why pairs like `[[module_alloc_real64_r2_bmod(module):find(interface)]]` and 
!! `[[module_alloc_real64_r2_bmod(module):sfind(interface)]]` have different public names.  
!!
!! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
!! Thus, since `[[module_alloc_real64_r2_bmod(module):find(interface)]]` returns an array, this code is correct:  
!! `integer, allocatable :: x(:,:)`  
!! `allocate(x, source=find(...))`  
!! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
!! `integer, allocatable :: x(:,:)`  
!! `x = find(...)`  
!! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
!! `[[module_alloc_real64_r2_bmod(module):sfind(interface)]]`.  
!
! PUBLIC PROCEDURES:
!   dealloc: deallocates memory
!   alloc: allocates memory
!   extend: extends the shape of the array
!   set: sets a scalar, section or matrix value in the array
!   add: adds a scalar, section or matrix value in the array
!   insert: inserts a section (row/column) in the array
!   reduce: reduces the shape of the array
!   insert_sorted: inserts a section (row/column) in a row-sorted array
!   sfind: finds all the occurrences of a scalar value or a vector of values in the array (subroutine version)
!   find:  finds all the occurrences of a scalar value or a vector of values in the array (function version)
!   find_sorted: finds the position of a section (row/column) in a (row/column)-sorted array
!
! This module has an older "check" version, with additional checks for the arguments. See source files in folder
! alloc/test/test_check.
!-----------------------------------------------------------------------
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error, info
use module_alloc_int_r2_bmod, only: alloc
use module_alloc_real64_r1_bmod, only: real64, alloc, bsearch
use module_alloc_common_bmod, only: DEFAULT_ALLOC, csize, search_multiple
implicit none

!Constants
private :: DEFAULT_ALLOC

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_section_prv, set_matrix_prv
private :: add_scalar_prv, add_section_prv, add_matrix_prv
private :: insert_section_prv, insert_section_sorted_prv
private :: find_sca_prv, find_vec_prv, sfind_sca_prv, sfind_vec_prv
private :: find_section_sorted_prv
private :: search_multiple

!Interfaces
interface dealloc
  !! Deallocates a real64 array of rank 2.  
  !!
  !! @note The advantages of using this procedure intead the intrinsic procedure deallocate are:  
  !!  - Testing the previous deallocation of `v` is automatically done.  
  !!  - Deallocation is called with error catching arguments and every error is shown.  
  module procedure dealloc_prv
end interface

interface alloc
  !! Allocates a real64 array of rank 2.  
  !!
  !! @note If `v` is allocated and its extensions are (`d1`,`d2`), allocation is not carried out; otherwise, `v` is deallocated 
  !! before allocation. 
  !! If there is not any error, `v` is set to 0.  
  module procedure alloc_prv
end interface

interface extend
  !! Extends the array to contain position (`d1`,`d2`). 
  !!
  !! @note Consider the case of `v` being not allocated. For each `i=1,2`:  
  !! - If `fit(i)` is present and `.false.`, `v` will be allocated with extension `2^n·DEFAULT_ALLOC` in the dimension `i`, where 
  !! `n` is the smallest exponent such that `2^n·DEFAULT_ALLOC > d(i)` (`DEFAULT_ALLOC` is a private parameter with value `1000`).
  !! - Otherwise, `v` is allocated with extension `d(i)` in the dimension `i`.  
  !!
  !! @note Consider the case of `v` being already allocated. For each `i=1,2`, if some `d(i)` is larger than `size(v, i)`, then  
  !! - If `fit(i)` is present and `.false.`, `v` will be reallocated with extension `2^n·size(v,i)` in the dimension `i`, where `n`
  !! is the smallest exponent such that `2^n·size(v,i) > d(i)`.  
  !! - Otherwise, `v` is reallocated with extension `max(size(v,i), d(i))` in the dimension `i`.  
  !!
  !! @note Argument `fit` is `.true.` by default.  
  !!
  !! @note New positions created in the extension are set to 0.  
  module procedure extend_prv
end interface

interface reduce
  !! Reduces the extension of an array.  
  !!
  !! @note If `v` is not allocated or some `d` is larger that its size, the procedure sends an info and returns. 
  !! If `v` is allocated and its extension in every dimension is its `d`, the procedure returns. 
  !! Otherwise, the procedure reduces the extension of `v` in each dimension to `d`.  
  module procedure reduce_prv
end interface

interface set
  !! Sets a scalar, section or matrix value in an array.  
  !!
  !! @note Function `extend` is initially called and then `val` is assigned.  
  !!
  !! @warning In the section and matrix case `val` must have the proper shape, although this check is not done due to performance 
  !! reasons.  
  module procedure set_scalar_prv
  module procedure set_section_prv
  module procedure set_matrix_prv
end interface

interface add
  !! Adds a scalar, section matrix value to an array.  
  !!
  !! @note Function `extend` is initially called and then `val` is assigned.  
  !!
  !! @warning In the section and matrix case `val` must have the proper shape, although this check is not done due to performance 
  !! reasons.  
  module procedure add_scalar_prv
  module procedure add_section_prv
  module procedure add_matrix_prv
end interface

interface insert
  !! Inserts a section (row/column) in the array.  
  !!
  !! @note The function `extend` is initially called. If `used` is present, `v` is extended in the given dimension until 
  !! `max(used+1, d)`; otherwise, it is extended until `max(size(v,dimen)+1` in the given dimension. Finally, sections after the 
  !! inserted one are shifted one position down and `val` is inserted.  
  !!
  !! @warning When `used` is not present, __the array is  efectively extended in every call__.  
  module procedure insert_section_prv
end interface

interface insert_sorted
  !! Inserts a section (row/column) in a sorted array.  
  !!
  !! @warning When `used` is not present, and the section must be inserted, __the array is efectively extended__.  
  module procedure insert_section_sorted_prv
end interface

interface sfind
  !! Find all the occurrences of a scalar value or an array of values in an array (subroutine version).
  module procedure sfind_sca_prv
  module procedure sfind_vec_prv
end interface

interface find
  !! Find all the occurrences of a scalar value or an array of values in an array (function version).  
  !!
  !! @warning The Fortran 2003 standard implements the clause `source` to allocate an array (included in ifort since version 13.0). 
  !! Thus, if `find` returns an array, this code is correct:  
  !! `integer, allocatable :: x(:,:)`  
  !! `allocate(x, source=find(...))`  
  !! GFortran 4.9 and below does not support this clause and implements an alternative non-standard way:  
  !! `integer, allocatable :: x(:,:)`  
  !! `x = find(...)`  
  !! If you use GFortran 4.9 or below, we do not recommend the latter non-standard code; use instead the subroutine version, 
  !! `sfind`.  
  module procedure find_sca_prv
  module procedure find_vec_prv
end interface

interface find_sorted
  !! Finds the position of a section (row/column) in a (row/column)-sorted array.
  module procedure find_section_sorted_prv
end interface

contains

!-----------------------------------------------------------------------
! dealloc_prv (dealloc)
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
!! Deallocates a real64 array of rank 2.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dealloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call dealloc(v)`  
!! `end program`  
real(real64), allocatable :: v(:,:)  !! Array to deallocate.
integer :: res
character(maxpath) :: cad

if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r2/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc_prv (alloc)
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d1, d2)
!! Allocates a real64 array of rank 2.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call alloc(v, 20, 30)`  
!! `end program`  
real(real64), allocatable :: v(:,:) !! Array to allocate.
integer, intent(in)       :: d1     !! Array extension in the first dimension.
integer, intent(in)       :: d2     !! Array extension in the second dimension.
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d1 .and. size(v,2) == d2) then; v = 0._real64; return; end if
  call dealloc(v)
end if
allocate(v(d1, d2), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r2/alloc) unable to allocate variable: '//trim(cad))
v = 0._real64
end subroutine

!-----------------------------------------------------------------------
! extend_prv (extend)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d1, d2, fit)
!! Extends the array to contain position (`d1`,`d2`).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, extend`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call alloc(v, 40, 20)`  
!! `call extend(v, 2, 40, fit=[.false.,.true.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:) !! Array to change dimensions.
integer, intent(in)           :: d1     !! Position in the first dimension to be included in the extension.
integer, intent(in)           :: d2     !! Position in the second dimension to be included in the extension.
logical, intent(in), optional :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`
real(real64), allocatable :: temp(:,:)
integer :: res, s1, s2, ns1, ns2
character(maxpath) :: cad

if (.not. allocated(v)) then
  !DIMENSIONS
  if (present(fit)) then
    if (fit(1)) then; ns1 = d1                     !we must fit to dimension given as argument
    else; ns1 = search_multiple(DEFAULT_ALLOC, d1) !a multiple of DEFAULT_ALLOC must be taken as new dimension
    end if
    if (fit(2)) then; ns2 = d2                     !we must fit to dimension given as argument
    else; ns2 = search_multiple(DEFAULT_ALLOC, d2) !a multiple of DEFAULT_ALLOC must be taken as new dimension
    end if
  else; ns1 = d1; ns2 = d2                         !fit is not present, the same as if it where .true.
  end if
  !ALLOCATION
  allocate(v(ns1, ns2), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_alloc_real64_r2/extend) unable to allocate variable v: '//trim(cad))
  v = 0._real64
else !v is already allocated
  s1 = size(v,1); s2 = size(v,2)
  if (d1 > s1 .or. d2 > s2) then !reallocation is mandatory
    !DIMENSIONS
    if (present(fit)) then
      if (fit(1)) then; ns1 = max(s1, d1)          !we must fit to dimension given as argument, if necessary
      else; ns1 = search_multiple(s1, d1)          !a multiple of the current size must be taken as new dimension
      end if
      if (fit(2)) then; ns2 = max(s2, d2)          !we must fit to dimension given as argument, if necessary
      else; ns2 = search_multiple(s2, d2)          !a multiple of the current size must be taken as new dimension
      end if
    else; ns1 = max(s1, d1); ns2 = max(s2, d2)     !fit is not present, the same as if it where .true.
    end if
    !REALLOCATION
    allocate(temp(ns1, ns2), stat = res, errmsg = cad)
    if (res /= 0) call error('(module_alloc_real64_r2/extend) unable to allocate variable temp: '//trim(cad))
    temp = 0._real64
    temp(1:s1,1:s2) = v
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce_prv (reduce)
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d1, d2)
!! Reduces the extension of an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, reduce`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call alloc(v, 40, 40)`  
!! `call reduce(v, 20, 20)`  
!! `end program`  
real(real64), allocatable :: v(:,:)   !! Array to reduce.
integer, intent(in)       :: d1       !! Array extension in the first dimension.
integer, intent(in)       :: d2       !! Array extension in the second dimension.
real(real64), allocatable :: temp(:,:)

!if (d1 == 0 .or. d2 == 0) then
!  call info('(module_alloc_real64_r2/reduce) Some given dimension is zero, variable will be deallocated')
!  call dealloc(v)
!else
if (.not. allocated(v)) then
  call info('(module_alloc_real64_r2/reduce) Variable not allocated'); return
end if
if (size(v,1) == d1 .and. size(v,2) == d2) return !rows and cols have the right size
if (size(v,1) <  d1 .or.  size(v,2) <  d2) then   !rows or cols are too large
  call info('(module_alloc_real64_r2/reduce) Some given dimension is too large to reduce'); return
end if
call alloc(temp, d1, d2)
temp(1:d1, 1:d2) = v(1:d1, 1:d2)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set_scalar_prv (set)
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d1, d2, fit)
!! Sets a scalar value in an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call set(v, 25._real64, 10, 20, fit=[.true.,.false.])`  
!! `end program`  
real(real64), allocatable :: v(:,:)     !! Array.
real(real64), intent(in)  :: val        !! Value.
integer, intent(in)       :: d1         !! Position in the first dimension.
integer, intent(in)       :: d2         !! Position in the second dimension.
logical, intent(in), optional :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, d1, d2, fit)
v(d1, d2) = val
end subroutine

!-----------------------------------------------------------------------
! set_section_prv (set)
!-----------------------------------------------------------------------
subroutine set_section_prv(dimen, v, val, d, fit)
!! Sets a section (row/column) in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call set(1, v, real([25,35,45,55],real64), 10, fit=[.true.,.false.])`  
!! `end program`  
integer,           intent(in) :: dimen  !! Dimension to work with.
real(real64), allocatable     :: v(:,:) !! Array.
real(real64),      intent(in) :: val(:) !! Section values.
integer,           intent(in) :: d      !! Position in the given dimension of the section to insert.
logical, optional, intent(in) :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`

select case(dimen)
case (1) ! Row
  call extend(v, d, size(val,1), fit)
  v(d, 1:size(val,1)) = val
case (2) ! Column
  call extend(v, size(val,1), d, fit)
  v(1:size(val,1), d) = val
case default
  call error('(module_alloc_real64_r2/set_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! set_matrix_prv (set)
!-----------------------------------------------------------------------
subroutine set_matrix_prv(v, val, d1, d2, fit)
!! Sets a matrix in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call set(v, reshape(real([25,35,45,55],real64),[2,2]), [10,20], [1,3], fit=[.true.,.false.])`  
!! `end program`  
real(real64), allocatable :: v(:,:)     !! Array.
real(real64), intent(in)  :: val(:,:)   !! Values.
integer, intent(in)       :: d1(:)      !! Positions in the first dimension.
integer, intent(in)       :: d2(:)      !! Positions in the second dimension.
logical, intent(in), optional :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, maxval(d1), maxval(d2), fit)
v(d1, d2) = val
end subroutine

!-----------------------------------------------------------------------
! add_scalar_prv (add)
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d1, d2, fit)
!! Adds a scalar in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call add(v, 25._real64, 10, 20, fit=[.true.,.false.])`  
!! `end program`  
real(real64), allocatable :: v(:,:)     !! Array.
real(real64), intent(in)  :: val        !! Value.
integer, intent(in)       :: d1         !! Position in the first dimension.
integer, intent(in)       :: d2         !! Position in the second dimension.
logical, intent(in), optional :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, d1, d2, fit)
v(d1, d2) = v(d1, d2) + val
end subroutine

!-----------------------------------------------------------------------
! add_section_prv (add)
!-----------------------------------------------------------------------
subroutine add_section_prv(dimen, v, val, d, fit)
!! Adds a section (row/column) in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call add(1, v, real([25,35,45,55],real64), 10, fit=[.true.,.false.])`  
!! `end program`  
integer,           intent(in) :: dimen  !! Dimension to work with.
real(real64), allocatable     :: v(:,:) !! Array.
real(real64),      intent(in) :: val(:) !! Section values.
integer,           intent(in) :: d      !! Position in the given dimension of the section to insert.
logical, optional, intent(in) :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`

select case(dimen)
case (1) !Row
  call extend(v, d, size(val,1), fit)
  v(d, 1:size(val,1)) = v(d, 1:size(val,1)) + val
case (2) !Column
  call extend(v, size(val,1), d, fit)
  v(1:size(val,1), d) = v(1:size(val,1), d) + val
case default
  call error('(module_alloc_real64_r2/set_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! add_matrix_prv (add)
!-----------------------------------------------------------------------
subroutine add_matrix_prv(v, val, d1, d2, fit)
!! Adds a matrix in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call add(v, reshape(real([25,35,45,55],real64),[2,2]), [10,20], [1,3], fit=[.true.,.false.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:)   !! Array.
real(real64),      intent(in) :: val(:,:) !! Values.
integer,           intent(in) :: d1(:)    !! Positions in the first dimension.
integer,           intent(in) :: d2(:)    !! Positions in the second dimension.
logical, optional, intent(in) :: fit(2)   !! Whether to fit the array extension to the maximum position; by default, `.true.`

call extend(v, maxval(d1), maxval(d2), fit)
v(d1, d2) = v(d1, d2) + val
end subroutine

!-----------------------------------------------------------------------
! insert_section_prv (insert)
!-----------------------------------------------------------------------
subroutine insert_section_prv(dimen, v, val, d, used, fit)
!! Inserts a row /column in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: insert`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `call insert(1, v, real([25,35,45,55],real64), 10, used=3, fit=[.true.,.false.])`  
!! `end program`  
integer,           intent(in) :: dimen  !! Dimension to work with.
real(real64), allocatable     :: v(:,:) !! Array.
real(real64),      intent(in) :: val(:) !! Section values.
integer,           intent(in) :: d      !! Position in the given dimension of the section to insert.
integer, optional, intent(in) :: used   !! Number of sections actually used in vector.
logical, optional, intent(in) :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`
integer :: s

if (dimen == 1 .or. dimen == 2) then
  if (present(used)) then; s = max(used+1, d)
  else; s = max(size(v,dimen)+1, d)
  end if
end if
select case(dimen)
case (1) ! Row
  call extend(v, s, size(val,1), fit)
  v(d+1:size(v,1),           :)                       = v(d:size(v,dimen)-1, :)
  v(d,                       1:size(val,1))           = val
  v(d,                       size(val,1)+1:size(v,2)) = 0._real64
case (2) ! Column
  call extend(v, size(val,1), s, fit)
  v(:,                       d+1:size(v,2))           = v(:, d:size(v,dimen)-1)
  v(1:size(val,1),           d)                       = val
  v(size(val,1)+1:size(v,1), d)                       = 0._real64
case default
  call error('(module_alloc_real64_r2/insert_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! insert_section_sorted_prv (insert_sorted)
!-----------------------------------------------------------------------
subroutine insert_section_sorted_prv(dimen, v, val, used, fit, pos)
!! Inserts a section (row/column) in a (row/column)-sorted array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: alloc, insert_sorted`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `integer :: i`  
!! `call alloc(v, 4, 4)`  
!! `v = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `call insert_sorted(1, v, [(3._real64*i, i=1,4)])`  
!! `end program`  
integer,           intent(in)    :: dimen  !! Dimension to work with.
real(real64), allocatable        :: v(:,:) !! Array.
real(real64),      intent(in)    :: val(:) !! Section values.
integer, optional, intent(inout) :: used   !! Number of sections actually used in vector.
logical, optional, intent(in)    :: fit(2) !! Whether to fit the array extension to the maximum position; by default, `.true.`
integer, optional, intent(out)   :: pos    !! Position in the given dimension of the inserted section.
integer :: n, a, b, anew, bnew, i, j,indx

!v not allocated
if (.not. allocated(v)) then
  indx = -1
  call set(dimen, v, val, -indx, fit)
  if (present(used)) used = 1
  if (present(pos)) pos = indx
  return
end if
!number of existing rows
if (present(used)) then; n = used
else;                    n = size(v,dimen)
end if
select case(dimen)
case (1) !Row
  !search among the first vertices
  indx = bsearch(v(1:n,1), val(1), n)
  if (indx < 0) then
    !insert and return
    call insert(1, v, val, -indx, used=n, fit=fit)
    if (present(used)) used = n+1
    if (present(pos)) pos = indx
    return
  end if
  a=1; b = n
  do j = 2, size(val,1)
    !determine the left extreme of the interval where to search the j-th vertex
    anew = indx
    do i = indx-1, a, -1
      if (v(i,j-1) /= v(indx,j-1)) exit
      anew = i
    end do
    !determine the right extreme of the interval where to search the j-th vertex
    bnew = indx
    do i = indx+1, b
      if (v(i,j-1) /= v(indx,j-1)) exit
      bnew = i
    end do
    a = anew; b = bnew
    indx = bsearch(v(a:b,j), val(j), b-a+1)
    if (indx < 0) then
      indx = indx - a+1
      !insert and return
      call insert(1, v, val, -indx, used=n, fit=fit)
      if (present(used)) used = n+1
      if (present(pos)) pos = indx
      return
    else
      indx = indx + a-1
    end if
  end do
case (2) !Column
  !search among the first vertices
  indx = bsearch(v(1,1:n), val(1), n)
  if (indx < 0) then
    !insert and return
    call insert(2, v, val, -indx, used=n, fit=fit)
    if (present(used)) used = n+1
    if (present(pos)) pos = indx
    return
  end if
  a=1; b = n
  do j = 2, size(val,1)
    !determine the left extreme of the interval where to search the j-th vertex
    anew = indx
    do i = indx-1, a, -1
      if (v(j-1,i) /= v(j-1,indx)) exit
      anew = i
    end do
    !determine the right extreme of the interval where to search the j-th vertex
    bnew = indx
    do i = indx+1, b
      if (v(j-1,i) /= v(j-1,indx)) exit
      bnew = i
    end do
    a = anew; b = bnew
    indx = bsearch(v(j,a:b), val(j), b-a+1)
    if (indx < 0) then
      indx = indx - a+1
      !insert and return
      call insert(2, v, val, -indx, used=n, fit=fit)
      if (present(used)) used = n+1
      if (present(pos)) pos = indx
      return
    else
      indx = indx + a-1
    end if
  end do
case default
  call error('(module_alloc_real64_r2/insert_section_sorted_prv) 1st argument dimen must be 1 or 2.')
end select
if (present(pos)) pos = indx
end subroutine

!-----------------------------------------------------------------------
! sfind_sca_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, res)
!! Finds all the occurrences of scalar value in a rank 2 array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: sfind`  
!! `implicit none`  
!! `integer :: i`  
!! `real(real64) :: v(4,4) = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `integer, allocatable :: res(:,:)`  
!! `call sfind(v, 2._real64, res)`  
!! `end program`  
real(real64),         intent(in)  :: v(:,:)   !! Array.
real(real64),         intent(in)  :: val      !! Value to search.
integer, allocatable, intent(out) :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.
integer :: i, j, n, p

if (csize(dbl2=v, d=1) <= 0 .or. csize(dbl2=v, d=2) <= 0) return
! allocate row, col
n = size(pack(v, abs(v-val) <= epsilon(1._real64)),1)
call alloc(res, 2, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(v,2)
    if ( abs(v(i,j)-val) <= epsilon(1._real64)) then
      res(1,p) = i
      res(2,p) = j
      p = p+1
    end if
    if (p > n) return
  end do
end do
end subroutine

!-----------------------------------------------------------------------
! sfind_vec_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_vec_prv(v, val, res)
!! Finds all the occurrences of the elements of a rank 1 array in a rank 2 array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: alloc, sfind`  
!! `implicit none`  
!! `integer:: i`  
!! `real(real64) :: v(4,4) = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `integer, allocatable :: res(:,:)`  
!! `call sfind(v, [(2._real64*i, i=1,8)], res)`  
!! `end program`  
real(real64),         intent(in)  :: v(:,:)   !! Array.
real(real64),         intent(in)  :: val(:)   !! Values to search.
integer, allocatable, intent(out) :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.
integer :: i, j, k, n, p

if (csize(dbl2=v, d=1) <= 0 .or. csize(dbl2=v, d=2) <= 0 .or. csize(dbl1=val, d=1) <= 0) return
! allocate row, col
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v, abs(v-val(j)) <= epsilon(1._real64)),1)
end do
call alloc(res, 2, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(v,2)
    do k = 1, size(val,1)
      if (abs(v(i,j)-val(k)) <= epsilon(1._real64)) then
        res(1,p) = i
        res(2,p) = j
        p = p+1
        exit
      end if
    end do
    if (p > n) return
  end do
end do
end subroutine

!-----------------------------------------------------------------------
! find_sca_prv (find)
!-----------------------------------------------------------------------
function find_sca_prv(v, val) result(res)
!! Finds all the occurrences of scalar value in a rank 2 array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: find`  
!! `implicit none`  
!! `integer :: i`  
!! `integer :: v(4,4) = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `integer, allocatable :: pos(:,:)`  
!! `allocate(pos, source=find(v, 4._real64))`  
!! `! the latter is not valid in Gfortan 4.9 and below; write the following code or use sfind instead:`  
!! `! pos = find(v, 4._real64)`  
!! `end program`  
real(real64), intent(in) :: v(:,:)   !! Array where the search is done.
real(real64), intent(in) :: val      !! Value to search.
integer, allocatable     :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.

call sfind(v, val, res)
end function

!-----------------------------------------------------------------------
! find_vec_prv (find)
!-----------------------------------------------------------------------
function find_vec_prv(v, val) result(res)
!! Finds all the occurrences of the elements of a rank 1 array in a rank 2 array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: find`  
!! `implicit none`  
!! `integer :: i`  
!! ``integer :: v(4,4) = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `integer, allocatable :: pos(:,:)`  
!! `allocate(pos, source=find(v, real([4,6],real64)))`  
!! `! the latter is not valid in Gfortan 4.9 and below; write the following code or use sfind instead:`  
!! `! pos = find(v, real([4,6],real64))`  
!! `end program`  
real(real64), intent(in) :: v(:,:)   !! Array where the search is done.
real(real64), intent(in) :: val(:)   !! Values to search.
integer, allocatable     :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.

call sfind(v, val, res)
end function

!-----------------------------------------------------------------------
! find_section_sorted_prv (find_sorted)
!-----------------------------------------------------------------------
function find_section_sorted_prv(dimen, v, val, used) result(pos)
!! Finds the position of a section (row/column) in a (row/column)-sorted array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: alloc, find_sorted`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:)`  
!! `integer:: i`  
!! `call alloc(v, 4, 4)`  
!! `v = reshape([(real(i,real64), i=1,16)], [4,4])`  
!! `print*, find_sorted(1, v, real([3, 7, 11, 15],real64))`  
!! `end program`  
integer,           intent(in) :: dimen  !! Dimension to work with.
real(real64),      intent(in) :: v(:,:) !! Array.
real(real64),      intent(in) :: val(:) !! Section values.
integer, optional, intent(in) :: used   !! Number of sections actually used in vector.
integer                       :: pos    !! Position in the given dimension of the found section.
integer :: n, a, b, anew, bnew, i, j

pos = -1
if (csize(dbl2=v, d=1) <= 0 .or. csize(dbl2=v, d=2) <= 0) return
!number of existing sections
if (present(used)) then; n = used
else;                    n = size(v,dimen)
end if
select case(dimen)
case (1) !Row
  !search among the first vertices
  pos = bsearch(v(1:n,1), val(1), n)
  if (pos <= 0) return
  a=1; b = n
  do j = 2, size(val,1)
    !determine the left extreme of the interval where to search the j-th vertex
    anew = pos
    do i = pos-1, a, -1
      if (v(i,j-1) /= v(pos,j-1)) exit
      anew = i
    end do
    !determine the right extreme of the interval where to search the j-th vertex
    bnew = pos
    do i = pos+1, b
      if (v(i,j-1) /= v(pos,j-1)) exit
      bnew = i
    end do
    a = anew; b = bnew
    pos = bsearch(v(a:b,j), val(j), b-a+1)
    if (pos <= 0) return
    pos = pos + a-1
  end do
case (2) !Column
  !search among the first vertices
  pos = bsearch(v(1,1:n), val(1), n)
  if (pos <= 0) return
  a=1; b = n
  do j = 2, size(val,1)
    !determine the left extreme of the interval where to search the j-th vertex
    anew = pos
    do i = pos-1, a, -1
      if (v(j-1,i) /= v(j-1,pos)) exit
      anew = i
    end do
    !determine the right extreme of the interval where to search the j-th vertex
    bnew = pos
    do i = pos+1, b
      if (v(j-1,i) /= v(j-1,pos)) exit
      bnew = i
    end do
    a = anew; b = bnew
    pos = bsearch(v(j,a:b), val(j), b-a+1)
    if (pos <= 0) return
    pos = pos + a-1
  end do
case default
  call error('(module_alloc_real64_r2/find_section_sorted_prv) 1st argument dimen must be 1 or 2.')
end select
end function
end module

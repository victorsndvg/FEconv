module module_alloc_real64_r3_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for memory allocation of_real64_rank 3 allocatable arrays.  
!!
!! @note Procedures `set`, `add` and `insert` allocate or extend the array when it is necessary.  
!!
! 
! PUBLIC PROCEDURES:
!   dealloc: deallocates memory
!   alloc: allocates memory
!   extend: extends the shape of the array
!   set: sets a scalar, section or matrix value in the array
!   add: adds a scalar, section or matrix value in the array
!   insert: inserts a section in the array
!   reduce: reduces the shape of the array
!   sfind: finds all the occurrences of a scalar value or a vector of values in the array (subroutine version)
!-----------------------------------------------------------------------
use module_compiler_dependant_bmod, only: real64
use module_os_dependant_bmod, only: maxpath
use module_report_bmod, only: error, info
use module_alloc_int_r2_bmod, only: alloc
use module_alloc_common_bmod, only: DEFAULT_ALLOC, csize, search_multiple
implicit none

!Constants
private :: DEFAULT_ALLOC

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_section_prv, set_matrix_prv
private :: add_scalar_prv, add_section_prv, add_matrix_prv
private :: insert_section_prv
private :: sfind_sca_prv, sfind_vec_prv
private :: search_multiple

!Interfaces
interface dealloc
  !! Deallocates a real64 array of rank 3.  
  !!
  !! @note The advantages of using this procedure intead the intrinsic procedure deallocate are:  
  !!  - Testing the previous deallocation of `v` is automatically done.  
  !!  - Deallocation is called with error catching arguments and every error is shown.  
  module procedure dealloc_prv
end interface

interface alloc
  !! Allocates a real64 array of rank 3.  
  !!
  !! @note If `v` is allocated and its extensions are (`d1`,`d2`,`d3`), allocation is not carried out; otherwise, `v` is 
  !! deallocated before allocation. 
  !! If there is not any error, `v` is set to 0.  
  module procedure alloc_prv
end interface

interface extend
  !! Extends the array to contain position (`d1`,`d2`,`d3`).  
  !!
  !! @note Consider the case of `v` being not allocated. For each `i=1,2,3`:  
  !! - If `fit(i)` is present and `.false.`, `v` will be allocated with extension `2^n·DEFAULT_ALLOC` in the dimension `i`, where 
  !! `n` is the smallest exponent such that `2^n·DEFAULT_ALLOC > d(i)` (`DEFAULT_ALLOC` is a private parameter with value `1000`).
  !! - Otherwise, `v` is allocated with extension `d(i)` in the dimension `i`.  
  !!
  !! @note Consider the case of `v` being already allocated. For each `i=1,2,3`, if some `d(i)` is larger than `size(v, i)`, then  
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
  !! Inserts a section in the array.  
  !!
  !! @note The function `extend` is initially called. If `used` is present, `v` is extended in the given dimension until 
  !! `max(used+1, d)`; otherwise, it is extended until `max(size(v,dimen)+1` in the given dimension. Finally, sections after the 
  !! inserted one are shifted one position down and `val` is inserted.  
  !!
  !! @warning When `used` is not present, __the array is  efectively extended in every call__.  
  module procedure insert_section_prv
end interface

interface sfind
  !! Find all the occurrences of a scalar value or an array of values in an array (subroutine version).
  module procedure sfind_sca_prv
  module procedure sfind_vec_prv
end interface

contains

!-----------------------------------------------------------------------
! dealloc_prv (dealloc)
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
!! Deallocates a real64 array of rank 3.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: dealloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call dealloc(v)`  
!! `end program`  
real(real64), allocatable :: v(:,:,:)  !! Array to deallocate.
integer :: res
character(maxpath) :: cad

if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r2/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc_prv (alloc)
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d1, d2, d3)
!! Allocates a real64 array of rank 3.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call alloc(v, 20, 30, 10)`  
!! `end program`  
real(real64), allocatable :: v(:,:,:) !! Array to allocate.
integer, intent(in)       :: d1       !! Array extension in the first dimension.
integer, intent(in)       :: d2       !! Array extension in the second dimension.
integer, intent(in)       :: d3       !! Array extension in the third dimension.
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d1 .and. size(v,2) == d2 .and. size(v,3) == d3) then; v = 0; return; end if
  call dealloc(v)
end if
allocate(v(d1, d2, d3), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_real64_r3/alloc) unable to allocate variable: '//trim(cad))
v = 0
end subroutine

!-----------------------------------------------------------------------
! extend_prv (extend)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d1, d2, d3, fit)
!! Extends the array to contain position (`d1`,`d2`,`d3`).  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, extend`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call alloc(v, 40, 20, 15)`  
!! `call extend(v, 2, 40, 20, fit=[.false.,.true., .true.])`  
!! `end program`
real(real64), allocatable     :: v(:,:,:) !! Array to change dimensions.
integer,           intent(in) :: d1       !! Position in the first dimension to be included in the extension.
integer,           intent(in) :: d2       !! Position in the first dimension to be included in the extension.
integer,           intent(in) :: d3       !! Position in the first dimension to be included in the extension.
logical, optional, intent(in) :: fit(3)   !! Whether or not to fit the array extension to the maximum position in each dimension.
real(real64), allocatable :: temp(:,:,:)
integer :: res, s1, s2, s3, ns1, ns2, ns3
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
    if (fit(3)) then; ns3 = d3                     !we must fit to dimension given as argument
    else; ns3 = search_multiple(DEFAULT_ALLOC, d3) !a multiple of DEFAULT_ALLOC must be taken as new dimension
    endif
  else; ns1 = d1; ns2 = d2; ns3 = d3               !fit is not present, the same as if it where .true.
  end if
  !ALLOCATION
  allocate(v(ns1, ns2, ns3), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_alloc_real64_r3/extend) unable to allocate variable v: '//trim(cad))
  v = 0
else !v is already allocated
  s1 = size(v,1); s2 = size(v,2); s3 = size(v,3)
  if (d1 > s1 .or. d2 > s2 .or. d3 > s3) then !reallocation is mandatory
    !DIMENSIONS
    if (present(fit)) then
      if (fit(1)) then; ns1 = max(s1, d1)          !we must fit to dimension given as argument, if necessary
      else; ns1 = search_multiple(s1, d1)          !a multiple of the current size must be taken as new dimension
      end if
      if (fit(2)) then; ns2 = max(s2, d2)          !we must fit to dimension given as argument, if necessary
      else; ns2 = search_multiple(s2, d2)          !a multiple of the current size must be taken as new dimension
      end if
      if (fit(3)) then; ns3 = max(s3, d3)          !we must fit to dimension given as argument, if necessary
      else; ns3 = search_multiple(s3, d3)          !a multiple of the current size must be taken as new dimension
      end if
    else; ns1 = max(s1, d1); ns2 = max(s2, d2); ns3 = max(s3, d3)!fit is not present, the same as if it where .true.
    end if
    !REALLOCATION
    allocate(temp(ns1, ns2, ns3), stat = res, errmsg = cad)
    if (res /= 0) call error('(module_alloc_real64_r3/extend) unable to allocate variable temp: '//trim(cad))
    temp = 0
    temp(1:s1,1:s2,1:s3) = v
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce_prv (reduce)
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d1, d2, d3)
!! Reduces the extension of an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: alloc, reduce`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call alloc(v, 40, 40, 40)`  
!! `call reduce(v, 20, 20, 20)`  
!! `end program`  
real(real64), allocatable :: v(:,:,:) !! Array to reduce.
integer, intent(in)       :: d1       !! Array extension in the first dimension.
integer, intent(in)       :: d2       !! Array extension in the second dimension.
integer, intent(in)       :: d3       !! Array extension in the third dimension.
real(real64), allocatable :: temp(:,:,:)

if (.not. allocated(v)) then
  call info('(module_alloc_real64_r3/reduce) Variable not allocated'); return
end if
if (size(v,1) == d1 .and. size(v,2) == d2 .and. size(v,3) == d3) return !rows and cols have the right size
if (size(v,1) <  d1 .or.  size(v,2) <  d2 .or. size(v,3) < d3) then     !rows or cols are too large
  call info('(module_alloc_real64_r3/reduce) Some given dimension is too large to reduce'); return
end if
call alloc(temp, d1, d2, d3)
temp(1:d1, 1:d2, 1:d3) = v(1:d1, 1:d2, 1:d3)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set_scalar_prv (set)
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d1, d2, d3, fit)
!! Sets a scalar value in an array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call set(v, 25._real64, 10, 20, 5, fit=[.true.,.false.,.true.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:,:) !! Array
real(real64),      intent(in) :: val      !! Value
integer,           intent(in) :: d1       !! Position in the first dimension
integer,           intent(in) :: d2       !! Position in the second dimension
integer,           intent(in) :: d3       !! Position in the third dimension
logical, optional, intent(in) :: fit(3)   !! Whether or not to fit the array extension to the maximum position in each dimension

call extend(v, d1, d2, d3, fit)
v(d1, d2, d3) = val
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
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call set(1, v, reshape(real([1,2,3,4],real64),[2,2]), 10, fit=[.true.,.false.,.true.])`  
!! `end program`  
integer,           intent(in) :: dimen    !! Dimension to work with.
real(real64), allocatable     :: v(:,:,:) !! Array.
real(real64),      intent(in) :: val(:,:) !! Section values.
integer,           intent(in) :: d        !! Position in the given dimension of the section to insert.
logical, optional, intent(in) :: fit(3)   !! Whether to fit the array extension to the maximum position; by default, `.true.`

select case(dimen)
case (1) ! Row
  call extend(v,           d,   size(val,1),   size(val,2), fit)
  v(                       d, 1:size(val,1), 1:size(val,2)) = val
case (2) ! Column
  call extend(v, size(val,1),             d,   size(val,2), fit)
  v(           1:size(val,1),             d, 1:size(val,2)) = val
case (3) ! Page
  call extend(v, size(val,1),   size(val,2),              d, fit)
  v(           1:size(val,1), 1:size(val,2),              d) = val
case default
  call error('(module_alloc_real64_r2/set_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! set_matrix_prv (set)
!-----------------------------------------------------------------------
subroutine set_matrix_prv(v, val, d1, d2, d3, fit)
!! Sets a matrix in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: set`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call set(v, reshape(real([1,2,3,4,5,6,7,8],real64),[2,2,2]), [10,20], [3,4], [1,3], fit=[.true.,.false.,.true.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:,:)   !! Array.
real(real64),      intent(in) :: val(:,:,:) !! Values.
integer,           intent(in) :: d1(:)      !! Positions in the first dimension.
integer,           intent(in) :: d2(:)      !! Positions in the second dimension.
integer,           intent(in) :: d3(:)      !! Positions in the third dimension.
logical, optional, intent(in) :: fit(3)     !! Whether to fit the array extension to the maximum position in each dimension.

call extend(v, maxval(d1), maxval(d2), maxval(d3), fit)
v(d1, d2, d3) = val
end subroutine

!-----------------------------------------------------------------------
! add_scalar_prv (add)
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d1, d2, d3, fit)
!! Adds a scalar in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call add(v, 25._real64, 10, 20, 7, fit=[.true.,.false.,.true.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:,:) !! Array.
real(real64),      intent(in) :: val      !! Value.
integer,           intent(in) :: d1       !! Position in the first dimension.
integer,           intent(in) :: d2       !! Position in the second dimension.
integer,           intent(in) :: d3       !! Position in the third dimension.
logical, optional, intent(in) :: fit(3)   !! Whether to fit the array extension to the maximum position in each dimension.

call extend(v, d1, d2, d3, fit)
v(d1, d2, d3) = v(d1, d2, d3) + val
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
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call add(1, v, reshape(real([1,2,3,4],real64),[2,2]), 10, fit=[.true.,.false.,.true.])`  
!! `end program`  
integer,           intent(in) :: dimen    !! Dimension to work with.
real(real64), allocatable     :: v(:,:,:) !! Array.
real(real64),      intent(in) :: val(:,:) !! Section values.
integer,           intent(in) :: d        !! Position in the given dimension of the section to insert.
logical, optional, intent(in) :: fit(3)   !! Whether to fit the array extension to the maximum position; by default, `.true.`

select case(dimen)
case (1) ! Row
  call extend(v,           d,   size(val,1),   size(val,2), fit)
  v(                       d, 1:size(val,1), 1:size(val,2)) = v(            d, 1:size(val,1), 1:size(val,2)) + val
case (2) ! Column
  call extend(v, size(val,1),             d,   size(val,2), fit)
  v(           1:size(val,1),             d, 1:size(val,2)) = v(1:size(val,1),             d, 1:size(val,2)) + val
case (3) ! Page
  call extend(v, size(val,1),   size(val,2),             d, fit)
  v(           1:size(val,1), 1:size(val,2),             d) = v(1:size(val,1), 1:size(val,2),             d) + val
case default
  call error('(module_alloc_real64_r2/add_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! add_matrix_prv (add)
!-----------------------------------------------------------------------
subroutine add_matrix_prv(v, val, d1, d2, d3, fit)
!! Adds a matrix in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: add`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call add(v, reshape(real([1,2,3,4,5,6,7,8],real64),[2,2,2]), [10,20], [3,4], [1,3], fit=[.true.,.false.,.true.])`  
!! `end program`  
real(real64), allocatable     :: v(:,:,:)   !! Array.
real(real64),      intent(in) :: val(:,:,:) !! Values.
integer,           intent(in) :: d1(:)      !! Positions in the first dimension.
integer,           intent(in) :: d2(:)      !! Positions in the second dimension.
integer,           intent(in) :: d3(:)      !! Positions in the third dimension.
logical, optional, intent(in) :: fit(3)     !! Whether to fit the array extension to the maximum position in each dimension.

call extend(v, maxval(d1), maxval(d2), maxval(d3), fit)
v(d1, d2, d3) = v(d1, d2, d3) + val
end subroutine

!-----------------------------------------------------------------------
! insert_section_prv (insert)
!-----------------------------------------------------------------------
subroutine insert_section_prv(dimen, v, val, d, used, fit)
!! Inserts a section in the array.  
!! __Example:__  
!! `program test`  
!! `use basicmod, only: insert`  
!! `implicit none`  
!! `real(real64), allocatable :: v(:,:,:)`  
!! `call insert(1, v, reshape(real([1,2,3,4],real64),[2,2]), 10, used=3, fit=[.true.,.false.,.true.])`  
!! `end program`  
integer,           intent(in) :: dimen    !! Dimension to work with.
real(real64), allocatable     :: v(:,:,:) !! Array.
real(real64),      intent(in) :: val(:,:) !! Section values.
integer,           intent(in) :: d        !! Position in the given dimension of the section to insert.
integer, optional, intent(in) :: used     !! Number of sections actually used in vector.
logical, optional, intent(in) :: fit(3)   !! Whether to fit the array extension to the maximum position; by default, `.true.`
integer :: s

if (1 <= dimen .and. dimen <= 3) then
  if (present(used)) then; s = max(used+1, d)
  else; s = max(size(v,dimen)+1, d)
  end if
end if
select case(dimen)
case (1) ! Row
  call extend(v, s, size(val,1), size(val,2), fit)
  v(d+1:size(v,dimen), :,                 :)                 = v(d:size(v,dimen)-1, :, :)
  v(d,                 :,                 :)                 = 0
  v(d,                 1:size(val,1),     1:size(val,2))     = val
case (2) ! Column
  call extend(v, size(val,1), s, size(val,2), fit)
  v(:,                 d+1:size(v,dimen), :)                 = v(:, d:size(v,dimen)-1, :)
  v(:,                 d,                 :)                 = 0
  v(1:size(val,1),     d,                 1:size(val,2))     = val
case (3) ! Page
  call extend(v, size(val,1), size(val,2), s, fit)
  v(:,                 :,                 d+1:size(v,dimen)) = v(:, :, d:size(v,dimen)-1)
  v(:,                 :,                 d)                 = 0
  v(1:size(val,1),     1:size(val,2),     d)                 = val
case default
  call error('(module_alloc_real64_r2/insert_section_prv) 1st argument dimen must be 1 or 2.')
end select
end subroutine

!-----------------------------------------------------------------------
! sfind_sca_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, res)
!! Finds all the occurrences of scalar value in a rank 3 array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: sfind`  
!! `implicit none`  
!! `integer :: i`  
!! `real(real64) :: v(4,4,4) = reshape([(real(i,real64), i=1,64)], [4,4,4])`  
!! `integer, allocatable :: res(:,:)`  
!! `call sfind(v, 2._real64, res)`  
!! `end program`  
real(real64),         intent(in)  :: v(:,:,:) !! Array.
real(real64),         intent(in)  :: val      !! Value to search.
integer, allocatable, intent(out) :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.
integer :: i, j, k, n, p

if (csize(dbl3=v, d=1) <= 0 .or. csize(dbl3=v, d=2) <= 0 .or. csize(dbl3=v, d=3) <= 0) return
! allocate row, col
n = size(pack(v,v==val),1)
call alloc(res, 3, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(v,2)
    do k = 1, size(v,3)
      if (v(i,j,k) == val) then
        res(1,p) = i
        res(2,p) = j
        res(3,p) = k
        p = p+1
      end if
      if (p > n) return
    enddo
  end do
end do
end subroutine

!-----------------------------------------------------------------------
! sfind_vec_prv (sfind)
!-----------------------------------------------------------------------
subroutine sfind_vec_prv(v, val, res)
!! Finds all the occurrences of the elements of a rank 1 array in a rank 3 array.  
!! __Example__  
!! `program test`  
!! `use basicmod, only: alloc, sfind`  
!! `implicit none`  
!! `integer :: i`  
!! `real(real64) :: v(4,4,4) = reshape([(real(i,real64), i=1,64)], [4,4,4]), val(8) = [(2._real64*i, i=1,8)]`  
!! `integer, allocatable :: res(:,:)`  
!! `call sfind(v, val, res)`  
!! `end program`  
real(real64),         intent(in)  :: v(:,:,:) !! Array.
real(real64),         intent(in)  :: val(:)   !! Values to search.
integer, allocatable, intent(out) :: res(:,:) !! Array of positions; a zero-sized array if nothing was found.
integer :: i, j, k, l, n, p

if (csize(dbl3=v, d=1) <= 0 .or. csize(dbl3=v, d=2) <= 0 .or. csize(dbl3=v, d=3) <= 0) return
! allocate row, col
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v,v==val(j)),1)
end do
call alloc(res, 3, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(v,2)
    do k = 1, size(v,3)
      do l = 1, size(val,1)
        if (v(i,j,k) == val(l)) then
          res(1,p) = i
          res(2,p) = j
          res(3,p) = k
          p = p+1
          exit
        end if
      end do
      if (p > n) return
    enddo
  end do
end do
end subroutine
end module

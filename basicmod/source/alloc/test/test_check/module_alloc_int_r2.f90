module module_alloc_int_r2
!-----------------------------------------------------------------------
! Module for memory allocation of integer rank 2 allocatable arrays in
! check mode: reports any problem
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 01/03/2012
!
! PUBLIC PROCEDURES:
!   dealloc: dealloc memory
!   alloc: alloc memory
!   extend: extend the extension of an array
!   set: set a scalar or a matrix in the array
!   set_row: set a row in the array
!   set_col: set a column in the array
!   add: add a scalar or a matrix in the array
!   add_row: add a row in the array
!   add_col: add a column in the array
!   insert_row: insert a row in the array
!   insert_col: insert a column in the array
!   reduce: reduce the array
!-----------------------------------------------------------------------
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_convers, only: string
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_row_prv, set_col_prv, set_matrix_prv
private :: add_scalar_prv, add_row_prv, add_col_prv, add_matrix_prv
private :: insert_row_prv, insert_col_prv
private :: search_multiple

!Interface
interface    dealloc; module procedure    dealloc_prv; end interface
interface      alloc; module procedure      alloc_prv; end interface
interface     extend; module procedure     extend_prv; end interface
interface     reduce; module procedure     reduce_prv; end interface
interface        set; module procedure set_scalar_prv; end interface
interface        set; module procedure set_matrix_prv; end interface
interface    set_row; module procedure    set_row_prv; end interface
interface    set_col; module procedure    set_col_prv; end interface
interface        add; module procedure add_scalar_prv; end interface
interface        add; module procedure add_matrix_prv; end interface
interface    add_row; module procedure    add_row_prv; end interface
interface    add_col; module procedure    add_col_prv; end interface
interface insert_row; module procedure insert_row_prv; end interface
interface insert_col; module procedure insert_col_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
integer, allocatable :: v(:,:)
integer :: res
character(maxpath) :: cad

if (.not. allocated(v)) then
  call info('(module_alloc_int_r2/dealloc) v is already not allocated')
  return
end if  
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_int_r2/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d1, d2)
integer, allocatable :: v(:,:)
integer, intent(in)  :: d1, d2
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  call info('(module_alloc_int_r2/alloc) v is already allocated')
  if (size(v,1) == d1 .and. size(v,2) == d2) then; v = 0; return; end if
  call dealloc(v)
end if
allocate(v(d1, d2), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_int_r2/alloc) unable to allocate variable: '//trim(cad))
v = 0
end subroutine

!-----------------------------------------------------------------------
! extend: extend the array to contain position (d1,d2)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d1, d2, fit)
integer, allocatable          :: v(:,:), temp(:,:)
integer, intent(in)           :: d1, d2 !new dimensions given by the user
logical, intent(in), optional :: fit(2) 
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
  if (res /= 0) call error('(module_alloc_int_r2/extend) unable to allocate variable v: '//trim(cad))
  v = 0
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
    if (res /= 0) call error('(module_alloc_int_r2/extend) unable to allocate variable temp: '//trim(cad))
    temp = 0
    temp(1:s1,1:s2) = v
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d1, d2)
integer, allocatable :: v(:,:), temp(:,:)
integer, intent(in)  :: d1, d2

if (.not. allocated(v)) then 
  call info('(module_alloc_int_r2/reduce) Variable not allocated'); return
end if
if (size(v,1) == d1 .and. size(v,2) == d2) return !rows and cols have the right size
if (size(v,1) <  d1 .or.  size(v,2) <  d2) then   !rows or cols are too large
  call info('(module_alloc_int_r2/reduce) Some given dimension is too large to reduce'); return
end if
call alloc(temp, d1, d2)
temp(1:d1, 1:d2) = v(1:d1, 1:d2)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set: set a scalar in the array
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d1, d2, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val, d1, d2
logical, intent(in), optional :: fit(2)

if (d1 < 1 .or. d1 > size(v,1) .or. d2 < 1 .or. d2 > size(v,2)) call info('(module_alloc_int_r2/set_scalar) Given dimensions: '&
//trim(string(d1))//', '//trim(string(d2))//', size(v): '//trim(string(size(v,1)))//', '//trim(string(size(v,2))))
call extend(v, d1, d2, fit)
v(d1, d2) = val
end subroutine

!-----------------------------------------------------------------------
! set_row: set a row in the array
!-----------------------------------------------------------------------
subroutine set_row_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,1)) call info('(module_alloc_int_r2/set_row) Given dimension: '&
//trim(string(d))//', size(v,1): '//trim(string(size(v,1))))
call extend(v, d, size(val,1), fit)
v(d, 1:size(val,1)) = val
end subroutine

!-----------------------------------------------------------------------
! set_col: set a column in the array
!-----------------------------------------------------------------------
subroutine set_col_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,2)) call info('(module_alloc_int_r2/set_col) Given dimension: '&
//trim(string(d))//', size(v,2): '//trim(string(size(v,2))))
call extend(v, size(val,1), d, fit)
v(1:size(val,1), d) = val
end subroutine

!-----------------------------------------------------------------------
! set: set a matrix in the array
!-----------------------------------------------------------------------
subroutine set_matrix_prv(v, val, d1, d2, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d1(:), d2(:)
logical, intent(in), optional :: fit(2)

if (minval(d1) < 1 .or. maxval(d1) > size(v,1)) call info('(module_alloc_int_r2/set_matrix) Given dimensions(1): ['&
//trim(string(minval(d1)))//', '//trim(string(maxval(d1)))//'], size(v,1): '//trim(string(size(v,1))))
if (minval(d2) < 1 .or. maxval(d2) > size(v,2)) call info('(module_alloc_int_r2/set_matrix) Given dimensions(2): ['&
//trim(string(minval(d2)))//', '//trim(string(maxval(d2)))//'], size(v,2): '//trim(string(size(v,2))))
if (size(d1,1)*size(d2,1) /= size(val,1)) call info('(module_alloc_int_r2/set_matrix) Given index dimensions: '&
//trim(string(size(d1,1)))//'x'//trim(string(size(d2,1)))//', size(val): '//trim(string(size(val,1))))
call extend(v, maxval(d1), maxval(d2), fit)
v(d1, d2) = reshape(val,[size(d1), size(d2)])
end subroutine

!-----------------------------------------------------------------------
! add: add a scalar in the array
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d1, d2, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val, d1, d2
logical, intent(in), optional :: fit(2)

if (d1 < 1 .or. d1 > size(v,1) .or. d2 < 1 .or. d2 > size(v,2)) call info('(module_alloc_int_r2/add_scalar) Given dimensions: '&
//trim(string(d1))//', '//trim(string(d2))//', size(v): '//trim(string(size(v,1)))//', '//trim(string(size(v,2))))
call extend(v, d1, d2, fit)
v(d1, d2) = v(d1, d2) + val
end subroutine

!-----------------------------------------------------------------------
! add_row: add a row in the array
!-----------------------------------------------------------------------
subroutine add_row_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,1)) call info('(module_alloc_int_r2/add_row) Given dimension: '&
//trim(string(d))//', size(v,1): '//trim(string(size(v,1))))
call extend(v, d, size(val,1), fit)
v(d, 1:size(val,1)) = v(d, 1:size(val,1)) + val
end subroutine

!-----------------------------------------------------------------------
! add_col: add a column in the array
!-----------------------------------------------------------------------
subroutine add_col_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,2)) call info('(module_alloc_int_r2/add_col) Given dimension: '&
//trim(string(d))//', size(v,2): '//trim(string(size(v,2))))
call extend(v, size(val,1), d, fit)
v(1:size(val,1), d) = v(1:size(val,1), d) + val
end subroutine

!-----------------------------------------------------------------------
! add: add a matrix in the array
!-----------------------------------------------------------------------
subroutine add_matrix_prv(v, val, d1, d2, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d1(:), d2(:)
logical, intent(in), optional :: fit(2)

if (minval(d1) < 1 .or. maxval(d1) > size(v,1)) call info('(module_alloc_int_r2/add_matrix) Given dimensions(1): ['&
//trim(string(minval(d1)))//', '//trim(string(maxval(d1)))//'], size(v,1): '//trim(string(size(v,1))))
if (minval(d2) < 1 .or. maxval(d2) > size(v,2)) call info('(module_alloc_int_r2/add_matrix) Given dimensions(2): ['&
//trim(string(minval(d2)))//', '//trim(string(maxval(d2)))//'], size(v,2): '//trim(string(size(v,2))))
if (size(d1,1)*size(d2,1) /= size(val,1)) call info('(module_alloc_int_r2/add_matrix) Given index dimensions: '&
//trim(string(size(d1,1)))//'x'//trim(string(size(d2,1)))//', size(val): '//trim(string(size(val,1))))
call extend(v, maxval(d1), maxval(d2), fit)
v(d1, d2) = v(d1, d2) + reshape(val,[size(d1), size(d2)])
end subroutine

!-----------------------------------------------------------------------
! insert_row: insert a row in the array
!-----------------------------------------------------------------------
subroutine insert_row_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,1)) call info('(module_alloc_int_r2/insert_row) Given dimension: '&
//trim(string(d))//', size(v,1): '//trim(string(size(v,1))))
call extend(v, max(size(v,1)+1, d), size(val,1), fit)
v(d+1:size(v,1), :) = v(d:size(v,1)-1, :)
v(d, 1:size(val,1)) = val
v(d, size(val,1)+1:size(v,2)) = 0
end subroutine

!-----------------------------------------------------------------------
! insert_col: insert a col in the array
!-----------------------------------------------------------------------
subroutine insert_col_prv(v, val, d, fit)
integer, allocatable :: v(:,:)
integer, intent(in)  :: val(:)
integer, intent(in)  :: d
logical, intent(in), optional :: fit(2)

if (d < 1 .or. d > size(v,2)) call info('(module_alloc_int_r2/insert_col) Given dimension: '&
//trim(string(d))//', size(v,2): '//trim(string(size(v,2))))
call extend(v, size(val,1), max(size(v,2)+1, d), fit)
v(:, d+1:size(v,2)) = v(:, d:size(v,2)-1)
v(1:size(val,1), d) = val
v(size(val,1)+1:size(v,1), d) = 0
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

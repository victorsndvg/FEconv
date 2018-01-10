module module_ALLOC_log_r2_fcnv
!-----------------------------------------------------------------------
! Module for memory allocation of logical rank 2 allocatable arrays
! Last update: 28/07/2009
! Programmer: fran.pena@usc.es
!
! PUBLIC PROCEDURES: dealloc, alloc, extend, reduce, enlarge and set
! - you only have to use "set" since it implies allocation, if necessary
! - "reduce" cuts off an array filled without fitting
! - "enlarge" ensures the allocation of a given position
!-----------------------------------------------------------------------
use basicmod, only: error, info
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation
integer, parameter, private :: DEFAULT_EXTEND = 1000 !additional size for extension

!Private procedures
private :: dealloc_prv,  alloc_prv,  extend_prv,  reduce_prv, enlarge_prv
private :: set_prv, setv_prv, enlarge1

!Interface
interface dealloc; module procedure dealloc_prv; end interface
interface   alloc; module procedure   alloc_prv; end interface
interface  extend; module procedure  extend_prv; end interface
interface  reduce; module procedure  reduce_prv; end interface
interface enlarge; module procedure enlarge_prv; end interface
interface     set; module procedure     set_prv; end interface
interface     set; module procedure    setv_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)

  logical, dimension(:,:), allocatable :: v
  integer :: res

  if (.not. allocated(v)) return
  deallocate(v, stat = res) !In f2003 add: errmsg = cad
  if (res /= 0) call error('(module_ALLOC_log_r2/dealloc) Unable to deallocate variable')

end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory
!-----------------------------------------------------------------------
subroutine alloc_prv(v, rows, cols)

  logical, dimension(:,:), allocatable :: v
  integer, intent(in), optional     :: rows, cols
  integer :: res, n, m

  n = DEFAULT_ALLOC; if (present(rows)) then
    if (rows <= 0) return; n = rows
  end if
  m = DEFAULT_ALLOC; if (present(cols)) then
    if (cols <= 0) return; m = cols
  end if
  if (allocated(v)) then
    if (size(v,1) == n .and. size(v,2) == m) then
      v(1:n, 1:m) = .false.; return
    end if
    call dealloc(v)
  end if
  allocate(v(n, m), stat = res) !In f2003 add: errmsg = cad
  if (res /= 0) call error('(module_alloc_log_r2/alloc) unable to allocate variable')
  v(1:n, 1:m) = .false.

end subroutine

!-----------------------------------------------------------------------
! extend: extend the array in the dimension dim
!-----------------------------------------------------------------------
subroutine extend_prv(v, dim, extension)

  logical, dimension(:,:), allocatable :: v, temp
  integer, intent(in)               :: dim
  integer, intent(in), optional     :: extension
  integer :: ext, n0, m0, n, m

  ext = DEFAULT_EXTEND !extension
  if (present(extension))then
    if (extension <= 0) return
    ext = extension
  end if
  if (.not. allocated(v)) then !v not allocated
     call alloc(v, ext, ext); return
  end if
  if (dim == 1) then; n = size(v,1) + ext; m = size(v,2) !selection of dimensions
  else;               n = size(v,1);       m = size(v,2) + ext
  end if
  n0 = size(v,1); m0 = size(v,2)
  call alloc(temp, n, m)        !temporal storage
  temp(1:n0, 1:m0) = v
  temp(n0+1:n, 1:m0) = .false.  !reset new positions
  temp(1:n0, m0+1:m) = .false.  !reset new positions
  call dealloc(v)
  call alloc(v, n, m)           !final allocation
  v(1:n, 1:m) = temp            !data copy
  call dealloc(temp)

end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, rows, cols)

  logical, dimension(:,:), allocatable :: v, temp
  integer, intent(in)               :: rows, cols

  if (.not. allocated(v)) then !v not allocated
    call info('(module_ALLOC_log_r2/reduce) Variable not allocated'); return
  end if
  if (rows == size(v,1) .and. cols == size(v,2)) return !rows and cols have the right size
  if (rows > size(v,1) .or. cols > size(v,2)) then !rows or cols too large
    call info('(module_ALLOC_log_r2/reduce) Given dimension is too large'); return
  end if
  call alloc(temp, rows, cols)             !temporal storage
  temp(1:rows, 1:cols) = v(1:rows, 1:cols) !data copy
  call dealloc(v)
  call alloc(v, rows, cols)                !final allocation
  v(1:rows, 1:cols) = temp                 !data copy
  call dealloc(temp)

end subroutine

!-----------------------------------------------------------------------
! enlarge: enlarge an array for ensuring a given position
!-----------------------------------------------------------------------
subroutine enlarge_prv(v, row, col, fit_row, fit_col)

  logical, dimension(:,:), allocatable :: v
  integer, intent(in), optional     :: row, col
  logical, intent(in), optional     :: fit_row, fit_col
  integer :: r, c

  r = enlarge1(v, 1, row, fit_row) !extension in rows, if necessary
  if (present(fit_col) .and. present(col)) then; if (fit_col) call reduce(v, size(v,1), col); end if
  c = enlarge1(v, 2, col, fit_col) !extension in cols, if necessary

end subroutine

!-----------------------------------------------------------------------
! set: set a value to the array
!-----------------------------------------------------------------------
subroutine set_prv(v, val, row, col, fit_row, fit_col, add)

  logical, dimension(:,:), allocatable :: v
  logical, intent(in)                  :: val
  integer, intent(in), optional     :: row, col
  logical, intent(in), optional     :: fit_row, fit_col, add
  integer :: r, c
  logical :: preval

  r = enlarge1(v, 1, row, fit_row) !extension in rows, if necessary
  c = enlarge1(v, 2, col, fit_col) !extension in cols, if necessary
  preval = .false. !previous value
  if (present(add)) then; if (add) preval = v(r, c); end if
  v(r, c) = preval .or. val

end subroutine

!-----------------------------------------------------------------------
! set: set a vector to the array
!-----------------------------------------------------------------------
subroutine setv_prv(v, val, row, col, fit_row, fit_col, add)

  logical, dimension(:,:), allocatable :: v
  logical, dimension(:), intent(in)    :: val
  integer, intent(in), optional     :: row, col
  logical, intent(in), optional     :: fit_row, fit_col, add
  logical, dimension(size(val,1)) :: preval
  integer :: r, c

  if (present(row) .eqv. present(col)) call error('(module_ALLOC_log_r2/set) &
  &Either row or col (not both) must appear when setting a vector')
  if (present(row)) then
    r = enlarge1(v, 1, row, fit_row) !extension in rows, if necessary
    c = enlarge1(v, 2, size(val,1), fit_col) !extension in cols, if necessary
    preval = .false. !previous value
    if (present(add)) then; if (add) preval = v(row, 1:size(val,1)); end if
    v(row, 1:size(val,1)) = preval .or. val
  elseif (present(col)) then
    r = enlarge1(v, 1, size(val,1), fit_row) !extension in rows, if necessary
    c = enlarge1(v, 2, col, fit_col) !extension in cols, if necessary
    preval = .false. !previous value
    if (present(add)) then; if (add) preval = v(1:size(val,1), col); end if
    v(1:size(val,1), col) = preval .or. val
  end if

end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! enlarge1: enlarge an array
! RESULT: p, position prepared to set data
!-----------------------------------------------------------------------
function enlarge1(v, dim, puser, fit) result(p)

  logical, dimension(:,:), allocatable :: v
  integer, intent(in)               :: dim
  integer, intent(in), optional     :: puser
  logical, intent(in), optional     :: fit
  integer :: p, n

  n = 0 !size in dimension
  if (allocated(v)) n = size(v, dim)
  p = n + 1 !position
  if (present(puser)) p = puser
  !ALLOCATE
  if (present(fit)) then
    if (fit) then; call extend(v, dim, extension=p-n)   !fit is true
    else; do while (n < p); call extend(v, dim); n = size(v,dim); end do !fit is false
    end if
  else; do while (n < p); call extend(v, dim); n = size(v,dim); end do   !fit is false
  end if

end function

end module

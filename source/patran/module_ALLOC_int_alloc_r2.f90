module module_ALLOC_int_alloc_r2_fcnv
!-----------------------------------------------------------------------
! Module for memory allocation of allocatable arrays composed of
! integer allocatable arrays
! Last update: 28/07/2009
! Programmer: fran.pena@usc.es
! 
! PUBLIC PROCEDURES: dealloc, alloc, extend, reduce, enlarge and set
! - you only have to use "set" since it implies allocation, if necessary
! - "reduce" cuts off an array filled without fitting
! - "enlarge" ensures the allocation of a given position
!-----------------------------------------------------------------------
use basicmod
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 10 !initial size for allocation
integer, parameter, private :: DEFAULT_EXTEND = 10 !additional size for extension

!Types
type int_alloc_r1 !every row is a rank 1 allocatable array
  integer, dimension(:), allocatable :: col
end type

type int_alloc_r2 !derived data type is composed of rows
  type(int_alloc_r1), dimension(:), allocatable :: row
end type

!Private procedures
private :: dealloc_prv,  alloc_prv,  extend_prv,  reduce_prv, enlarge_prv, set_prv
private :: enlarge_rows, enlarge_cols, transf

!Interface
interface dealloc; module procedure dealloc_prv; end interface
interface   alloc; module procedure   alloc_prv; end interface
interface  extend; module procedure  extend_prv; end interface
interface  reduce; module procedure  reduce_prv; end interface
interface enlarge; module procedure enlarge_prv; end interface
interface     set; module procedure     set_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc rows and cols
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)

  type(int_alloc_r2) :: v
  integer :: res, i

  if (.not. allocated(v%row)) return
  do i = 1, size(v%row, 1)
    if (.not. allocated(v%row(i)%col)) cycle
    call dealloc(v%row(i)%col)
  end do
  deallocate(v%row, stat = res) !In f2003 add: errmsg = cad
  if (res /= 0) call error('(module_ALLOC_int_alloc_r2/dealloc) Unable to deallocate variable')

end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory for rows, but dealloc cols
!-----------------------------------------------------------------------
subroutine alloc_prv(v, rows)

  type(int_alloc_r2) :: v
  integer, optional, intent(in) :: rows
  integer :: res, n, i

  n = DEFAULT_ALLOC; if (present(rows)) then
    if (rows <= 0) return; n = rows
  end if
  if (allocated(v%row)) call dealloc(v)
  allocate(v%row(n), stat = res) !In f2003 add: errmsg = cad
  if (res /= 0) call error('(module_ALLOC_int_alloc_r2/alloc) Unable to allocate variable')
  do i = 1, n
    if (allocated(v%row(i)%col)) call dealloc(v%row(i)%col)
  end do

end subroutine

!-----------------------------------------------------------------------
! extend: extend the array
!-----------------------------------------------------------------------
subroutine extend_prv(v, extension)

  type(int_alloc_r2) :: v
  integer, intent(in), optional :: extension
  type(int_alloc_r2) :: temp
  integer :: ext, m

  ext = DEFAULT_EXTEND !extension
  if (present(extension))then
    if (extension <= 0) return
    ext = extension
  end if
  if (.not. allocated(v%row)) then
    call alloc(v, ext); return
  end if
  m = size(v%row, 1)
  call alloc(temp, m)    !temporal storage
  call transf(v, temp)   !data copy
  call dealloc(v)
  call alloc(v, m + ext) !final allocation
  call transf(temp, v)   !data copy
  call dealloc(temp)

end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, rows)

  type(int_alloc_r2) :: v, temp
  integer, intent(in) :: rows

  if (.not. allocated(v%row)) then
    call info('(module_ALLOC_int_alloc_r2/reduce) Variable not allocated'); return
  end if
  if (rows > size(v%row, 1)) then
    call info('(module_ALLOC_int_alloc_r2/reduce) Given dimension is too large'); return
  end if
  call alloc(temp, rows)  !temporal storage
  call transf(v, temp)    !data copy
  call dealloc(v)
  call alloc(v, rows)     !final allocation
  call transf(temp, v)    !data copy
  call dealloc(temp)

end subroutine

!-----------------------------------------------------------------------
! enlarge: enlarge an array for ensuring a given position
!-----------------------------------------------------------------------
subroutine enlarge_prv(v, row, col, fit_row, fit_col)

  type(int_alloc_r2)           :: v
  integer, intent(in), optional :: row, col
  logical, intent(in), optional :: fit_row, fit_col
  integer :: r, c

  r = enlarge_rows(v, row, fit_row)    !extension in rows, if necessary
  c = enlarge_cols(v, r, col, fit_col) !extension in cols, if necessary

end subroutine

!-----------------------------------------------------------------------
! set: set a value in the array
!-----------------------------------------------------------------------
subroutine set_prv(v, val, row, col, fit_row, fit_col, add)

  type(int_alloc_r2)           :: v
  integer, intent(in)              :: val
  integer, intent(in), optional :: row, col
  logical, intent(in), optional :: fit_row, fit_col, add
  integer :: r, c
  integer :: preval

  r = enlarge_rows(v, row, fit_row)    !extension in rows, if necessary
  c = enlarge_cols(v, r, col, fit_col) !extension in cols, if necessary
  preval = 0 !previous value
  if (present(add)) then; if (add) preval = v%row(r)%col(c); end if
  v%row(r)%col(c) = preval + val

end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! enlarge_rows: enlarge the rows
! RESULT: p, row prepared to set data
!-----------------------------------------------------------------------
function enlarge_rows(v, row, fit) result(p)

  type(int_alloc_r2)           :: v
  integer, intent(in), optional :: row
  logical, intent(in), optional :: fit
  integer :: p, n

  n = 0 !size in dimension
  if (allocated(v%row)) n = size(v%row, 1)
  p = n + 1 !position
  if (present(row)) p = row
  !ALLOCATE
  if (present(fit)) then
    if (fit) then; call extend(v, p-n)  !fit is true
    else !fit is false
      do while (n < p)
        call extend(v)
        n = size(v%row, 1)
      end do
    end if
  else !fit is false
    do while (n < p)
      call extend(v)
      n = size(v%row, 1)
    end do
  end if

end function

!-----------------------------------------------------------------------
! enlarge_cols: enlarge the cols of the row r
! RESULT: p, row prepared to set data
!-----------------------------------------------------------------------
function enlarge_cols(v, r, col, fit) result(p)

  type(int_alloc_r2)            :: v
  integer, intent(in)           :: r
  integer, intent(in), optional :: col
  logical, intent(in), optional :: fit
  integer :: p, n

  n = 0 !size in dimension
  if (allocated(v%row(r)%col)) n = size(v%row(r)%col, 1)
  p = n + 1 !position
  if (present(col)) p = col
  !ALLOCATE
  if (present(fit)) then
    if (fit) then; call extend(v%row(r)%col, p)   !fit is true
    else !fit is false
      do while (n < p)
        call extend(v%row(r)%col, p, fit=.false.)
        n = size(v%row(r)%col, 1)
      end do
    end if
  else !fit is false
    do while (n < p)
      call extend(v%row(r)%col, p, fit=.false.)
      n = size(v%row(r)%col, 1)
    end do
  end if

end function

!-----------------------------------------------------------------------
! transf: transfer data between to int_alloc_r2 objects
!-----------------------------------------------------------------------
subroutine transf(u, v)

  type(int_alloc_r2) :: u, v
  integer :: i

  do i = 1, size(u%row, 1)
    if (allocated(u%row(i)%col)) then
      call alloc(v%row(i)%col, size(u%row(i)%col, 1))
      v%row(i)%col(1:size(u%row(i)%col, 1)) = u%row(i)%col
    end if
  end do

end subroutine

end module

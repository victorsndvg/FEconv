module module_alloc_real64_r3
!-----------------------------------------------------------------------
! Module for memory allocation of integer rank 2 allocatable arrays
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
! Last update: 02/10/2012
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
!   insert_row_sorted: insert a row in a row-sorted array
!   insert_col_sorted: insert a col in a col-sorted array
!   find_row_sorted: find the position of a row in a row-sorted array
!   find_col_sorted: find the position of a col in a col-sorted array
!   sfind: subroutine to find all the occurrences of a value in an array
!     Private functions under 'sfind' are:
!     sfind_sca: find a scalar value in an array
!     sfind_vec: find an array of values in an array
!
! REMARK:
!   find (as function) is not implemented since functions in Fortran cannot
!   return several output entities
!-----------------------------------------------------------------------
use module_compiler_dependant, only: real64
use module_os_dependant, only: maxpath
use module_report, only: error, info
use module_alloc_real64_r1, only: alloc, bsearch 
use module_alloc_real64_r2, only: alloc, bsearch 
implicit none

!Constants
integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation

!Private procedures
private :: dealloc_prv, alloc_prv, extend_prv, reduce_prv
private :: set_scalar_prv, set_row_prv, set_col_prv, set_matrix_prv
private :: add_scalar_prv, add_row_prv, add_col_prv, add_matrix_prv
private :: insert_row_prv, insert_col_prv
!private :: insert_col_sorted_prv, insert_row_sorted_prv
!private :: find_row_sorted_prv, find_col_sorted_prv
private :: sfind_sca_prv, sfind_vec_prv, search_multiple

!Interface
interface           dealloc; module procedure           dealloc_prv; end interface
interface             alloc; module procedure             alloc_prv; end interface
interface            extend; module procedure            extend_prv; end interface
interface            reduce; module procedure            reduce_prv; end interface
interface               set; module procedure        set_scalar_prv; end interface
interface               set; module procedure        set_matrix_prv; end interface
interface           set_row; module procedure           set_row_prv; end interface
interface           set_col; module procedure           set_col_prv; end interface
interface               add; module procedure        add_scalar_prv; end interface
interface               add; module procedure        add_matrix_prv; end interface
interface           add_row; module procedure           add_row_prv; end interface
interface           add_col; module procedure           add_col_prv; end interface
interface        insert_row; module procedure        insert_row_prv; end interface
interface        insert_col; module procedure        insert_col_prv; end interface
!interface insert_row_sorted; module procedure insert_row_sorted_prv; end interface
!interface insert_col_sorted; module procedure insert_col_sorted_prv; end interface
!interface   find_row_sorted; module procedure   find_row_sorted_prv; end interface
!interface   find_col_sorted; module procedure   find_col_sorted_prv; end interface
interface             sfind; module procedure         sfind_sca_prv; end interface
interface             sfind; module procedure         sfind_vec_prv; end interface

contains

!-----------------------------------------------------------------------
! dealloc: dealloc memory
!-----------------------------------------------------------------------
subroutine dealloc_prv(v)
real(real64), allocatable :: v(:,:,:)
integer :: res
character(maxpath) :: cad

if (.not. allocated(v)) return
deallocate(v, stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_int_r2/dealloc) Unable to deallocate variable: '//trim(cad))
end subroutine

!-----------------------------------------------------------------------
! alloc: alloc memory
!-----------------------------------------------------------------------
subroutine alloc_prv(v, d1, d2, d3)
real(real64), allocatable :: v(:,:,:)
integer,      intent(in)  :: d1, d2, d3
integer :: res
character(maxpath) :: cad

if (allocated(v)) then
  if (size(v,1) == d1 .and. size(v,2) == d2 .and. size(v,3) == d3) then; v = 0; return; end if
  call dealloc(v)
end if
allocate(v(d1, d2, d3), stat = res, errmsg = cad)
if (res /= 0) call error('(module_alloc_int_r3/alloc) unable to allocate variable: '//trim(cad))
v = 0
end subroutine

!-----------------------------------------------------------------------
! extend: extend the array to contain position (d1,d2,d3)
!-----------------------------------------------------------------------
subroutine extend_prv(v, d1, d2, d3, fit)
real(real64), allocatable     :: v(:,:,:), temp(:,:,:)
integer, intent(in)           :: d1, d2, d3 !new dimensions given by the user
logical, intent(in), optional :: fit(3) 
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
    if (fit(3)) then; ns3 = d2                     !we must fit to dimension given as argument
    else; ns3 = search_multiple(DEFAULT_ALLOC, d3) !a multiple of DEFAULT_ALLOC must be taken as new dimension
    endif
  else; ns1 = d1; ns2 = d2; ns3 = d3               !fit is not present, the same as if it where .true.
  end if
  !ALLOCATION
  allocate(v(ns1, ns2, ns3), stat = res, errmsg = cad)
  if (res /= 0) call error('(module_alloc_int_r3/extend) unable to allocate variable v: '//trim(cad))
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
    if (res /= 0) call error('(module_alloc_int_r3/extend) unable to allocate variable temp: '//trim(cad))
    temp = 0
    temp(1:s1,1:s2,1:s3) = v
    call move_alloc(from=temp, to=v)
  end if
end if
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce the array
!-----------------------------------------------------------------------
subroutine reduce_prv(v, d1, d2, d3)
real(real64), allocatable :: v(:,:,:), temp(:,:,:)
integer,      intent(in)  :: d1, d2, d3

if (.not. allocated(v)) then 
  call info('(module_alloc_int_r3/reduce) Variable not allocated'); return
end if
if (size(v,1) == d1 .and. size(v,2) == d2 .and. size(v,3) == d3) return !rows and cols have the right size
if (size(v,1) <  d1 .or.  size(v,2) <  d2 .or. size(v,3) < d3) then   !rows or cols are too large
  call info('(module_alloc_int_r3/reduce) Some given dimension is too large to reduce'); return
end if
call alloc(temp, d1, d2, d3)
temp(1:d1, 1:d2, 1:d3) = v(1:d1, 1:d2, 1:d3)
call move_alloc(from=temp, to=v)
end subroutine

!-----------------------------------------------------------------------
! set: set a scalar in the array
!-----------------------------------------------------------------------
subroutine set_scalar_prv(v, val, d1, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val
integer,      intent(in)  :: d1, d2, d3
logical, intent(in), optional :: fit(3)

call extend(v, d1, d2, d3, fit)
v(d1, d2, d3) = val
end subroutine

!-----------------------------------------------------------------------
! set_row: set a row in the array
!-----------------------------------------------------------------------
subroutine set_row_prv(v, val, d1, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d1, d3
logical, intent(in), optional :: fit(3)

call extend(v, d1, size(val,1), d3, fit)
v(d1, 1:size(val,1), d3) = val
end subroutine

!-----------------------------------------------------------------------
! set_col: set a column in the array
!-----------------------------------------------------------------------
subroutine set_col_prv(v, val, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d2, d3
logical, intent(in), optional :: fit(3)

call extend(v, size(val,1), d2, d3, fit)
v(1:size(val,1), d2, d3) = val
end subroutine

!-----------------------------------------------------------------------
! set: set a matrix in the array
!-----------------------------------------------------------------------
subroutine set_matrix_prv(v, val, d1, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d1(:), d2(:), d3(:)
logical, intent(in), optional :: fit(3)

call extend(v, maxval(d1), maxval(d2), maxval(d3), fit)
v(d1, d2, d3) = reshape(val,[size(d1), size(d2), size(d3)])
end subroutine

!-----------------------------------------------------------------------
! add: add a scalar in the array
!-----------------------------------------------------------------------
subroutine add_scalar_prv(v, val, d1, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val
integer,      intent(in)  :: d1, d2, d3
logical, intent(in), optional :: fit(3)

call extend(v, d1, d2, d3, fit)
v(d1, d2, d3) = v(d1, d2, d3) + val
end subroutine

!-----------------------------------------------------------------------
! add_row: add a row in the array
!-----------------------------------------------------------------------
subroutine add_row_prv(v, val, d1, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d1, d3
logical, intent(in), optional :: fit(3)

call extend(v, d1, size(val,1), d3, fit)
v(d1, 1:size(val,1), d3) = v(d1, 1:size(val,1), d3) + val
end subroutine

!-----------------------------------------------------------------------
! add_col: add a column in the array
!-----------------------------------------------------------------------
subroutine add_col_prv(v, val, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d2, d3
logical, intent(in), optional :: fit(3)

call extend(v, size(val,1), d2, d3, fit)
v(1:size(val,1), d2, d3) = v(1:size(val,1), d2, d3) + val
end subroutine

!-----------------------------------------------------------------------
! add: add a matrix in the array
!-----------------------------------------------------------------------
subroutine add_matrix_prv(v, val, d1, d2, d3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d1(:), d2(:), d3(:)
logical, intent(in), optional :: fit(3)

call extend(v, maxval(d1), maxval(d2), maxval(d3), fit)
v(d1, d2, d3) = v(d1, d2, d3) + reshape(val,[size(d1), size(d2), size(d3)])
end subroutine

!-----------------------------------------------------------------------
! insert_row: insert a row in the array
!-----------------------------------------------------------------------
subroutine insert_row_prv(v, val, d1, d3, maxd1, maxd3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d1, d3
integer, intent(in), optional :: maxd1, maxd3
logical, intent(in), optional :: fit(3)
integer :: s1, s3

if (present(maxd1)) then; s1 = max(maxd1+1, d1)
else; s1 = max(size(v,1)+1, d1)
end if
if (present(maxd3)) then; s3 = max(maxd3+1, d3)
else; s3 = max(size(v,1)+1, d3)
end if
call extend(v, s1, size(val,1), s3, fit)
v(d1+1:size(v,1), : ,d3+1:size(v,3)) = v(d1:size(v,1)-1, : , d3:size(v,3)-1)
v(d1, 1:size(val,1), d3) = val
v(d1, size(val,1)+1:size(v,2), d3) = 0
end subroutine

!-----------------------------------------------------------------------
! insert_col: insert a col in the array
!-----------------------------------------------------------------------
subroutine insert_col_prv(v, val, d2, d3, maxd2, maxd3, fit)
real(real64), allocatable :: v(:,:,:)
real(real64), intent(in)  :: val(:)
integer,      intent(in)  :: d2, d3
integer, intent(in), optional :: maxd2, maxd3
logical, intent(in), optional :: fit(3)
integer :: s2, s3

if (present(maxd2)) then; s2 = max(maxd2+1, d2)
else; s2 = max(size(v,2)+1, d2)
end if
if (present(maxd3)) then; s2 = max(maxd3+1, d3)
else; s3 = max(size(v,3)+1, d3)
end if
call extend(v, size(val,1), s2, s3, fit)
v(:, d2+1:size(v,2), d3+1:size(v,3)) = v(:, d2:size(v,2)-1, d3:size(v,3)-1)
v(1:size(val,1), d2, d3) = val
v(size(val,1)+1:size(v,1), d2, d3) = 0
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FALTA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!-----------------------------------------------------------------------
!! insert_row_sorted: insert a row in a row-sorted array
!!-----------------------------------------------------------------------
!subroutine insert_row_sorted_prv(v, val, used, fit, pos)
!integer, allocatable             :: v(:,:)
!integer,           intent(in)    :: val(:)
!integer, optional, intent(inout) :: used
!logical, optional, intent(in)    :: fit(2)
!integer, optional, intent(out)   :: pos
!integer :: n, a, b, anew, bnew, i, j,indx
!
!!v not allocated
!if (.not. allocated(v)) then
!  indx = -1
!  call set_row(v, val, -indx, fit)
!  if (present(used)) used = 1
!  if (present(pos)) pos = indx
!  return
!end if
!!number of existing rows
!if (present(used)) then; n = used
!else;                    n = size(v,1)
!end if
!!search among the first vertices
!indx = bsearch(v(1:n,1), val(1), n)
!if (indx < 0) then
!! insert and return
!  call insert_row(v, val, -indx, maxrow=n, fit=fit)
!  if (present(used)) used = n+1
!  if (present(pos)) pos = indx
!  return
!end if
!a=1; b = n
!do j = 2, size(val,1)
!! determine the left extreme of the interval where to search the j-th vertex
!  anew = indx
!  do i = indx-1, a, -1
!    if (v(i,j-1) /= v(indx,j-1)) exit
!    anew = i
!  end do
!! determine the right extreme of the interval where to search the j-th vertex
!  bnew = indx
!  do i = indx+1, b
!    if (v(i,j-1) /= v(indx,j-1)) exit
!    bnew = i
!  end do
!  a = anew; b = bnew
!  indx = bsearch(v(a:b,j), val(j), b-a+1)
!  if (indx < 0) then
!    indx = indx - a+1
!!   insert and return
!    call insert_row(v, val, -indx, maxrow=n, fit=fit)
!    if (present(used)) used = n+1
!    if (present(pos)) pos = indx
!    return
!  else
!    indx = indx + a-1
!  end if
!end do
!if (present(pos)) pos = indx
!end subroutine
!
!!-----------------------------------------------------------------------
!! insert_col_sorted: insert a col in a col-sorted array
!!-----------------------------------------------------------------------
!subroutine insert_col_sorted_prv(v, val, used, fit, pos)
!integer, allocatable             :: v(:,:,:)
!integer,           intent(in)    :: val(:)
!integer, optional, intent(inout) :: used
!logical, optional, intent(in)    :: fit(2)
!integer, optional, intent(out)   :: pos
!integer :: n, a, b, anew, bnew, i, j,indx
!
!!v not allocated
!if (.not. allocated(v)) then
!  indx = -1
!  call set_col(v, val, -indx, fit)
!  if (present(used)) used = 1
!  if (present(pos)) pos = indx
!  return
!end if
!!number of existing rows
!if (present(used)) then; n = used
!else;                    n = size(v,2)
!end if
!!search among the first vertices
!indx = bsearch(v(1,1:n), val(1), n)
!if (indx < 0) then
!! insert and return
!  call insert_col(v, val, -indx, maxcol=n, fit=fit) 
!  if (present(used)) used = n+1
!  if (present(pos)) pos = indx
!  return
!end if
!a=1; b = n
!do j = 2, size(val,1)
!! determine the left extreme of the interval where to search the j-th vertex
!  anew = indx
!  do i = indx-1, a, -1
!    if (v(j-1,i) /= v(j-1,indx)) exit
!    anew = i
!  end do
!! determine the right extreme of the interval where to search the j-th vertex
!  bnew = indx
!  do i = indx+1, b
!    if (v(j-1,i) /= v(j-1,indx)) exit
!    bnew = i
!  end do
!  a = anew; b = bnew
!  indx = bsearch(v(j,a:b), val(j), b-a+1)
!  if (indx < 0) then
!    indx = indx - a+1
!!   insert and return
!    call insert_col(v, val, -indx, maxcol=n, fit=fit)
!    if (present(used)) used = n+1
!    if (present(pos)) pos = indx
!    return
!  else
!    indx = indx + a-1
!  end if
!end do
!if (present(pos)) pos = indx
!end subroutine
!
!!-----------------------------------------------------------------------
!! find_row_sorted: find the position of a row in a row-sorted array
!!-----------------------------------------------------------------------
!function find_row_sorted_prv(v, val, used) result(pos)
!integer, allocatable          :: v(:,:)
!integer,           intent(in) :: val(:)
!integer, optional, intent(in) :: used
!integer :: n, pos, a, b, anew, bnew, i, j
!
!pos = -1
!if (.not. allocated(v)) return
!!number of existing rows
!if (present(used)) then; n = used
!else;                    n = size(v,1)
!end if
!!search among the first vertices
!pos = bsearch(v(1:n,1), val(1), n)
!if (pos <= 0) return
!a=1; b = n
!do j = 2, size(val,1)
!! determine the left extreme of the interval where to search the j-th vertex
!  anew = pos
!  do i = pos-1, a, -1
!    if (v(i,j-1) /= v(pos,j-1)) exit
!    anew = i
!  end do
!! determine the right extreme of the interval where to search the j-th vertex
!  bnew = pos
!  do i = pos+1, b
!    if (v(i,j-1) /= v(pos,j-1)) exit
!    bnew = i
!  end do
!  a = anew; b = bnew
!  pos = bsearch(v(a:b,j), val(j), b-a+1)
!  if (pos <= 0) return
!  pos = pos + a-1
!end do
!end function
!
!!-----------------------------------------------------------------------
!! find_col_sorted: find the position of a col in a col-sorted array
!!-----------------------------------------------------------------------
!function find_col_sorted_prv(v, val, used) result(pos)
!integer, allocatable          :: v(:,:)
!integer,           intent(in) :: val(:)
!integer, optional, intent(in) :: used
!integer :: n, pos, a, b, anew, bnew, i, j
!
!pos = -1
!if (.not. allocated(v)) return
!!number of existing cols
!if (present(used)) then; n = used
!else;                    n = size(v,2)
!end if
!!search among the first vertices
!pos = bsearch(v(1,1:n), val(1), n)
!if (pos <= 0) return
!a=1; b = n
!do j = 2, size(val,1)
!! determine the left extreme of the interval where to search the j-th vertex
!  anew = pos
!  do i = pos-1, a, -1
!    if (v(j-1,i) /= v(j-1,pos)) exit
!    anew = i
!  end do
!! determine the right extreme of the interval where to search the j-th vertex
!  bnew = pos
!  do i = pos+1, b
!    if (v(j-1,i) /= v(j-1,pos)) exit
!    bnew = i
!  end do
!  a = anew; b = bnew
!  pos = bsearch(v(j,a:b), val(j), b-a+1)
!  if (pos <= 0) return
!  pos = pos + a-1
!end do
!end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FALTA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------------
! sfind (subroutine): find all the occurrences of scalar value in an array
!-----------------------------------------------------------------------
subroutine sfind_sca_prv(v, val, d1, d2, d3)
real(real64),         intent(in)  :: v(:,:,:)            ! array
real(real64),         intent(in)  :: val                 ! value to search
integer, allocatable, intent(out) :: d1(:), d2(:), d3(:) ! results
integer :: i, j, k, n, p

! allocate row, col
n = size(pack(v,v==val),1)
call alloc(d1, n)
call alloc(d2, n)
call alloc(d3, n)
! find positions
p = 1
do i = 1, size(v,1)
  do j = 1, size(v,2)
    do k = 1, size(v,3)
      if (v(i,j,k) == val) then
        d1(p) = i; d2(p) = j; d3(p) = k
        p = p+1
      end if  
      if (p > n) return
    enddo
  end do  
end do
end subroutine

!-----------------------------------------------------------------------
! sfind (subroutine): find all the occurrences of an array of values in an array
!-----------------------------------------------------------------------
subroutine sfind_vec_prv(v, val, d1, d2, d3)
real(real64),         intent(in)  :: v(:,:,:)            ! array
real(real64),         intent(in)  :: val(:)              ! values to search
integer, allocatable, intent(out) :: d1(:), d2(:), d3(:) ! results
integer :: i, j, k, l, n, p

! allocate row, col
n = 0
do j = 1, size(val,1)
  n = n + size(pack(v,v==val(j)),1)
end do  
call alloc(d1, n)
call alloc(d2, n)
call alloc(d3, n)
! find positions
p = 1
do i = 1, size(v,1)
  do k = 1, size(v,2)
    do l = 1, size(v,3)
      do j = 1, size(val,1)
        if (v(i,k,l) == val(j)) then
          d1(p) = i; d2(p) = k; d3(p) = l
          p = p+1
        end if  
      end do  
      if (p > n) return
    enddo
  end do
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

module module_groups_fcnv
!-----------------------------------------------------------------------
! Module for group management
! Last update: 04/04/2010
!
! groups is a matrix where:
!  - groups(d+1)%ncells is the number of cells that have dimension d
!  - groups(d+1)%cell(:,1:d+1) are the numbers of vertices of the cells
!  - groups(d+1)%cell(:,  d+2) are the reference numbers of the cells
!-----------------------------------------------------------------------
use basicmod
implicit none

!Types
type :: groups_type
  integer                              :: ncells = 0
  integer, allocatable, dimension(:,:) :: cell
end type

!Variables
type(groups_type), public, allocatable, dimension(:) :: groups

!Interface
interface alloc;  module procedure alloc_prv;     end interface
interface insert; module procedure insert_groups; end interface

!Private procedures
private :: binary_search, insert_prv, alloc_prv, insert_groups

contains

!-----------------------------------------------------------------------
! alloc: alloc memory
!-----------------------------------------------------------------------
subroutine alloc_prv(v, n)

  type(groups_type), dimension(:), allocatable :: v
  integer, intent(in)                          :: n
  integer :: res, i

  if (allocated(v)) then
    do i = 1, size(v,1)
      if (allocated(v(i)%cell)) then
        deallocate(v(i)%cell, stat=res)
        if (res /= 0) call error('(module_groups/alloc/deallocate) Unable to deallocate variable v(i)%cell')
      end if
    end do
    deallocate(v, stat = res)
    if (res /= 0) call error('(module_groups/alloc/deallocate) Unable to deallocate variable v')
  end if
  allocate(v(n), stat = res) !In f2003 add: errmsg = cad
  if (res /= 0) call error('(module_groups/alloc/allocate) Unable to allocate variable v')

end subroutine

!-----------------------------------------------------------------------
! search
!-----------------------------------------------------------------------
function search(this, v) result(pos)
type (groups_type),    intent(in) :: this !group of entities
integer, dimension(:), intent(in) :: v    !vector to find
integer :: pos, a, b, anew, bnew, i, j

!search among the first vertices
pos = binary_search(this%ncells, this%cell(1:this%ncells,1), v(1))
if (pos <= 0) return
a=1; b = this%ncells
do j = 2, size(v,1)
! determine the left extreme of the interval where to search the j-th vertex
  anew = pos
  do i = pos-1, a, -1
    if (this%cell(i,j-1) /= this%cell(pos,j-1)) exit
    anew = i
  end do
! determine the right extreme of the interval where to search the j-th vertex
  bnew = pos
  do i = pos+1, b
    if (this%cell(i,j-1) /= this%cell(pos,j-1)) exit
    bnew = i
  end do
  a = anew; b = bnew
  pos = binary_search(b-a+1, this%cell(a:b,j), v(j))
  if (pos <= 0) return
  pos = pos + a-1
end do
end function

!-----------------------------------------------------------------------
! insert_groups
!-----------------------------------------------------------------------
subroutine insert_groups(this, v)
type (groups_type),    intent(inout) :: this !group of entities
integer, dimension(:), intent(in)    :: v    !vector to find
integer :: pos, a, b, anew, bnew, i, j

!cell not allocated
if (.not. allocated(this%cell)) then
  call insert_prv(this%ncells, this%cell, v, 1)
  return
end if

!search among the first vertices
pos = binary_search(this%ncells, this%cell(1:this%ncells,1), v(1))
if (pos < 0) then
! insert and return
  call insert_prv(this%ncells, this%cell, v, -pos)
  return
end if
a=1; b = this%ncells
do j = 2, size(v,1)
! determine the left extreme of the interval where to search the j-th vertex
  anew = pos
  do i = pos-1, a, -1
    if (this%cell(i,j-1) /= this%cell(pos,j-1)) exit
    anew = i
  end do
! determine the right extreme of the interval where to search the j-th vertex
  bnew = pos
  do i = pos+1, b
    if (this%cell(i,j-1) /= this%cell(pos,j-1)) exit
    bnew = i
  end do
  a = anew; b = bnew
  pos = binary_search(b-a+1, this%cell(a:b,j), v(j))
  if (pos < 0) then
    pos = -pos + a-1
!   insert and return
    call insert_prv(this%ncells, this%cell, v, pos)
    return
  else
    pos = pos + a-1
  end if
end do

end subroutine

!***********************************************************************
! PRIVATE PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! binary_search: If 'u' is in 'x' then its index is returned.
! Otherwise − ind is returned, where ind is the index where 'u'
! should be inserted in 'x'
!-----------------------------------------------------------------------
function binary_search(n, x, u) result(pos)
integer,               intent(in) :: n !components of seq
integer, dimension(:), intent(in) :: x !array, qhere to search
integer,               intent(in) :: u !value to search
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

!  right = n
!  left = 1
!  previous_center = 0
!  do
!    center = ceiling((left + right) / 2.)
!    candidate = seq(center)
!    if (search == candidate) then
!      pos = center; return
!    end if
!    if (center == previous_center) then
!      pos = -center-1; return
!    elseif (search < candidate) then
!      right = center
!    else
!      left = center
!    end if
!    previous_center = center
!  end do
end function

!-----------------------------------------------------------------------
! insert_prv
!-----------------------------------------------------------------------
subroutine insert_prv(n, cell, v, row)
integer, intent(inout) :: n !total rows in cell
integer, allocatable, dimension(:,:) :: cell !group cells
integer, dimension(:), intent(in) :: v !new cell
integer, intent(in) :: row !position for new cell

! enlarge cell if necessary
call extend(cell, n+1, size(v,1), fit=[.true., .true.])
! copy to an upper position
cell((row+1):(n+1),:) = cell(row:n,:)
!insert the new cell
cell(row,:) = v
n = n + 1
end subroutine

end module

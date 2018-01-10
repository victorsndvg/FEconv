module module_manage_unv_fcnv
!-----------------------------------------------------------------------
! Module for UNV file management
! Last update: 04/04/2010
!-----------------------------------------------------------------------
use basicmod
!use module_mesh
use module_pmh_fcnv, only: pmh_mesh
use module_dataset_2411_fcnv
use module_dataset_2412_fcnv
use module_dataset_2467_fcnv
use module_dataset_2414_fcnv
implicit none

!Types
type unv
  private
  character(maxpath) :: filename = ' ' !file name
  integer            :: UNIT     = -1  !associated unit number
end type

!Private module procedures
private :: search_dataset_delimiter, search_dataset_type, create_vertex_data

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! open: open universal file
!-----------------------------------------------------------------------
subroutine open_unv(this, filename)
type(unv)                    :: this     !unv object
character(len=*), intent(in) :: filename !unv file
integer :: ios

! open file
this%filename = filename
this%unit = get_unit()
open (unit=this%unit, file=this%filename, form='formatted', iostat=ios, &
status='old', position='rewind')
if (ios /= 0) call error('unv/open, #'//trim(string(ios)))

end subroutine

!-----------------------------------------------------------------------
! read: read universal file
!-----------------------------------------------------------------------
subroutine read_unv(this, pmh, padval, infield, ca_opt)
type(unv),                intent(in)    :: this       !universal object
type(pmh_mesh),           intent(inout) :: pmh        !mesh
real(real64),             intent(in)    :: padval     !field padding value
character(*),allocatable, intent(in)    :: infield(:) !input field names
logical,                  intent(in)    :: ca_opt     !-ca Code Aster option
integer :: maxdim !dimension detected
integer, allocatable :: eloc(:,:) !elements location in pmh structure
integer :: ios, n, i, pgroup(6), fgroup(4),nfield

if(.not. allocated(pmh%pc)) allocate(pmh%pc(1))
n = size(pmh%pc,1)

!dataset 2411, node coordinates
rewind(unit=this%unit, iostat=ios)
if (ios /= 0) call error('unv/read/rewind, #'//trim(string(ios)))
if (search_dataset_type(this,2411) /= 0) then
  if (search_dataset_type(this,781) /= 0) then
    call error('unv/read, dataset 2411 or 781 not found')
  else
    continue
  end if
end if
call read_2411(this%unit, pmh%pc(n))

!dataset 2412, elements
rewind(unit=this%unit, iostat=ios)
if (ios /= 0) call error('unv/read/rewind, #'//trim(string(ios)))
if (search_dataset_type(this,2412) /= 0) call error('unv/read, dataset 2412 not found')
call read_2412(this%unit, pmh%pc(n), maxdim, eloc)

!datasets for permanent groups, several numbers (see pgroup)
pgroup = [2467, 2477, 2452, 2435, 2432, 2430]
do i = 1, size(pgroup,1)
  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('unv/read/rewind, #'//trim(string(ios)))
  if (search_dataset_type(this,pgroup(i)) == 0) then
    call read_2467(this%unit, pmh%pc(n), eloc)
    exit
  end if
end do

!datasets for fields, several numbers (see fgroup)
fgroup = [2414, 55, 56, 57]
nfield = 0
do i = 1, size(fgroup,1)
  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('unv/read/rewind, #'//trim(string(ios)))
  do while (search_dataset_type(this,fgroup(i)) == 0)
    call read_2414(this%unit, pmh, n, nfield, eloc, fgroup(i), infield, ca_opt,  padval )
    !exit not needed, reading multiple fields
  end do
enddo

!close file
close(unit=this%unit, iostat=ios)
if (ios /= 0) call error('univ_file/close, #'//trim(string(ios)))
end subroutine

!***********************************************************************
! PRIVATE  PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! search_dataset_delimiter: searchs the next dataset delimiter
! RETURN:  0 if the delimiter is found
!         <0 if the end of file is reached
!         >0 if an error occurs
!-----------------------------------------------------------------------
function search_dataset_delimiter(this) result(res)
type(unv), intent(in) :: this !universal file
integer :: res, val

do
  read (unit=this%unit, fmt='(I6)', iostat=res) val
  !find EOF or '-1', then return
  if (res<0 .or. val==-1) return
  !res>0 when reading chars, then cycle
end do
end function

!-----------------------------------------------------------------------
! search_dataset_type: searchs a specific dataset type
! RETURN:  0 if the delimiter is found
!         <0 if the end of file is reached
!         >0 if an error occurs
!-----------------------------------------------------------------------
function search_dataset_type(this, ds) result(res)
type(unv), intent(in) :: this !universal object
integer,   intent(in) :: ds   !dataset type
integer :: res, dsnumber

do
  res = search_dataset_delimiter(this)
  if (res /= 0) return
  !read the type of dataset
  read (unit=this%unit, fmt='(I6)', iostat=res) dsnumber
  if (res /= 0) call error('unv/search_dataset_type, #'//trim(string(res)))
  if (dsnumber == ds) exit
  res = search_dataset_delimiter(this)
  if (res /= 0) return
enddo
end function

!-----------------------------------------------------------------------
! create_vertex_data (unused): create iv and xv from id and xd
!-----------------------------------------------------------------------
subroutine create_vertex_data(m)
type(mfm_mesh), intent(inout) :: m !mesh
integer, allocatable :: vert2node(:), node2vert(:)
integer :: nv2d, pos, maxv, i, j, k

if (m%LNN /= m%LNV) then
  !construct vert2node
  nv2d = 0
  do k = 1, m%nl
    do j = 1, m%LNV
      if (m%id(j,k) == 0) cycle
      pos = bsearch(vert2node, m%id(j,k), nv2d)
      if (pos < 0) then
        call insert(vert2node, m%id(j,k), -pos, nv2d, fit=.false.)
        nv2d = nv2d + 1
      endif
    end do
  end do
  call reduce(vert2node, nv2d)
  m%nv = nv2d
  !construct node2vert
  maxv = 0
  do i = 1, size(vert2node, 1)
    maxv = max(maxv, vert2node(i))
  enddo
  allocate(node2vert(maxv))
  node2vert = 0
  do i = 1, size(vert2node, 1)
    if (vert2node(i) == 0) cycle
    node2vert(vert2node(i)) = i
  enddo
  !vertices renumbering
  call alloc(m%iv, m%LNV, m%nl)
  do k = 1, m%nl
    do j = 1, m%LNV
      m%iv(j,k) = node2vert(m%id(j,k))
    end do
  end do
  !correct dimension
  if (any(abs(m%xd(3,:))>1e3*epsilon(m%xd))) m%DIM = 3
  print'(a,i1)', 'Maximum detected dimension: ',m%DIM
  !vertices coordinates
  call alloc(m%xv, m%DIM, m%nv)
  do i = 1, m%nv
    m%xv(1:m%DIM,i) = m%xd(1:m%DIM, vert2node(i))
  end do
else
   m%nv = m%nd
   call alloc(m%iv, m%LNV, m%nl)
   m%iv(1:m%LNV,1:m%nl) = m%id(1:m%LNV,1:m%nl)
   if (any(abs(m%xd(3,:))>1e3*epsilon(m%xd))) m%DIM = 3
   print'(a,i1)', 'Maximum detected dimension: ',m%DIM
   !vertices coordinates
   call alloc(m%xv, m%DIM, m%nv)
   m%xv(1:m%DIM,1:m%nv) = m%xd(1:m%DIM,1:m%nv)
end if
end subroutine

end module

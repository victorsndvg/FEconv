module module_mum_fcnv
!-----------------------------------------------------------------------
! Module to manage MUM (Modulef Unformatted Mesh) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2013
!
! PUBLIC PROCEDURES:
!   load_mum: loads a mesh from a MUM format file
!   save_mum: saves a mesh in a MUM format file
!-----------------------------------------------------------------------
use basicmod, only: real64, maxpath, error, string, int, feed, empty
implicit none

contains

!-----------------------------------------------------------------------
! load_mum: mesh load
!-----------------------------------------------------------------------
subroutine load_mum(filename, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, &
nn, mm, nrc, nra, nrv, z, nsd)
character(*), intent(in) :: filename !mesh filename
integer, intent(in)    :: iu         !file unit
integer, intent(inout) :: nel        !global number of elements
integer, intent(inout) :: nnod       !global number of nodes
integer, intent(inout) :: nver       !global number of vertices
integer, intent(inout) :: dim        !space dimension
integer, intent(inout) :: lnv        !local number of vertices
integer, intent(inout) :: lne        !local number of edges
integer, intent(inout) :: lnf        !local number of faces
integer, intent(inout) :: lnn        !local number of nodes
integer, allocatable   :: nn(:,:)    !nodes index array
integer, allocatable   :: mm(:,:)    !vertices index array
integer, allocatable   :: nrv(:,:)   !vertices reference array
integer, allocatable   :: nra(:,:)   !edge reference array
integer, allocatable   :: nrc(:,:)   !face reference array
real(real64), allocatable :: z(:,:)  !vertices coordinates array
integer, allocatable   :: nsd(:)     !subdomain index array
integer :: i, j, k, ln2, ios

!open file
open (unit=iu, file=filename, form='unformatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('load/open, #'//trim(string(ios)))

!try read scalar data from file: nel, nnod, nver, dim, lnn, lnv, lne, lnf
read (unit=iu, iostat=ios) nel, nnod, nver, dim, lnn, lnv, lne, lnf
if (ios /= 0) then
  !try read scalar data from file: nel, nnod, nver
  backspace(unit=iu, iostat=ios)
  if (ios /= 0) call error('backspace, #'//trim(string(ios)))
  read (unit=iu, iostat=ios) nel, nnod, nver
  if (ios /= 0) call error('read (nel, nnod, nver), #'//trim(string(ios)))
  if (dim <= 0 .or. lnn <= 0 .or. lnv <= 0) &
  call error('RECONVXX format does not include dimensions and the user did not provide them, unable to allocate')
end if

!allocation
if (allocated(mm))  then
deallocate( mm,           stat=ios); if (ios /= 0) call error('dealloc (mm), unable to deallocate');  end if
  allocate( mm(lnv, nel), stat=ios); if (ios /= 0) call error(  'alloc (mm), unable to allocate')
if (allocated(nn))  then
deallocate( nn,           stat=ios); if (ios /= 0) call error('dealloc (nn), unable to deallocate');  end if
  allocate( nn(lnn, nel), stat=ios); if (ios /= 0) call error(  'alloc (nn), unable to allocate')
if (allocated(nrc)) then
deallocate(nrc,           stat=ios); if (ios /= 0) call error('dealloc (nrc), unable to deallocate'); end if
  allocate(nrc(lnf, nel), stat=ios); if (ios /= 0) call error(  'alloc (nrc), unable to allocate')
if (allocated(nra)) then
deallocate(nra,           stat=ios); if (ios /= 0) call error('dealloc (nra), unable to deallocate'); end if
  allocate(nra(lne, nel), stat=ios); if (ios /= 0) call error(  'alloc (nra), unable to allocate')
if (allocated(nrv)) then
deallocate(nrv,           stat=ios); if (ios /= 0) call error('dealloc (nrv), unable to deallocate'); end if
  allocate(nrv(lnv, nel), stat=ios); if (ios /= 0) call error(  'alloc (nrv), unable to allocate')
if (allocated(z))   then
deallocate( z,            stat=ios); if (ios /= 0) call error('dealloc (z), unable to deallocate');   end if
  allocate( z(dim, nver), stat=ios); if (ios /= 0) call error(  'alloc (z), unable to allocate')
if (allocated(nsd)) then
deallocate(nsd,           stat=ios); if (ios /= 0) call error('dealloc (nsd), unable to deallocate'); end if
  allocate(nsd(nel),      stat=ios); if (ios /= 0) call error(  'alloc (nsd), unable to allocate')
!ln2: lecture of nodes, only if (nnod /= nver)
ln2 = 0; if (nnod /= nver) ln2 = lnn
!arrays from file ([nn,if nnod/=nver], mm, [nrc,if dim==3], [nra,if dim==2], nrv, x)
select case(dim)
case(1)
  read (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
                             ((mm(i,k),  i=1,lnv), k=1,nel), &
                             ((nrv(i,k), i=1,lnv), k=1,nel), &
                             ((z(i,j),   i=1,dim), j=1,nver)
case(2)
  read (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
                             ((mm(i,k),  i=1,lnv), k=1,nel), &
                             ((nra(i,k), i=1,lne), k=1,nel), &
                             ((nrv(i,k), i=1,lnv), k=1,nel), &
                             ((z(i,j),   i=1,dim), j=1,nver)
case(3)
  read (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
                             ((mm(i,k),  i=1,lnv), k=1,nel), &
                             ((nrc(i,k), i=1,lnf), k=1,nel), &
                             ((nra(i,k), i=1,lne), k=1,nel), &
                             ((nrv(i,k), i=1,lnv), k=1,nel), &
                             ((z(i,j),   i=1,dim), j=1,nver)
case default
  call error('load, invalid dimension DIM')
end select
if (ios /= 0) call error('read (mm,...) #'//trim(string(ios)))
!subdomain index array (nsd)
read (unit=iu, iostat=ios) (nsd(k), k=1,nel)
if (ios /= 0) call error('read (nsd), #'//trim(string(ios)))
close(iu)
! nn was read only if (nnod /= nver)
!if (nnod == nver) nn(1:lnv,1:nel) = mm

print'(a,i9)','MUM file loaded!'
print'(a,i9)','Global number of elements: ', nel
print'(a,i9)','Global number of nodes:    ', nnod
print'(a,i9)','Global number of vertices: ', nver
print'(a,i9)','Space dimension :          ', dim
print'(a,i9)','Local number of nodes :    ', lnn
print'(a,i9)','Local number of vertices : ', lnv
print'(a,i9)','Local number of edges :    ', lne
print'(a,i9)','Local number of faces :    ', lnf
end subroutine

!-----------------------------------------------------------------------
! save_mum: save mesh
!-----------------------------------------------------------------------
subroutine save_mum(filename, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, &
nn, mm, nrc, nra, nrv, z, nsd)
character(*), intent(in) :: filename !mesh filename
integer, intent(in) :: iu            !file unit
integer, intent(in) :: nel           !global number of elements
integer, intent(in) :: nnod          !global number of nodes
integer, intent(in) :: nver          !global number of vertices
integer, intent(in) :: dim           !space dimension
integer, intent(in) :: lnv           !local number of vertices
integer, intent(in) :: lne           !local number of edges
integer, intent(in) :: lnf           !local number of faces
integer, intent(in) :: lnn           !local number of nodes
integer, intent(in) :: nn(:,:)       !nodes index array
integer, intent(in) :: mm(:,:)       !vertices index array
integer, intent(in) :: nrv(:,:)      !vertices reference array
integer, intent(in) :: nra(:,:)      !edge reference array
integer, intent(in) :: nrc(:,:)      !face reference array
real(real64), intent(in) :: z(:,:)   !vertices coordinates array
integer, intent(in) :: nsd(:)        !subdomain index array
integer :: i, j, k, ln2, ios

!open file
open (unit=iu, file=filename, form='unformatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))
!write data (nel, nnod, nver, dim, ...)
write (unit=iu, iostat=ios) nel, nnod, nver, dim, lnn, lnv, lne, lnf
if (ios /= 0) call error('write (nel,...), #'//trim(string(ios)))
!ln2: write nodes, only if (nnod /= nver)
ln2 = lnn; if (nnod == nver) ln2 = 0
!arrays from file ([nn,if nnod/=nver], mm, [nrc,if dim==3], [nra,if dim==2], nrv, x)

select case(dim)
case(1)

write (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
((mm(i,k),  i=1,lnv), k=1,nel), &
((nrv(i,k), i=1,lnv), k=1,nel), &
((z(i,j),   i=1,dim), j=1,nver)
case(2)
write (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
((mm(i,k),  i=1,lnv), k=1,nel), &
((nra(i,k), i=1,lne), k=1,nel), &
((nrv(i,k), i=1,lnv), k=1,nel), &
((z(i,j),   i=1,dim), j=1,nver)
case(3)
write (unit=iu, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
((mm(i,k),  i=1,lnv), k=1,nel), &
((nrc(i,k), i=1,lnf), k=1,nel), &
((nra(i,k), i=1,lne), k=1,nel), &
((nrv(i,k), i=1,lnv), k=1,nel), &
((z(i,j),   i=1,dim), j=1,nver)
case default
call error('save, invalid dimension DIM')
end select
if (ios /= 0) call error('write (mm,...) #'//trim(string(ios)))
!subdomain index array (nsd)
write (unit=iu, iostat=ios) (nsd(k), k=1,nel)
if (ios /= 0) call error('write (nsd), #'//trim(string(ios)))
close(iu)
end subroutine

end module

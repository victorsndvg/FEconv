module module_mfm_fcnv
!-----------------------------------------------------------------------
! Module to manage MFM (Modulef Formatted Mesh) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 15/05/2013
!
! PUBLIC PROCEDURES:
!   load_mfm: loads a mesh from a MFM format file
!   save_mfm: saves a mesh in a MFM format file
!-----------------------------------------------------------------------
use basicmod, only: real64, maxpath, error, string, int, feed, empty, info
implicit none

contains

!-----------------------------------------------------------------------
! load_mfm: mesh load
!-----------------------------------------------------------------------
subroutine load_mfm(filename, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, &
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
character(maxpath) :: str

!open file
open (unit=iu, file=filename, form='formatted', status='old', position='rewind', iostat=ios)
if (ios /= 0) call error('load/open, #'//trim(string(ios)))

!try read scalar data from file: nel, nnod, nver, dim, lnn, lnv, lne, lnf
read (unit=iu, fmt='(a)', iostat=ios) str
if (ios /= 0) call error('read (str), #'//trim(string(ios)))
read (str, fmt=*, iostat=ios) nel, nnod, nver, dim, lnn, lnv, lne, lnf
if (ios /= 0) then
!try read scalar data from file: nel, nnod, nver
read (str, fmt=*, iostat=ios) nel, nnod, nver
if (ios /= 0) call error('read (nel, nnod, nver), #'//trim(string(ios)))
if (dim <= 0 .or. lnn <= 0 .or. lnv <= 0) &
call error('RECONVXX format does not include dimensions and the user did not provide them, unable to allocate')
end if

!allocation
if (allocated(nn))  then
deallocate( nn,           stat=ios); if (ios /= 0) call error('dealloc (nn), unable to deallocate');  end if
allocate( nn(lnn, nel), stat=ios); if (ios /= 0) call error(  'alloc (nn), unable to allocate')
if (allocated(mm))  then
deallocate( mm,           stat=ios); if (ios /= 0) call error('dealloc (mm), unable to deallocate');  end if
allocate( mm(lnv, nel), stat=ios); if (ios /= 0) call error(  'alloc (mm), unable to allocate')
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
read (unit=iu, fmt=*, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
((mm(i,k),  i=1,lnv), k=1,nel), &
((nrv(i,k), i=1,lnv), k=1,nel), &
((z(i,j),   i=1,dim), j=1,nver)
case(2)
read (unit=iu, fmt=*, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
((mm(i,k),  i=1,lnv), k=1,nel), &
((nra(i,k), i=1,lne), k=1,nel), &
((nrv(i,k), i=1,lnv), k=1,nel), &
((z(i,j),   i=1,dim), j=1,nver)
case(3)
read (unit=iu, fmt=*, iostat=ios) ((nn(i,k),  i=1,ln2), k=1,nel), &
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
read (unit=iu, fmt=*, iostat=ios) (nsd(k), k=1,nel)
if (ios /= 0) call error('read (nsd), #'//trim(string(ios)))
close(iu)
! nn was read only if (nnod /= nver)
!if (nnod == nver) nn(1:lnv,1:nel) = mm

call info('(feconv::load_mfm) File '//trim(filename)//' loaded!')
call info('(feconv::load_mfm) Global number of elements: '//trim(string(nel)))
call info('(feconv::load_mfm) Global number of nodes:    '//trim(string(nnod)))
call info('(feconv::load_mfm) Global number of vertices: '//trim(string(nver)))
call info('(feconv::load_mfm) Space dimension :          '//trim(string(dim)))
call info('(feconv::load_mfm) Local number of nodes :    '//trim(string(lnn)))
call info('(feconv::load_mfm) Local number of vertices : '//trim(string(lnv)))
call info('(feconv::load_mfm) Local number of edges :    '//trim(string(lne)))
call info('(feconv::load_mfm) Local number of faces :    '//trim(string(lnf)))
end subroutine

!-----------------------------------------------------------------------
! save_mfm: save mesh
!-----------------------------------------------------------------------
subroutine save_mfm(filename, iu, nel, nnod, nver, dim, lnn, lnv, lne, lnf, &
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
integer :: i, j, k, ln2, le2, lf2, ios

!open file
open (unit=iu, file=filename, form='formatted', position='rewind', iostat=ios)
if (ios /= 0) call error('save/open, #'//trim(string(ios)))
print'(a,i9)','Writing data ...'
!write data (nel, nnod, nver, dim, ...)
call feed(iu, string(nel)); call feed(iu, string(nnod)); call feed(iu, string(nver))
call feed(iu, string(dim)); call feed(iu, string(lnn));  call feed(iu, string(lnv))
call feed(iu, string(lne)); call feed(iu, string(lnf));  call empty(iu)
!ln2: write nodes, only if (nnod /= nver)
ln2 = lnn; if (nnod == nver) ln2 = 0
le2 = lne; if (dim < 2)      le2 = 0
lf2 = lnf; if (dim < 3)      lf2 = 0
!arrays from file ([nn,if nnod/=nver], mm, [nrc,if dim==3], [nra,if dim==2], nrv, x)
do k = 1, nel;  do i = 1, ln2; call feed(iu, string(nn(i,k)));  end do; end do
do k = 1, nel;  do i = 1, lnv; call feed(iu, string(mm(i,k)));  end do; end do
do k = 1, nel;  do i = 1, lf2; call feed(iu, string(nrc(i,k))); end do; end do
do k = 1, nel;  do i = 1, le2; call feed(iu, string(nra(i,k))); end do; end do
do k = 1, nel;  do i = 1, lnv; call feed(iu, string(nrv(i,k))); end do; end do
do j = 1, nver; do i = 1, dim; call feed(iu, string(z(i,j)));   end do; end do
call empty(iu)
!subdomain index array (nsd)
do k = 1, nel; call feed(iu, string(nsd(k)));  end do
call empty(iu)
close(iu)

call info('(feconv::load_mfm) File '//trim(filename)//' saved!')
call info('(feconv::load_mfm) Global number of elements: '//trim(string(nel)))
call info('(feconv::load_mfm) Global number of nodes:    '//trim(string(nnod)))
call info('(feconv::load_mfm) Global number of vertices: '//trim(string(nver)))
call info('(feconv::load_mfm) Space dimension :          '//trim(string(dim)))
call info('(feconv::load_mfm) Local number of nodes :    '//trim(string(lnn)))
call info('(feconv::load_mfm) Local number of vertices : '//trim(string(lnv)))
call info('(feconv::load_mfm) Local number of edges :    '//trim(string(lne)))
call info('(feconv::load_mfm) Local number of faces :    '//trim(string(lnf)))
end subroutine

end module

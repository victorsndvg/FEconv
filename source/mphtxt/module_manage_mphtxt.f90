module module_manage_mphtxt
!-----------------------------------------------------------------------
! Module for MPHTXT file management
! Last update: 04/04/2010
!-----------------------------------------------------------------------
use module_ALLOC
use module_files, only: get_unit
use module_mesh
use module_read_mphtxt
implicit none

!Types
type mphtxt
  private
  character(len=MAXPATH)              :: filename = ' '    !file name
  integer                             :: UNIT     = -1     !associated unit number
end type

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! open: open mphtxt file
!-----------------------------------------------------------------------
subroutine open_mphtxt(this, filename)

type(mphtxt),     intent(inout) :: this     !mphtxt object
character(len=*), intent(in)    :: filename !mphtxt file
integer                         :: ios

! open file
this%filename = filename
this%unit = get_unit()
open (unit=this%unit, file=this%filename, form='formatted', iostat=ios, &
status='old', position='rewind')
if (ios /= 0) call error('mphtxt/open, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! close: close mphtxt file
!-----------------------------------------------------------------------
subroutine close_mphtxt(this)

type(mphtxt),     intent(inout) :: this     !mphtxt object
integer                         :: ios

  ! closes mphtxt file
  close(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt_file/close, #'//trim(string(ios)))

end subroutine

!-----------------------------------------------------------------------
! read: read MPHTXT file
!-----------------------------------------------------------------------
subroutine read_mphtxt(this, m, mphtxt_m, maxdim)

  type(mphtxt),                             intent(inout) :: this     ! mphtxt object
  type(mfm_mesh), dimension(:), allocatable,intent(inout) :: m        ! mfm mesh
  type(pmh_mesh),                           intent(inout) :: mphtxt_m ! pmh_mesh
  integer,                                  intent(inout) :: maxdim   ! dimension detected
  integer :: ios, n, j, i
  logical :: fit(2)

  maxdim = 0

  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt/read/rewind, #'//trim(string(ios)))

  ! Reads the mphtxt file header and allocates the number of pieces
  call read_mphtxt_header(this%unit, mphtxt_m)

  ! Reads every piece of the mesh and calculate the max of its space dimension
  if (.not. allocated(mphtxt_m%pc)) call error('mphtxt/read/object, objects not allocated')
  do i = 1, size(mphtxt_m%pc,1)
      call info('Reading piece '//trim(string(i))//' ...')
      call read_mphtxt_object(this%unit, mphtxt_m%pc(i))
      if (maxdim < mphtxt_m%pc(i)%dim) maxdim = mphtxt_m%pc(i)%dim
  enddo

end subroutine

end module

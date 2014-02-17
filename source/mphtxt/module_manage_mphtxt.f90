module module_manage_mphtxt
!-----------------------------------------------------------------------
! Module for MPHTXT file management
! Last update: 04/04/2010
!-----------------------------------------------------------------------
use module_ALLOC
use module_files, only: get_unit
use module_mesh
use module_read_mphtxt
use module_write_mphtxt
implicit none

!Types
type mphtxt
  private
character(len=MAXPATH) :: filename = ' ' !file name
  integer :: UNIT = -1 !associated unit number
end type

contains

!***********************************************************************
! INPUT PROCEDURES
!***********************************************************************
!-----------------------------------------------------------------------
! open: open mphtxt file
!-----------------------------------------------------------------------
subroutine open_mphtxt(this, filename, st)

  type(mphtxt), intent(inout) :: this !mphtxt object
  character(len=*), intent(in) :: filename !mphtxt file
  integer :: ios
  character(len=*), optional, intent(in) :: st
  character(len=MAXPATH) :: aux
  ! open file
    this%filename = filename
    this%unit = get_unit()
    aux = 'old'
    if(present(st)) then
aux = trim(st)
    endif

open (unit=this%unit, file=this%filename, form='formatted', iostat=ios, &
    status=trim(aux), position='rewind')
    if (ios /= 0) call error('mphtxt/open, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! close: close mphtxt file
!-----------------------------------------------------------------------
subroutine close_mphtxt(this)

type(mphtxt), intent(inout) :: this !mphtxt object
integer :: ios

  ! closes mphtxt file
  close(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt_file/close, #'//trim(string(ios)))

end subroutine

!-----------------------------------------------------------------------
! read: read MPHTXT file
!-----------------------------------------------------------------------
subroutine read_mphtxt(this, m, pmh, maxdim)

  type(mphtxt), intent(inout) :: this ! mphtxt object
  type(mfm_mesh), dimension(:), allocatable,intent(inout) :: m ! mfm mesh
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer, intent(inout) :: maxdim ! dimension detected
  integer :: ios, i

  maxdim = 0

  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt/read/rewind, #'//trim(string(ios)))

  ! Reads the mphtxt file header and allocates the number of pieces
  call read_mphtxt_header(this%unit, pmh)

  ! Reads every piece of the mesh and calculate the max of its space dimension
  if (.not. allocated(pmh%pc)) call error('mphtxt/read/object, objects not allocated')
  do i = 1, size(pmh%pc,1)
      call info('Reading piece '//trim(string(i))//' ...')
      call read_mphtxt_object(this%unit, pmh%pc(i))
      if (maxdim < pmh%pc(i)%dim) maxdim = pmh%pc(i)%dim
  enddo

  ! Build mm in modulef style
  call build_vertices(pmh)

end subroutine

!-----------------------------------------------------------------------
! write: write MPHTXT file
!-----------------------------------------------------------------------
subroutine write_mphtxt(this, pmh)

  type(mphtxt), intent(inout) :: this ! mphtxt object
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer :: i, ios


  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt/read/rewind, #'//trim(string(ios)))

  ! Reads the mphtxt file header and allocates the number of pieces
  call write_mphtxt_header(this%unit, pmh)

  ! Reads every piece of the mesh and calculate the max of its space dimension
  if (.not. allocated(pmh%pc)) call error('mphtxt/read/object, objects not allocated')
  do i = 1, size(pmh%pc,1)
      call info('Writing piece '//trim(string(i))//' ...')
      call write_mphtxt_object(this%unit, pmh%pc(i),i)
  enddo

end subroutine

end module

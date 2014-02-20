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
use module_utils_mphtxt
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
  integer :: ios, i, j
  integer, dimension(size(FEDB,1)) :: minelindx

  maxdim = 0
  minelindx = 1

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
      ! PMH min reference number si 1
      do j = 1, size(pmh%pc(i)%el,1)
        if (minelindx(pmh%pc(i)%el(j)%type) > minval(pmh%pc(i)%el(j)%ref)) then
          minelindx(pmh%pc(i)%el(j)%type) = minval(pmh%pc(i)%el(j)%ref)
        endif
      enddo
  enddo
  do i = 1, size(pmh%pc,1)
      do j = 1, size(pmh%pc(i)%el,1)
        if(minelindx(pmh%pc(i)%el(j)%type) <=0) then
          pmh%pc(i)%el(j)%ref(:) = pmh%pc(i)%el(j)%ref(:) + abs(minelindx(pmh%pc(i)%el(j)%type)) + 1 !PMH min value = 1
        endif
      enddo
  enddo

  ! Build mm in modulef style
  call build_vertices(pmh)

end subroutine

!-----------------------------------------------------------------------
! write: write MPHTXT file
!-----------------------------------------------------------------------
subroutine write_mphtxt(this, pmh)

  type(mphtxt), intent(inout)   :: this ! mphtxt object
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer                       :: i, j, ios, tp, mphlnn
  logical                       :: all_P1
  real(real64), allocatable     :: znod(:,:)


  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt/read/rewind, #'//trim(string(ios)))

  ! Reads the mphtxt file header and allocates the number of pieces
  call write_mphtxt_header(this%unit, pmh)

  ! Reads every piece of the mesh and calculate the max of its space dimension
  if (.not. allocated(pmh%pc)) call error('mphtxt/read/object, objects not allocated')
  do i = 1, size(pmh%pc, 1)
    call build_node_coordinates(pmh%pc(i), i, all_P1, znod) ! build edge midpoints
    do j = 1, size(pmh%pc(i)%el, 1)
      tp = pmh%pc(i)%el(j)%type
      mphlnn = mphtxt_get_lnn(tp)
      if(mphlnn >= FEDB(tp)%lnn + 1) then                   ! build element baricenters
        if(mphlnn == 9) call build_elements_baricenter(pmh%pc(i),i,j,znod)
        if((FEDB(tp)%lnv + FEDB(tp)%lnf) > (FEDB(tp)%lnv))  then    ! build face baricenters
          call build_faces_baricenter(pmh%pc(i),i,j,znod)
        endif
      endif
    enddo
  end do
  do i = 1, size(pmh%pc,1)
      call info('Writing piece '//trim(string(i))//' ...')
      call write_mphtxt_object(this%unit, pmh%pc(i),i,znod)
  enddo

end subroutine

end module

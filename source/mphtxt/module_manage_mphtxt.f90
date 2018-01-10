module module_manage_mphtxt_fcnv

!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! open_mphtxt:  Open a text file associated with a MPHTXT mesh
! close_mphtxt: Close a text file associated with a MPHTXT mesh
! read_mphtxt:  Read a MPHTXT file
! write_mphtxt: Write a MPHTXT file
!-----------------------------------------------------------------------

use basicmod
!use module_mesh
use module_read_mphtxt_fcnv
use module_write_mphtxt_fcnv
use module_utils_mphtxt_fcnv

implicit none

!Types
type mphtxt
  private
  character(len=MAXPATH) :: filename = ' ' !file name
  integer                :: UNIT = -1      !associated unit number
end type

contains


!-----------------------------------------------------------------------
! open_mphtxt(this, filename, st): open mphtxt file
!-----------------------------------------------------------------------
! this:     mphtxt type containing name a unit number
! filename: file name
! st:       status returned in open statement
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
! close_mphtxt(this): close mphtxt file
!-----------------------------------------------------------------------
! this: mphtxt type containing name a unit number
!-----------------------------------------------------------------------

subroutine close_mphtxt(this)

type(mphtxt), intent(inout) :: this !mphtxt object
integer :: ios

  ! closes mphtxt file
  close(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('mphtxt_file/close, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! read_mphtxt(this, pmh, maxdim): read MPHTXT file
!-----------------------------------------------------------------------
! this:   mphtxt type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
! maxdim: max dimension of the PMH pieces
!-----------------------------------------------------------------------
subroutine read_mphtxt(this, pmh, maxdim)
  type(mphtxt), intent(inout) :: this ! mphtxt object
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
  if (.not. allocated(pmh%pc)) call error('mphtxt/read/object, objects not allocated')
  ! Reads every piece of the mesh
  do i = 1, size(pmh%pc,1)
      call info('Reading piece '//trim(string(i))//' ...')
      call read_mphtxt_object(this%unit, pmh%pc(i))
      ! Calculate the max of its space dimension
      if (maxdim < pmh%pc(i)%dim) maxdim = pmh%pc(i)%dim
      ! Calculate the minimum reference number for each FE type, minelindx
      do j = 1, size(pmh%pc(i)%el,1)
        if (minelindx(pmh%pc(i)%el(j)%type) > minval(pmh%pc(i)%el(j)%ref)) then
          minelindx(pmh%pc(i)%el(j)%type) = minval(pmh%pc(i)%el(j)%ref)
        endif
      enddo
  enddo
  ! Set the minimum reference number for each FE type to 1
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
! write_mphtxt(this, pmh): write MPHTXT file
!-----------------------------------------------------------------------
! this:   mphtxt type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
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
    ! build edge midpoints
    call build_node_coordinates(pmh%pc(i), i, all_P1, znod)
    do j = 1, size(pmh%pc(i)%el, 1)
      tp = pmh%pc(i)%el(j)%type
      mphlnn = mphtxt_get_lnn(tp)
      if(mphlnn >= FEDB(tp)%lnn + 1) then
        ! build element baricenters
        call build_elements_baricenter(pmh%pc(i),i,j,znod)
        if((FEDB(tp)%lnv + FEDB(tp)%lnf) > (FEDB(tp)%lnv))  then
          ! build face baricenters
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

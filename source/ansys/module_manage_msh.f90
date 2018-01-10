module module_manage_msh_fcnv

!-----------------------------------------------------------------------
! Module to manage MSH (Ansys) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! open_msh:  Open a text file associated with a MSH mesh
! close_msh: Close a text file associated with a MSH mesh
! read_msh:  Read a MSH file
! write_msh: Write a MSH file
!-----------------------------------------------------------------------

use basicmod
use module_transform_fcnv, only: to_l1
!use module_mesh
use module_pmh_fcnv
use module_read_msh_fcnv
use module_write_msh_fcnv
use module_utils_msh_fcnv


implicit none

!Types
type msh
  private
  character(len=MAXPATH) :: filename = ' ' !file name
  integer                :: UNIT = -1      !associated unit number
end type


contains


!-----------------------------------------------------------------------
! open_msh(this, filename, st): open msh file
!-----------------------------------------------------------------------
! this:     msh type containing name a unit number
! filename: file name
! st:       status returned in open statement
!-----------------------------------------------------------------------

subroutine open_msh(this, filename, st)

  type(msh), intent(inout) :: this !msh object
  character(len=*), intent(in) :: filename !msh file
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
    if (ios /= 0) call error('msh/open, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! close_msh(this): close msh file
!-----------------------------------------------------------------------
! this: msh type containing name a unit number
!-----------------------------------------------------------------------

subroutine close_msh(this)

type(msh), intent(inout) :: this !msh object
integer :: ios

  ! closes msh file
  close(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('msh_file/close, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! read_msh(this, pmh, maxdim): read msh file
!-----------------------------------------------------------------------
! this:   msh type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
! maxdim: max dimension of the PMH pieces
!-----------------------------------------------------------------------

subroutine read_msh(this, pmh)

  type(msh),      intent(inout) :: this ! msh object
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  character(len=MAXPATH)        :: line
  type(msh_faces)               :: faces ! msh faces
  type(msh_cells)               :: cells ! msh faces
  type(msh_zone)                :: izones ! interior zones
  integer                       :: ios, i, indx

  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('msh/read/rewind, #'//trim(string(ios)))

  if(allocated(pmh%pc)) deallocate(pmh%pc)
  allocate(pmh%pc(1)) ! Only one piece

  do
    read(this%unit,fmt='(A)', iostat=ios) line
    if(ios /= 0) exit
    if(is_blank_line(line)) cycle
    indx = get_section_index(line)
    if(indx<0) then
      call read_msh_section(this%unit, line)
!      print*, "ERROR: Unknown section index # "//string(indx)//trim(line)
!      exit
    endif

    select case(indx)

      case(0) !Comments
        call read_msh_comment(this%unit, line)
      case(1) !header
        call read_msh_header(this%unit, line)
      case(2) !Dimension
        call read_msh_dimensions(pmh, line)
      case(10) !Nodes
        call read_msh_nodes(this%unit, pmh, line)
      case(11) !Edges
!        call read_msh_dimensions(pmh, line)
      case(12) !Cells
        call read_msh_cells(this%unit, cells, line)
      case(13) !Faces
        call read_msh_faces(this%unit, faces, line)
      case(39) !Zones & BC
        call read_msh_zones(izones, line)
      case(45) !Zones
        call read_msh_zones(izones, line)
    end select
  enddo

  call check_interface_names(izones)
  call add_faces_to_pmh(faces,izones,pmh)
  call build_cells(faces,cells, pmh)

  do i=1,size(pmh%pc(1)%el,1)
    call reduce(pmh%pc(1)%el(i)%mm, FEDB(pmh%pc(1)%el(i)%type)%lnv ,pmh%pc(1)%el(i)%nel)
    call reduce(pmh%pc(1)%el(i)%ref,pmh%pc(1)%el(i)%nel)
  enddo

  !Build mm in modulef style

  call build_vertices(pmh)

end subroutine


!-----------------------------------------------------------------------
! write_msh(this, pmh): write MSH file
!-----------------------------------------------------------------------
! this:   msh type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine write_msh(this, pmh)

  type(msh), intent(inout)   :: this ! msh object
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer                       :: ios, maxtopdim, maxref
  real(real64), allocatable     :: znod(:,:)


  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('msh/read/rewind, #'//trim(string(ios)))

  maxref = 0
  call write_msh_header(this%unit, pmh, maxtopdim)
  call write_msh_nodes(this%unit, pmh, maxref, z=znod)
  call pmh2msh(this%unit, pmh, maxref, z=znod)
  call write_msh_cells(this%unit, pmh, maxtopdim, maxref)

!  call build_msh_faces(pmh, maxtopdim)

end subroutine


end module

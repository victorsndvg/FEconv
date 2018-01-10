module module_manage_pf3_fcnv

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! open_pf3:  Open a text file associated with a PF3 mesh
! close_pf3: Close a text file associated with a PF3 mesh
! read_pf3:  Read a PF3 file
! write_pf3: Write a PF3 file
!-----------------------------------------------------------------------

use basicmod
!use module_mesh
use module_read_pf3_fcnv
use module_write_pf3_fcnv

implicit none

!Types
type pf3
  private
  character(len=MAXPATH) :: filename = ' ' !file name
  integer                :: UNIT = -1      !associated unit number
end type

contains


!-----------------------------------------------------------------------
! open_pf3(this, filename, st): open pf3 file
!-----------------------------------------------------------------------
! this:     pf3 type containing name a unit number
! filename: file name
! st:       status returned in open statement
!-----------------------------------------------------------------------

subroutine open_pf3(this, filename, st)

  type(pf3), intent(inout) :: this !pf3 object
  character(len=*), intent(in) :: filename !pf3 file
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
    if (ios /= 0) call error('pf3/open, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! close_pf3(this): close PF3 file
!-----------------------------------------------------------------------
! this: pf3 type containing name a unit number
!-----------------------------------------------------------------------

subroutine close_pf3(this)

type(pf3), intent(inout) :: this !pf3 object
integer :: ios

  ! closes PF3 file
  close(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('pf3_file/close, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! read_pf3this, pmh, maxdim): read PF3 file
!-----------------------------------------------------------------------
! this:   pf3 type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
! maxdim: max dimension of the PMH pieces
!-----------------------------------------------------------------------

subroutine read_pf3(this, pmh, maxdim)

  type(pf3), intent(inout) :: this ! PF3 object
  type(pmh_mesh), intent(inout) :: pmh ! pmh_mesh
  integer, intent(inout) :: maxdim ! dimension detected
  character(len=MAXPATH) :: line
  integer :: ios
  integer :: nel, nelvol, nelsur, neledge, nelpoint

  maxdim = 0

  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('pf3/read/rewind, #'//trim(string(ios)))

  ! Reads the pf3 file header and allocates the number of pieces
  call read_pf3_header(this%unit, pmh, nel, nelvol, nelsur, neledge, nelpoint)
  do
    read (unit=this%unit, fmt='(a)', iostat = ios) line
    if(is_iostat_end(ios)) then; exit;
    elseif (ios /= 0) then; call error('pf3/read, #'//trim(string(ios)));
    endif
    if(index(lcase(line),'descripteur de topologie des elements') /= 0) then
      call read_pf3_elements(this%unit, pmh, nel, nelvol, nelsur, neledge, nelpoint)
    elseif(index(lcase(line),'coordonnees des noeuds') /= 0) then
      call read_pf3_coordinates(this%unit, pmh%pc(1))
    elseif(index(lcase(line),'table of the values of') /= 0) then
      call read_pf3_field(this%unit, pmh%pc(1), line)
    endif
   enddo

  call build_vertices(pmh)

end subroutine


!-----------------------------------------------------------------------
! write_pf3(this, pmh): write PF3 file
!-----------------------------------------------------------------------
! this:   pf3 type containing name a unit number
! pmh:    PMH structure storing the piecewise mesh
!-----------------------------------------------------------------------

subroutine write_pf3(this, pmh, infield, outfield, path, param)
  type(pf3),            intent(inout)   :: this !   PF3 object
  type(pmh_mesh),         intent(inout) :: pmh ! pmh_mesh
  character(*), allocatable, intent(in) :: infield(:)  ! In field names
  character(*), allocatable, intent(in) :: outfield(:) ! Out field names
  character(*),              intent(in) :: path !file names
  real(real64), optional,    intent(in) :: param
  integer                       :: i, ios, prevnnod
  logical                       :: all_P1
  real(real64), allocatable     :: znod(:,:)

  ! Write the pf3 file

  rewind(unit=this%unit, iostat=ios)
  if (ios /= 0) call error('pf3/read/rewind, #'//trim(string(ios)))
  call write_pf3_header(this%unit, pmh)
  call write_pf3_elements(this%unit, pmh)
  prevnnod = 0
  write(unit=this%unit, fmt='(a)', iostat = ios) ' COORDONNEES DES NOEUDS'
  if (ios /= 0) call error('module_write_pf3/write_coordinates # write error #'//trim(string(ios)))
  do i=1, size(pmh%pc,1)
    call build_node_coordinates(pmh%pc(i), i, all_P1, znod)
    call write_pf3_coordinates(this%unit, pmh%pc(i), all_P1, znod, prevnnod)
    if(all_P1) then
      prevnnod = prevnnod + pmh%pc(i)%nver
    else
      prevnnod = prevnnod + pmh%pc(i)%nnod
    endif
  enddo
  write(unit=this%unit, fmt='(a)', iostat = ios) ' ==== DECOUPAGE  TERMINE'
  if (ios /= 0) call error('module_write_pf3/write_coordinates # write error #'//trim(string(ios)))
  call write_pf3_node_field(this%unit, pmh, infield, outfield, param)


end subroutine

end module

module module_write_mphtxt_fcnv

!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! write_mphtxt_header: write the header of the MPHTXT file
! write_mphtxt_object: write a pieze of the MPHTXT mesh
! write_mphtxt_etype:  write a element group of the MPHTXT mesh
! write_line:          write a line in the MPHTXT file
! write_comment:       write a comment in the MPHTXT file
! write_empty_line:    write an empty line in the MPHTXT file
! write_string:        write string data in the MPHTXT file
!-----------------------------------------------------------------------

use basicmod
use module_pmh_fcnv
use module_utils_mphtxt_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! write_mphtxt_header(iu, pmh): write mphtxt file header
!-----------------------------------------------------------------------
! iu:  unit number of the MPHTXT file
! pmh: PMH structure for store the mphtxt mesh
!-----------------------------------------------------------------------
subroutine write_mphtxt_header(iu, pmh)
  integer, intent(in)        :: iu  ! File unit number
  type(pmh_mesh), intent(in) :: pmh ! PMH mesh
  integer :: i


  call write_comment(iu,                           '#','Converted with FEconv')
  call write_empty_line(iu)
  call write_comment(iu,                           '#','Major & minor version')
  call write_line(iu, '0 1')
  call write_line(iu, trim(string(size(pmh%pc,1))),'#','number of tags')
  call write_comment(iu,                           '#','Tags')
  do i=1, size(pmh%pc,1)
    call write_string(iu, 'mesh'//trim(string(i)))
  enddo
  ! Types of objects in mphtxt mesh
  call write_line(iu, trim(string(size(pmh%pc,1))),'#','number of types')
  call write_comment(iu,                           '#','Types')
  do i=1, size(pmh%pc,1)
    call write_string(iu, 'obj')
  enddo
  call write_empty_line(iu)

end subroutine


!-----------------------------------------------------------------------
! write_mphtxt_object(iu, pmh_o, n, znod): write a pieze of the MPHTXT mesh
!-----------------------------------------------------------------------
! iu:    unit number of the MPHTXT file
! pmh_o: piece of the PMH structure
! n:     piece number
! znod:  array of node coordinates
!-----------------------------------------------------------------------
subroutine write_mphtxt_object(iu, pmh_o, n, znod)
  integer, intent(in)                      :: iu    ! File unit number
  type(piece), intent(inout)               :: pmh_o ! PMH piece
  integer, intent(in   )                   :: n     ! Piece number
  real(real64), dimension(:,:), intent(in) :: znod
  integer                                  :: i

  call write_comment(iu,                        '#','--------- Object '//trim(string(n))//' ----------')
  call write_empty_line(iu)
  call write_line(iu, '0 0 1')
  call write_string(iu, 'Mesh',                 '#', 'class')
  call write_line(iu, '2',                      '#', 'version')
  call write_line(iu, string(pmh_o%dim),        '#', 'sdim')

  ! Node coords
  if(pmh_o%nnod /= 0 .and. pmh_o%nver /= pmh_o%nnod) then                     ! z = nodes
    call write_line(iu, string(pmh_o%nnod),     '#', 'number of mesh points') ! nnod
    call write_line(iu, '1',                    '#', 'lowest mesh point index')
    call write_empty_line(iu)
    call write_comment(iu,                      '#', 'Mesh point coordinates')

    do i=1,pmh_o%nnod
      call write_line(iu, string(znod(:,i)))
    enddo
  else                                                                        ! z = vertices
    call write_line(iu, string(pmh_o%nver),     '#', 'number of mesh points') ! nver
    call write_line(iu, '1',                    '#', 'lowest mesh point index')
    call write_empty_line(iu)
    call write_comment(iu,                      '#', 'Mesh point coordinates')
    do i=1,pmh_o%nver
      call write_line(iu, string(pmh_o%z(:,i)))
    enddo
  endif

  ! Element groups types
  call write_empty_line(iu)
  call write_line(iu, string(size(pmh_o%el,1)), '#', 'number of element types')
  call write_empty_line(iu)

  do i=1, size(pmh_o%el,1)
    call write_mphtxt_etype(iu, pmh_o%el(i),i)
  enddo


end subroutine


!-----------------------------------------------------------------------
! write_mphtxt_etype(iu, pmh_t,n): write a element group of the MPHTXT mesh
!-----------------------------------------------------------------------
! iu:    unit number of the MPHTXT file
! pmh_t: element group of the PMH structure
! n:     piece number
!-----------------------------------------------------------------------

subroutine write_mphtxt_etype(iu, pmh_t,n)
  integer, intent(in)       :: iu       ! File unit number
  type(elgroup), intent(inout) :: pmh_t ! PMH elgroup
  integer, intent(in)       :: n        ! Piece number
  integer                   :: i

  ! Element type
  call write_comment(iu,                           '#', 'Type #'//trim(string(n)))
  call write_empty_line(iu)
  call write_string(iu,mphtxt_get_desc(pmh_t%type),'#', 'type name')
  call write_empty_line(iu)
  ! Node conectivity
  if(allocated(pmh_t%nn)) then
    call write_line(iu, string(size(pmh_t%nn,1)),  '#', 'number of nodes per element')
    call write_line(iu, string(pmh_t%nel),         '#', 'number of elements')
    call write_comment(iu,                         '#', 'Elements')

    do i=1, size(pmh_t%nn,2)
      call mphtxt_node_ordering(pmh_t%nn(:,i), pmh_t%type)
      call write_line(iu, string(pmh_t%nn(:,i)))
    enddo
    call write_empty_line(iu)
    call write_line(iu, string(size(pmh_t%nn,1)),  '#', 'number of parameter values per element')
  elseif(allocated(pmh_t%mm)) then
    call write_line(iu, string(size(pmh_t%mm,1)),  '#', 'number of nodes per element')
    call write_line(iu, string(pmh_t%nel),         '#', 'number of elements')
    call write_comment(iu,                         '#', 'Elements')
    do i=1, size(pmh_t%mm,2)
      call mphtxt_node_ordering(pmh_t%mm(:,i), pmh_t%type)
      call write_line(iu, string(pmh_t%mm(:,i)))
    enddo
    call write_empty_line(iu)
    call write_line(iu, string(size(pmh_t%nn,1)),  '#', 'number of parameter values per element')

  else
    call error('mphtxt/write_etype# Wrong node definition: No allocated')
  endif
  ! Parameters
  call write_line(iu, '0',                         '#', 'number of parameters')
  call write_comment(iu,                           '#', 'Parameters')
  call write_empty_line(iu)
  ! Geometric entity indices
  call write_line(iu, string(size(pmh_t%ref,1)),   '#', 'number of geometric entity indices')
  call write_comment(iu,                           '#', 'Geometric entity indices')
  do i=1, size(pmh_t%ref,1)
    call write_line(iu, string(pmh_t%ref(i)))
  enddo
  ! Up/down pairs
  call write_empty_line(iu)
  call write_line(iu, '0',                         '#', 'number of up/down pairs')
  call write_comment(iu,                           '#', 'Up/down')
  call write_empty_line(iu)


end subroutine


!-----------------------------------------------------------------------
! write_line(iu,line,ch,comm): write a line in the MPHTXT file
!-----------------------------------------------------------------------
! iu:   unit number of the MPHTXT file
! line: text included in one line
! ch:   comments character (Optional)
! comm: commentary (Optional)
!-----------------------------------------------------------------------

subroutine write_line(iu,line,ch,comm)
  integer, intent(in) :: iu ! File unit number
  character(len=*), intent(in) :: line ! String
  character(len=*), optional, intent(in) :: ch ! String: Comment character
  character(len=*), optional, intent(in) :: comm ! String: Comment
  character(len=MAXPATH) :: aux
  integer :: ios

  if(present(comm)) then
    if(present(ch)) then
      aux = trim(ch)//' '//trim(comm)
    else
      aux = '# '//trim(comm)
    endif
  else
    aux = ''
  endif


  write(unit=iu, fmt='(a)', iostat = ios) trim(line)//' '//trim(aux)
  if (ios /= 0) call error('write_mphtxt/header, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! write_comment(iu,ch,line): write a comment in the MPHTXT file
!-----------------------------------------------------------------------
! iu:   unit number of the MPHTXT file
! ch:   comments character
! line: commentary
!-----------------------------------------------------------------------

subroutine write_comment(iu,ch,line)
  integer, intent(in) :: iu            ! File unit number
  character(len=*), intent(in) :: ch   ! String: Comment character
  character(len=*), intent(in) :: line ! String

  call write_line(iu,trim(ch)//' '//trim(line))

end subroutine


!-----------------------------------------------------------------------
! write_empty_line(iu): write an empty line in the MPHTXT file
!-----------------------------------------------------------------------
! iu:   unit number of the MPHTXT file
!-----------------------------------------------------------------------

subroutine write_empty_line(iu)
  integer, intent(in) :: iu ! File unit number

  call write_line(iu,'')

end subroutine


!-----------------------------------------------------------------------
! write_string(iu,str,ch,comm): write a string in the MPHTXT file
!-----------------------------------------------------------------------
! iu:   unit number of the MPHTXT file
! str:  string to write to the file
! ch:   comments character (Optional)
! line: commentary (Optional)
!-----------------------------------------------------------------------

subroutine write_string(iu,str,ch,comm)
  integer, intent(in) :: iu ! File unit number
  character(len=*), intent(in) :: str ! String
  character(len=*), optional, intent(in) :: ch ! String: Comment character
  character(len=*), optional, intent(in) :: comm ! String: Comment
  character(len=MAXPATH) :: aux1, aux2

  if(present(ch)) then
    aux1 = trim(ch)
  else
    aux1 = ''
  endif

  if(present(comm)) then
    aux2 = trim(comm)
  else
    aux2 = ''
  endif

  call write_line(iu,trim(string(len_trim(str)))//' '//trim(str),trim(aux1),trim(aux2))

end subroutine


end module

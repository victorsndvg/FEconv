module module_write_pf3

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! write_line:          write a line in the PF3 file
! write_comment:       write a comment in the PF3 file
! write_empty_line:    write an empty line in the PF3 file
! write_string:        write string data in the PF3 file
!-----------------------------------------------------------------------

use module_COMPILER_DEPENDANT, only: real64
use module_os_dependant, only: maxpath
use module_report, only:error
use module_convers
use module_mesh
use module_pmh

implicit none


contains


subroutine write_pf3_header(iu, pmh)
  integer, intent(in)        :: iu  ! Unit number for PF3 file
  type(pmh_mesh), intent(in) :: pmh
  integer                    :: i, j, ios
  integer                    :: npoints = 0
  integer                    :: nel = 0
  integer                    :: npointel = 0
  integer                    :: nedgeel = 0
  integer                    :: nsurfel = 0
  integer                    :: nvolel = 0
  integer                    :: nregs = 0
  integer                    :: npointregs = 0
  integer                    :: nedgeregs = 0
  integer                    :: nsurfregs = 0
  integer                    :: nvolregs = 0

    if(.not. allocated(pmh%pc)) call error('module_write_pf3/write_header # Piece not allocated')


    ! Count the number of elements and regions of each topological dimension
    do i = 1, size(pmh%pc,1)
      if(.not. allocated(pmh%pc(i)%el)) call error('module_write_pf3/write_header # Element group not allocated')
      do j = 1, size(pmh%pc(i)%el,1)
        nel = nel + pmh%pc(i)%el(j)%nel
        nregs = nregs + 1
        if(FEDB(pmh%pc(i)%el(j)%type)%tdim == 0) then
          npointel = npointel + pmh%pc(i)%el(j)%nel
          npointregs = npointregs + 1
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 1) then
          nedgeel = nedgeel + pmh%pc(i)%el(j)%nel
          nedgeregs = nedgeregs + 1
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 2) then
          nsurfel = nsurfel + pmh%pc(i)%el(j)%nel
          nsurfregs = nsurfregs + 1
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 3) then
          nvolel = nvolel + pmh%pc(i)%el(j)%nel
          nvolregs = npointregs + 1
        endif
     enddo
      if(pmh%pc(i)%nnod /= 0) then 
        npoints = npoints + pmh%pc(i)%nnod
      else 
        npoints = npoints + pmh%pc(i)%nver
      endif
    enddo

    ! Write Flux PF3 header
    write(unit=iu, fmt='(a)', iostat = ios)                     'File converted with FEconv'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(pmh%pc(1)%dim), 'NOMBRE DE DIMENSIONS DU DECOUPAGE'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nel),        'NOMBRE  D''ELEMENTS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nvolel),     'NOMBRE  D''ELEMENTS VOLUMIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nsurfel),    'NOMBRE  D''ELEMENTS SURFACIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nedgeel),    'NOMBRE  D''ELEMENTS LINEIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(npointel),   'NOMBRE  D''ELEMENTS PONCTUELS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(0),          'NOMBRE DE MACRO-ELEMENTS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(npoints),    'NOMBRE DE POINTS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nregs),       'NOMBRE DE REGIONS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nvolregs),    'NOMBRE DE REGIONS VOLUMIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nsurfregs),   'NOMBRE DE REGIONS SURFACIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(nedgeregs),   'NOMBRE DE REGIONS LINEIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(npointregs),  'NOMBRE DE REGIONS PONCTUELLES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(0),          'NOMBRE DE REGIONS MACRO-ELEMENTAIRES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(20),         'NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    write(unit=iu, fmt='(a)', iostat = ios) string(20),         'NOMBRE DE POINTS D''INTEGRATION / ELEMENT (MAX)'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))

    ! Write Flux PF3 regions
    write(unit=iu, fmt='(a)', iostat = ios) 'NOMS DES REGIONS'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))

    ! Write Flux PF3 volume region names
    if(nvolregs /= 0) write(unit=iu, fmt='(a)', iostat = ios) 'REGIONS VOLUMIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    do i = 1, nvolregs
      write(unit=iu, fmt='(a)', iostat = ios) 'REGIONVOLUME_'//string(i)
      if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    enddo

    ! Write Flux PF3 surface region names
    if(nsurfregs /= 0) write(unit=iu, fmt='(a)', iostat = ios) 'REGIONS SURFACIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    do i = 1, nsurfregs
      write(unit=iu, fmt='(a)', iostat = ios) 'REGIONFACE_'//string(i)
      if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    enddo

    ! Write Flux PF3 edge region names
    if(nedgeregs /= 0) write(unit=iu, fmt='(a)', iostat = ios) 'REGIONS LINEIQUES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    do i = 1, nedgeregs
      write(unit=iu, fmt='(a)', iostat = ios) 'REGIONLINE_'//string(i)
      if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    enddo

    ! Write Flux PF3 point region names
    if(npointregs /= 0) write(unit=iu, fmt='(a)', iostat = ios) 'REGIONS PONCTUELLES'
    if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    do i = 1, npointregs
      write(unit=iu, fmt='(a)', iostat = ios) 'REGIONPOINT_'//string(i)
      if (ios /= 0) call error('module_write_pf3/write_header # write error #'//trim(string(ios)))
    enddo


    ! Write Flux PF3 region descriptors
    do i = 1, nvolregs; 
      write(unit=iu, fmt='(a)', iostat = ios) string(0)//' '//string(4)//' '//string(5)
    enddo
    do i = 1, nsurfregs
      write(unit=iu, fmt='(a)', iostat = ios) string(0)//' '//string(3)//' '//string(5)
    enddo
    do i = 1, nedgeregs
      write(unit=iu, fmt='(a)', iostat = ios) string(0)//' '//string(2)//' '//string(5)
    enddo
    do i = 1, npointregs
      write(unit=iu, fmt='(a)', iostat = ios) string(0)//' '//string(1)//' '//string(5)
    enddo

end subroutine

subroutine write_pf3_elements(iu, pmh)
  integer, intent(in)        :: iu  ! Unit number for PF3 file
  type(pmh_mesh), intent(in) :: pmh
  integer                    :: i, j, k, ios, desc1, desc2, desc3, ref, tdim, lnn

    ! Write Flux PF3 header
    write(unit=iu, fmt='(a)', iostat = ios) 'DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS'
    if (ios /= 0) call error('module_write_pf3/write_elements # write error #'//trim(string(ios)))

    do i = 1, size(pmh%pc,1)
      do j = 1, size(pmh%pc(i)%el,1)
          desc1 = pf3_get_element_desc1(pmh%pc(i)%el(j)%type)
          desc2 = pf3_get_element_desc2(pmh%pc(i)%el(j)%type)
          desc3 = pf3_get_element_desc3(pmh%pc(i)%el(j)%type)
          tdim = FEDB(pmh%pc(i)%el(j)%type)%tdim
          if(allocated(pmh%pc(i)%el(j)%nn)) then
            lnn = size(pmh%pc(i)%el(j)%nn,1)
          else
            lnn = size(pmh%pc(i)%el(j)%mm,1)
          endif
        do k = 1, pmh%pc(i)%el(j)%nel
          ref = pmh%pc(i)%el(j)%ref(k)
          write(unit=iu, fmt='(12I10)', iostat = ios) j,desc1,desc2,ref,tdim,0,desc3,lnn,0,0,0,0
          if(allocated(pmh%pc(i)%el(j)%nn)) then
            write(unit=iu, fmt='(a)', iostat = ios) pmh%pc(i)%el(j)%nn(:,k)
            if (ios /= 0) call error('module_write_pf3/write_elements # write error #'//trim(string(ios)))
          else
            write(unit=iu, fmt='(a)', iostat = ios) pmh%pc(i)%el(j)%mm(:,k)
            if (ios /= 0) call error('module_write_pf3/write_elements # write error #'//trim(string(ios)))
          endif
        enddo
      enddo
    enddo

end subroutine

subroutine write_pf3_coordinates(iu, pmh)
  integer, intent(in)           :: iu  ! Unit number for PF3 file
  type(pmh_mesh), intent(inout) :: pmh
end subroutine


function pf3_get_element_desc1(tp) result(res)
  integer, intent(in)  :: tp
  integer              :: res

  res = 0

  if(tp ==  check_fe(.true., 1, 1, 0, 0)) then ! Vertex
    res = 1
  elseif(tp == check_fe(.true., 2, 2, 1, 0) .or. tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge P1-P2  
    res = 2
  elseif(tp == check_fe(.true., 3, 3, 3, 0) .or. tp == check_fe(.false., 6, 3, 3, 0)) then ! Triangle P1-P2
    res = 3
  elseif(tp == check_fe(.true., 4, 4, 6, 4) .or. tp == check_fe(.true., 10, 4, 6, 4)) then ! Tetrahedron  P1-P2
    res = 5
  endif

end function


function pf3_get_element_desc2(tp) result(res)
  integer, intent(in)  :: tp
  integer              :: res

  res = 0

  if(tp ==  check_fe(.true., 1, 1, 0, 0)) then ! Vertex
    res = 1
  elseif(tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge P1
    res = 2
  elseif(tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge P2  
    res = 3
  elseif(tp == check_fe(.true., 3, 3, 3, 0) .or. tp == check_fe(.false., 6, 3, 3, 0)) then  ! Triangle P1-P2  
    res = 7
  elseif(tp == check_fe(.true., 4, 4, 6, 4)) then  ! Tetrahedron P1
    res = 4
  elseif(tp == check_fe(.true., 10, 4, 6, 4)) then ! Tetrahedron P2
    res = 15
  endif

end function



function pf3_get_element_desc3(tp) result(res)
  integer, intent(in)  :: tp
  integer              :: res

  res = 0

  if(tp ==  check_fe(.true., 1, 1, 0, 0)) then ! Vertex
    res = 2
  elseif(tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge P1
    res = 3
  elseif(tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge P2  
    res = 4
  elseif(tp == check_fe(.true., 3, 3, 3, 0)) then 
    res = 5
  elseif(tp == check_fe(.false., 6, 3, 3, 0)) then  ! Triangle P1-P2  
    res = 6
  elseif(tp == check_fe(.true., 4, 4, 6, 4)) then  ! Tetrahedron P1
    res = 10
  elseif(tp == check_fe(.true., 10, 4, 6, 4)) then ! Tetrahedron P2
    res = 11
  endif

end function


!-----------------------------------------------------------------------
! write_line(iu,line,ch,comm): write a line in the PF3 file
!-----------------------------------------------------------------------
! iu:   unit number of the PF3 file
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
  if (ios /= 0) call error('write_pf3/header, #'//trim(string(ios)))

end subroutine


!-----------------------------------------------------------------------
! write_comment(iu,ch,line): write a comment in the PF3 file
!-----------------------------------------------------------------------
! iu:   unit number of the PF3 file
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
! write_empty_line(iu): write an empty line in the PF3 file
!-----------------------------------------------------------------------
! iu:   unit number of the PF3 file
!-----------------------------------------------------------------------

subroutine write_empty_line(iu)
  integer, intent(in) :: iu ! File unit number

  call write_line(iu,'')

end subroutine


!-----------------------------------------------------------------------
! write_string(iu,str,ch,comm): write a string in the PF3 file
!-----------------------------------------------------------------------
! iu:   unit number of the PF3 file
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

module module_write_pf3_fcnv

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! write_pf3_header: write PF3 header
! write_pf3_elements: write PF3 elements
! write_pf3_coordinates: write PF3 coordinates
! write_pf3_node_field: write PF3 node field
!-----------------------------------------------------------------------

use basicmod
!use module_mesh
use module_pmh_fcnv
use module_utils_pf3_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! write_pf3_header(iu, pmh): write PF3 header
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! pmh: pmh mesh
!-----------------------------------------------------------------------

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
        if(FEDB(pmh%pc(i)%el(j)%type)%tdim == 0) then
          npointel = npointel + pmh%pc(i)%el(j)%nel
          npointregs = npointregs + size(unique(pmh%pc(i)%el(j)%ref),1)
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 1) then
          nedgeel = nedgeel + pmh%pc(i)%el(j)%nel
          nedgeregs = nedgeregs + size(unique(pmh%pc(i)%el(j)%ref),1)
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 2) then
          nsurfel = nsurfel + pmh%pc(i)%el(j)%nel
          nsurfregs = nsurfregs + size(unique(pmh%pc(i)%el(j)%ref),1)
        elseif(FEDB(pmh%pc(i)%el(j)%type)%tdim == 3) then
          nvolel = nvolel + pmh%pc(i)%el(j)%nel
          nvolregs = nvolregs + size(unique(pmh%pc(i)%el(j)%ref),1)
        endif
      enddo
      nregs = npointregs + nedgeregs + nsurfregs + nvolregs
      if(pmh%pc(i)%nnod /= 0) then
        npoints = npoints + pmh%pc(i)%nnod
      else
        npoints = npoints + pmh%pc(i)%nver
      endif
    enddo

    ! Write Flux PF3 header
    call write_header_line(iu,'File converted with FEconv')
    call write_header_line(iu,'NOMBRE DE DIMENSIONS DU DECOUPAGE',pmh%pc(1)%dim)
    call write_header_line(iu,'NOMBRE  D''ELEMENTS', nel)
    call write_header_line(iu,'NOMBRE  D''ELEMENTS VOLUMIQUES', nvolel)
    call write_header_line(iu,'NOMBRE  D''ELEMENTS SURFACIQUES', nsurfel)
    call write_header_line(iu,'NOMBRE  D''ELEMENTS LINEIQUES', nedgeel)
    call write_header_line(iu,'NOMBRE  D''ELEMENTS PONCTUELS', npointel)
    call write_header_line(iu,'NOMBRE DE MACRO-ELEMENTS', 0)
    call write_header_line(iu,'NOMBRE DE POINTS', npoints)
    call write_header_line(iu,'NOMBRE DE REGIONS', nregs)
    call write_header_line(iu,'NOMBRE DE REGIONS VOLUMIQUES', nvolregs)
    call write_header_line(iu,'NOMBRE DE REGIONS SURFACIQUES', nsurfregs)
    call write_header_line(iu,'NOMBRE DE REGIONS LINEIQUES', nedgeregs)
    call write_header_line(iu,'NOMBRE DE REGIONS PONCTUELLES', npointregs)
    call write_header_line(iu,'NOMBRE DE REGIONS MACRO-ELEMENTAIRES', 0)
    call write_header_line(iu,'NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)', 20)
    call write_header_line(iu,'NOMBRE DE POINTS D''INTEGRATION / ELEMENT (MAX)', 20)

    ! Write Flux PF3 regions
    call write_header_line(iu,'NOMS DES REGIONS')

    ! Write Flux PF3 volume region names
    if(nvolregs>0) call write_header_line(iu, '    REGIONS VOLUMIQUES')
    do i = 1, nvolregs
      call write_header_line(iu,'REGIONVOLUME_'//string(i))
    enddo

    ! Write Flux PF3 surface region names
    if(nsurfregs>0) call write_header_line(iu, '    REGIONS SURFACIQUES')
    do i = 1, nsurfregs
      call write_header_line(iu, 'REGIONFACE_'//string(i))
    enddo

    ! Write Flux PF3 edge region names
    if(nedgeregs>0) call write_header_line(iu, '    REGIONS LINEIQUES')
    do i = 1, nedgeregs
      call write_header_line(iu, 'REGIONLINE_'//string(i))
    enddo

    ! Write Flux PF3 point region names
    if(npointregs>0) call write_header_line(iu, '    REGIONS PONCTUELLES')
    do i = 1, npointregs
      call write_header_line(iu, 'REGIONPOINT_'//string(i))
    enddo


    ! Write Flux PF3 region descriptors
    do i = 1, nvolregs;
      write(unit=iu, fmt='(a,3I6)', iostat = ios)  ' ',0,4,5
      if (ios /= 0) call error('module_utils_pf3/write_header # '//trim(string(ios)))
    enddo
    do i = 1, nsurfregs
      write(unit=iu, fmt='(a,3I6)', iostat = ios)  ' ',0,3,5
      if (ios /= 0) call error('module_utils_pf3/write_header # '//trim(string(ios)))
    enddo
    do i = 1, nedgeregs
      write(unit=iu, fmt='(a,3I6)', iostat = ios)  ' ',0,2,5
      if (ios /= 0) call error('module_utils_pf3/write_header # '//trim(string(ios)))
    enddo
    do i = 1, npointregs
      write(unit=iu, fmt='(a,3I6)', iostat = ios)  ' ',0,1,5
      if (ios /= 0) call error('module_utils_pf3/write_header # '//trim(string(ios)))
    enddo

end subroutine


!-----------------------------------------------------------------------
! write_pf3_elements(iu, pmh): write PF3 elements
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! pmh: pmh mesh
!-----------------------------------------------------------------------

subroutine write_pf3_elements(iu, pmh)
  integer, intent(in)        :: iu  ! Unit number for PF3 file
  type(pmh_mesh), intent(inout) :: pmh
  integer                    :: i, j, k, ios, desc1, desc2, desc3, ref, tdim, lnn, prevnnod, prevnel
  character(len=MAXPATH)     :: fm

    ! Write Flux PF3 header
    write(unit=iu, fmt='(a)', iostat = ios) ' DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS'
    if (ios /= 0) call error('module_utils_pf3/write_elements # write error #'//trim(string(ios)))

    prevnnod = 0
    prevnel = 0

    do i = 1, size(pmh%pc,1)
      do j = 1, size(pmh%pc(i)%el,1)
          desc1 = pf3_get_element_desc1(pmh%pc(i)%el(j)%type)
          desc2 = pf3_get_element_desc2(pmh%pc(i)%el(j)%type)
          desc3 = pf3_get_element_desc3(pmh%pc(i)%el(j)%type)
          tdim = FEDB(pmh%pc(i)%el(j)%type)%tdim
          lnn = FEDB(pmh%pc(i)%el(j)%type)%lnn
        do k = 1, pmh%pc(i)%el(j)%nel
          ref = pmh%pc(i)%el(j)%ref(k)
          write(unit=iu, fmt='(a,12I9)', iostat = ios) '      ',k+prevnel,desc1,desc2,ref,tdim+1,0,desc3,lnn,0,0,0,0
          if(.not. FEDB(pmh%pc(i)%el(j)%type)%nver_eq_nnod .and. allocated(pmh%pc(i)%el(j)%nn)) then
            write(unit=fm, fmt='(a)', iostat = ios)  '(a,'//trim(string(FEDB(pmh%pc(i)%el(j)%type)%lnn))//'I9)'
            write(unit=iu, fmt=fm, iostat = ios) '      ', &
              pmh2pf3_ordering(pmh%pc(i)%el(j)%nn(:,k), pmh%pc(i)%el(j)%type, prevnnod)
            if (ios /= 0) call error('module_write_pf3/write_elements # write error #'//trim(string(ios)))
          elseif(FEDB(pmh%pc(i)%el(j)%type)%nver_eq_nnod .and. allocated(pmh%pc(i)%el(j)%mm)) then
            write(unit=fm, fmt='(a)', iostat = ios)  '(a,'//trim(string(FEDB(pmh%pc(i)%el(j)%type)%lnn))//'I9)'
            write(unit=iu, fmt=fm, iostat = ios) '      ', &
              pmh2pf3_ordering(pmh%pc(i)%el(j)%mm(:,k), pmh%pc(i)%el(j)%type, prevnnod)
            if (ios /= 0) call error('module_write_pf3/write_elements # write error #'//trim(string(ios)))
          else
            call error('module_write_pf3/write_elements # connectivity array not allocated')
          endif
        enddo
        prevnel = prevnel + pmh%pc(i)%el(j)%nel
      enddo
      if(pmh%pc(i)%nnod /= 0) then
        prevnnod = prevnnod + pmh%pc(i)%nnod
      else
        prevnnod = prevnnod + pmh%pc(i)%nver
      endif
    enddo

end subroutine


!-----------------------------------------------------------------------
! write_pf3_coordinates(iu, pc, all_P1, znod, prevnnod): write PF3 coordinates
!-----------------------------------------------------------------------
! iu:       unit number of the PF3 file
! pc:       piece structure of the mesh
! all_P1:   logical flag, true if all element groups are Lagrange P1
! znod:     array of node coordinates, necessary if not all element groups are Lagrange P1
! prevnnod: number of nodes writed in previous pieces
!-----------------------------------------------------------------------

subroutine write_pf3_coordinates(iu, pc, all_P1, znod, prevnnod)
  integer, intent(in)           :: iu  ! Unit number for PF3 file
  type(piece), intent(inout)    :: pc
  integer, intent(in)                   ::prevnnod
  logical, intent(in)                   :: all_P1
  real(real64), allocatable, intent(in) :: znod(:,:)
  integer                       :: j, ios

!    do i = 1, size(pmh%pc,1)

      if(all_P1) then
        do j = 1, pc%nver
          write(unit=iu, fmt='(I8,a)', iostat = ios) j+prevnnod, trim(string(pc%z(:,j)))
          if (ios /= 0) call error('module_write_pf3/write_coordinates # write error #'//trim(string(ios)))
        enddo
      else
        do j = 1, size(znod,2)
          write(unit=iu, fmt='(I8,a)', iostat = ios) j+prevnnod, trim(string(znod(:,j)))
          if (ios /= 0) call error('module_write_pf3/write_coordinates # write error #'//trim(string(ios)))
        enddo
      endif

!    enddo

end subroutine


!-----------------------------------------------------------------------
! write_pf3_node_field(iu, pmh): write PF3 node field
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! pmh: pmh mesh
!-----------------------------------------------------------------------

subroutine write_pf3_node_field(iu, pmh, infield, outfield, param)
  integer,                   intent(in) :: iu  ! Unit number for PF3 file
  type(pmh_mesh),         intent(inout) :: pmh
  character(*), allocatable, intent(in) :: infield(:)  ! In field names
  character(*), allocatable, intent(in) :: outfield(:) ! Out field names
  real(real64), optional,    intent(in) :: param
  character(len=maxpath)                :: fieldname
  integer                               :: i, j, k, ios
  integer                               :: tnnod, ncomp, pindex
  logical                               :: singlefield, exists

  singlefield = .false.
  exists = .false.

!  if(.not. allocated(infield)) return

  if(allocated(infield)) then
    if(size(infield,1) /= 1) then
      call info('Only one field can be saved in PF3 file.')
      return
    endif
    elseif(allocated(outfield)) then
      if(size(outfield,1) /= 1) then
        call info('Only one field can be saved in PF3 file.')
        return
      endif
  endif

  singlefield = .true.
  tnnod = 0

  do i=1, size(pmh%pc,1)
    if(.not. allocated(infield) .and. get_piece_num_fields(pmh%pc(i), 'node') /= 1) return
    if(allocated(pmh%pc(i)%fi)) then
      do j=1, size(pmh%pc(i)%fi,1)

        if(.not. allocated(infield) .and. get_piece_num_fields(pmh%pc(i), 'node') /= 1) return
        if(allocated(infield)) then
          if(trim(infield(1)) /= trim(pmh%pc(i)%fi(j)%name)) cycle
        endif
        if(.not. allocated(pmh%pc(i)%fi(j)%val)) call error('write/pf3/node/field # Not allocated!')
        if(present(param)) then
          do k=1,size(pmh%pc(i)%fi(j)%param)
            if(pmh%pc(i)%fi(j)%param(k)-param<pmh%ztol) pindex = k
          enddo
        else
          pindex = 1
        endif
        if(.not. allocated(outfield)) then
          fieldname = trim(pmh%pc(i)%fi(j)%name)
        else
          fieldname = trim(outfield(1))
        endif
        call replace(fieldname, ' ', '_')
        tnnod = tnnod + pmh%pc(i)%nnod
        ncomp = size(pmh%pc(i)%fi(j)%val,1)
        exists = .true.
      enddo
    endif
  enddo


  if(.not. exists) return

!  if(trim(outfield(1)) /= '*') fieldname = trim(outfield))
  call info('Writing node field: '//trim(fieldname))

  write(unit=iu, fmt=*, iostat = ios) ''
  write(unit=iu, fmt='(a,a50,a)', iostat = ios) &
     & 'Table of the values of ',fieldname,' (Local values on the nodes all over the mesh)'
  write(unit=iu, fmt=*, iostat = ios) ncomp, tnnod, ' (Nb of components, nb of points)'
  write(unit=iu, fmt=*, iostat = ios) ''
  if (ios /= 0) call error('write/pf3/node/field # write error #'//trim(string(ios)))

  do i=1, size(pmh%pc,1)
    if(allocated(pmh%pc(i)%fi)) then
      do j=1, size(pmh%pc(i)%fi,1)
        if(.not. allocated(infield) .and. get_piece_num_fields(pmh%pc(i), 'node') /= 1) return
        if(allocated(infield)) then
          if(trim(infield(1)) /= trim(pmh%pc(i)%fi(j)%name)) cycle
        endif
        if(.not. allocated(pmh%pc(i)%fi(j)%val)) call error('write/pf3/node/field # Not allocated!')
        do k=1, size(pmh%pc(i)%fi(j)%val,2)
          write(unit=iu, fmt='(a)', iostat = ios) trim(string(pmh%pc(i)%fi(j)%val(:,k,pindex)))
          if (ios /= 0) call error('write/pf3/node/field # write error #'//trim(string(ios)))
        enddo
      enddo
    endif
  enddo


end subroutine


end module

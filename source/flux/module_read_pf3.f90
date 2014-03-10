module module_read_pf3

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! read_pf3_header(iu): read PF3 file header
!-----------------------------------------------------------------------

use module_COMPILER_DEPENDANT, only: real64
use module_os_dependant, only: maxpath
use module_report, only:error
use module_convers
use module_ALLOC
use module_mesh
use module_pmh

implicit none


contains


!-----------------------------------------------------------------------
! read_pf3_header(iu): read PF3 file header
!-----------------------------------------------------------------------
! iu: unit number of the PF3 file
!-----------------------------------------------------------------------

subroutine read_pf3_header(iu, pmh, nel, nelvol, nelsur, neledge, nelpoint)

  integer, intent(in)           :: iu ! Unit number for PF3 file
  type(pmh_mesh), intent(inout) :: pmh
  integer, intent(inout)        :: nel,nelvol,nelsur,neledge,nelpoint
  character(len=MAXPATH)        :: line
  integer                       :: i, ngroups, ios
  logical                       :: readed = .false.

    nel=0;nelvol=0;nelsur=0;neledge=0;nelpoint=0;ngroups=0



    if(allocated(pmh%pc)) deallocate(pmh%pc)
    allocate(pmh%pc(1)) ! 1 piece mesh

    do i=1, 20
      read (unit=iu, fmt='(a)', iostat = ios) line
      if (ios /= 0) call error('read_pf3/header, #'//trim(string(ios)))

      if(index(line,'NOMBRE DE DIMENSIONS DU DECOUPAGE') /= 0) then
        read(line,*) pmh%pc(1)%dim
      elseif(index(line,'NOMBRE  D''ELEMENTS VOLUMIQUES') /= 0) then
        read(line,*) nelvol
        if(nelvol /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS SURFACIQUES') /= 0) then
        read(line,*) nelsur
        if(nelsur /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS LINEIQUES') /= 0) then
        read(line,*) neledge
        if(neledge /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS PONCTUELS') /= 0) then
        read(line,*) nelpoint
        if(nelpoint /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS') /= 0) then
        read(line,*) nel
      elseif(index(line,'NOMBRE DE MACRO-ELEMENTS') /= 0) then
!        read(line,*) nelme
!        if(nelme /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE DE REGIONS VOLUMIQUES') /= 0) then
!        read(line,*) nregvol
      elseif(index(line,'NOMBRE DE REGIONS SURFACIQUES') /= 0) then
!        read(line,*) nregsur
      elseif(index(line,'NOMBRE DE REGIONS LINEIQUES') /= 0) then
!        read(line,*) nregedge
      elseif(index(line,'NOMBRE DE REGIONS PONCTUELLES') /= 0) then
!        read(line,*) nregpoint
      elseif(index(line,'NOMBRE DE REGIONS MACRO-ELEMENTAIRES') /= 0) then
!        read(line,*) nregvol
      elseif(index(line,'NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)') /= 0) then
!        read(line,*) 
      elseif(index(line,'NOMBRE DE POINTS D''INTEGRATION / ELEMENT (MAX)') /= 0) then
!        read(line,*) 
        readed = .true.  ! Force exit. Solo si la cabecera está ordenada
      elseif(index(line,'NOMBRE DE POINTS') /= 0) then
        read(line,*) pmh%pc(1)%nnod
      elseif(index(line,'NOMBRE DE REGIONS') /= 0) then
!        read(line,*) nreg
      endif

      if(readed) exit

    enddo

    if(allocated(pmh%pc(1)%el)) deallocate(pmh%pc(1)%el)
    allocate(pmh%pc(1)%el(ngroups)) ! 1 piece mesh

    if(allocated(pmh%pc(1)%z)) deallocate(pmh%pc(1)%z)
    allocate(pmh%pc(1)%z(pmh%pc(1)%dim,pmh%pc(1)%nnod)) ! 1 piece mesh


end subroutine


!-----------------------------------------------------------------------
! read_pf3_elements(iu, nel): read PF3 elements
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! nel: number of elements
!-----------------------------------------------------------------------

subroutine read_pf3_elements(iu, pmh, nel, nelvol, nelsur, neledge, nelpoint)

  integer, intent(in)           :: iu  ! Unit number for PF3 file
  type(pmh_mesh), intent(inout) :: pmh
  integer, intent(in)           :: nel, nelvol, nelsur, neledge, nelpoint
  character(len=MAXPATH)        :: line
  integer                       :: nodnum, lnn, ref, ngroup, ttype, aux
  integer                       :: desc1, desc2, desc3 ! element type descriptors
  integer                       :: i, k, ios
  integer,dimension(4,2)        :: tgroup = 0 ! 4:3D(vol),3:2D(sur),2:1D(edge),0:0D(point)


    ! Associates each group to an element type
    ngroup = 0
    if(nelvol /= 0) then
      ngroup = ngroup + 1
      pmh%pc(1)%el(ngroup)%nel = nelvol
      tgroup(4,1) = ngroup
    endif
    if(nelsur /= 0) then
      ngroup = ngroup + 1
      pmh%pc(1)%el(ngroup)%nel = nelsur
      tgroup(3,1) = ngroup
    endif
    if(neledge /= 0) then
      ngroup = ngroup + 1
      pmh%pc(1)%el(ngroup)%nel = neledge
      tgroup(2,1) = ngroup
    endif
    if(nelpoint /= 0) then
      ngroup = ngroup + 1
      pmh%pc(1)%el(ngroup)%nel = nelpoint
      tgroup(1,1) = ngroup
    endif

    do i = 1, nel
      ! Read element description
      read (unit=iu, fmt=*, iostat = ios) nodnum,desc1,desc2,ref,ttype,aux,desc3,lnn
      if (ios /= 0) call error('pf3/read/elements, #'//trim(string(ios)))

      if(.not. allocated(pmh%pc(1)%el(tgroup(ttype,1))%nn)) then
        allocate(pmh%pc(1)%el(tgroup(ttype,1))%nn(lnn,pmh%pc(1)%el(tgroup(ttype,1))%nel))
        allocate(pmh%pc(1)%el(tgroup(ttype,1))%ref(pmh%pc(1)%el(tgroup(ttype,1))%nel))
        pmh%pc(1)%el(tgroup(ttype,1))%type = pf3_assign_element_type(desc1, desc2, desc3)
      endif
      ! Counts number of elements readed for each element group
      tgroup(ttype,2) = tgroup(ttype,2) + 1

      ! Assign reference number
      pmh%pc(1)%el(tgroup(ttype,1))%ref(tgroup(ttype,2)) = ref

      ! Read list of nodes
      read (unit=iu, fmt=*, iostat = ios) (pmh%pc(1)%el(tgroup(ttype,1))%nn(k,tgroup(ttype,2)), k=1,lnn)
      if (ios /= 0) call error('pf3/read/elements, #'//trim(string(ios)))

    enddo

end subroutine

!-----------------------------------------------------------------------
! read_pf3_coordinates(iu, nnodk, dim): read PF3 elements
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! pc:  piece
!-----------------------------------------------------------------------

subroutine read_pf3_coordinates(iu, pc)

  integer, intent(in) :: iu   ! Unit number for PF3 file
  type(piece), intent(inout) :: pc
  character(len=MAXPATH) :: line
  integer :: indx
  integer :: i, j, ios

  if(allocated(pc%z)) deallocate(pc%z)
  allocate(pc%z(pc%dim,pc%nnod))

  do i = 1, pc%nnod
    read (unit=iu, fmt=*, iostat = ios)  indx, (pc%z(j,i), j=1,pc%dim)
    if (ios /= 0) call error('pf3/read/coordinates, #'//trim(string(ios)))
  enddo

end subroutine



function pf3_assign_element_type(desc1, desc2, desc3) result(res)
  integer, intent(in)  :: desc1, desc2, desc3
  integer              :: res

  if(desc1 == 1 .and. desc2 == 1 .and. desc3 == 2) then ! Vertex
    res = check_fe(.true., 1, 1, 0, 0)
    call info('Element type: Vertex')
  elseif(desc1 == 2 .and. desc2 == 2 .and. desc3 == 3) then ! Edge P1  
    res = check_fe(.true., 2, 2, 1, 0)
    call info('Element type: Edge Lagrange P1')
  elseif(desc1 == 2 .and. desc2 == 3 .and. desc3 == 4) then ! Edge P2
    res = check_fe(.false., 3, 2, 1, 0)
    call info('Element type: Edge Lagrange P2')
  elseif(desc1 == 3 .and. desc2 == 7 .and. desc3 == 5) then ! Triangle P1
    res = check_fe(.true., 3, 3, 3, 0)
    call info('Element type: Triangle Lagrange P1')
  elseif(desc1 == 3 .and. desc2 == 7 .and. desc3 == 6) then ! Triangle P2
    res = check_fe(.false., 6, 3, 3, 0)
    call info('Element type: Triangle Lagrange P2')
  elseif(desc1 == 5 .and. desc2 == 4 .and. desc3 == 10) then ! Tetrahedron  P1
    res = check_fe(.true., 4, 4, 6, 4)
    call info('Element type: Tetrahedron Lagrange P1')
  elseif(desc1 == 5 .and. desc2 == 15 .and. desc3 == 11) then ! Tetrahedron  P2
    res = check_fe(.false., 10, 4, 6, 4)
    call info('Element type: Tetrahedron Lagrange P2')
  else
    call info('Element type: Unknown element type')
    res = 0
  endif

end function

end module

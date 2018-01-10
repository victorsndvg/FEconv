module module_read_pf3_fcnv

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/03/2014
!
! PUBLIC PROCEDURES:
! read_pf3_header: read PF3 file header
! read_pf3_elements: read PF3 elements
! read_pf3_coordinates: read PF3 coordinates
!-----------------------------------------------------------------------

use basicmod
!use module_mesh
use module_pmh_fcnv
use module_utils_pf3_fcnv

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
!        if(nelvol /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS SURFACIQUES') /= 0) then
        read(line,*) nelsur
!        if(nelsur /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS LINEIQUES') /= 0) then
        read(line,*) neledge
!        if(neledge /= 0) ngroups=ngroups+1
      elseif(index(line,'NOMBRE  D''ELEMENTS PONCTUELS') /= 0) then
        read(line,*) nelpoint
!        if(nelpoint /= 0) ngroups=ngroups+1
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

!    if(allocated(pmh%pc(1)%el)) deallocate(pmh%pc(1)%el)
!    allocate(pmh%pc(1)%el(ngroups)) ! 1 piece mesh

    if(allocated(pmh%pc(1)%z)) deallocate(pmh%pc(1)%z)
    allocate(pmh%pc(1)%z(pmh%pc(1)%dim,pmh%pc(1)%nnod)) ! 1 piece mesh

    if(pmh%pc(1)%nnod == 0 .or. pmh%pc(1)%dim == 0 .or. nel == 0) then
      call error('Empty mesh or wrong file format.')
    endif


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
  integer                       :: nodnum, lnn, ref, minngroup, ttype, aux
  integer                       :: desc1, desc2, desc3 ! element type descriptors
  integer                       :: i, k, ios, pos
  integer,dimension(15)         :: tel = -1 ! 4:3D(vol),3:2D(sur),2:1D(edge),0:0D(point)


!    ! Associates each group to an element type
    minngroup = 0
    if(nelvol /= 0) then
      minngroup = minngroup + 1
    endif
    if(nelsur /= 0) then
      minngroup = minngroup + 1
    endif
    if(neledge /= 0) then
      minngroup = minngroup + 1
    endif
    if(nelpoint /= 0) then
      minngroup = minngroup + 1
    endif

    ! Allocate pmh%pc(1)%el to the number of different topology dimensions
    if(.not. allocated(pmh%pc(1)%el)) allocate(pmh%pc(1)%el(minngroup))

    do i = 1, nel
      ! Read element description
      read (unit=iu, fmt=*, iostat = ios) nodnum,desc1,desc2,ref,ttype,aux,desc3,lnn
      if (ios /= 0) call error('pf3/read/elements, #'//trim(string(ios)))

      ! Search the position in pmh%pc(1)%el for each type of element
      pos = linear_search(size(tel,1), tel, desc3)

      ! If necessary, allocates pmh%pc(1)%el
      if(pos <= 0) then
        tel(abs(pos)) = desc3
        if(size(pmh%pc(1)%el,1) < abs(pos)) then
          call extend_elgroup(pmh%pc(1)%el, abs(pos))
        endif
      endif

      if(.not. allocated(pmh%pc(1)%el(abs(pos))%nn)) then
        allocate(pmh%pc(1)%el(abs(pos))%nn(lnn,nel))
        if(.not. allocated(pmh%pc(1)%el(abs(pos))%ref)) allocate(pmh%pc(1)%el(abs(pos))%ref(nel))
        pmh%pc(1)%el(abs(pos))%type = pf3_assign_element_type(desc1, desc2, desc3)
      endif

      ! Counts the number of element of each type
      pmh%pc(1)%el(abs(pos))%nel = pmh%pc(1)%el(abs(pos))%nel + 1

      ! Assign reference number
      pmh%pc(1)%el(abs(pos))%ref(pmh%pc(1)%el(abs(pos))%nel) = ref

      ! Read list of nodes
      read (unit=iu, fmt=*, iostat = ios) (pmh%pc(1)%el(abs(pos))%nn(k,pmh%pc(1)%el(abs(pos))%nel), k=1,lnn)
      ! Reorder pf3 nodes to pmh
      call pf32pmh_ordering(pmh%pc(1)%el(abs(pos))%nn(:,pmh%pc(1)%el(abs(pos))%nel), pmh%pc(1)%el(abs(pos))%type)
      if (ios /= 0) call error('pf3/read/elements, #'//trim(string(ios)))

    enddo

    do i = 1, size(pmh%pc(1)%el,1)
      if(.not. allocated(pmh%pc(1)%el(i)%nn) .or. .not. allocated(pmh%pc(1)%el(i)%ref)) then
        call error('Element group without nodes or references.')
      endif
      call reduce(pmh%pc(1)%el(i)%nn,size(pmh%pc(1)%el(i)%nn,1), pmh%pc(1)%el(i)%nel)
      call reduce(pmh%pc(1)%el(i)%ref, pmh%pc(1)%el(i)%nel)
    enddo

end subroutine

!-----------------------------------------------------------------------
! read_pf3_coordinates(iu, nnodk, dim): read PF3 coordinates
!-----------------------------------------------------------------------
! iu:  unit number of the PF3 file
! pc:  piece
!-----------------------------------------------------------------------

subroutine read_pf3_coordinates(iu, pc)

  integer, intent(in) :: iu   ! Unit number for PF3 file
  type(piece), intent(inout) :: pc
  integer :: indx
  integer :: i, j, ios

  if(allocated(pc%z)) deallocate(pc%z)
  allocate(pc%z(pc%dim,pc%nnod))

  do i = 1, pc%nnod
    read (unit=iu, fmt=*, iostat = ios)  indx, (pc%z(j,i), j=1,pc%dim)
    if (ios /= 0) call error('pf3/read/coordinates, #'//trim(string(ios)))
  enddo

  if(.not. allocated(pc%z)) call error('Piece without coordinates.')

end subroutine


subroutine read_pf3_field(iu, pc, line, param)

  integer,                intent(in) :: iu     ! Unit number for PF3 file
  type(piece),         intent(inout) :: pc     ! PMH Piece
  character(len=*),       intent(in) :: line   ! Field header
  real(real64), optional, intent(in) :: param  ! Field shot parameter
  integer                            :: ncomp, npoint
  integer                            :: i, j, ios




  read (unit=iu, fmt=*, iostat = ios)  ncomp, npoint
  if (ios /= 0) call error('pf3/read/field, #'//trim(string(ios)))

  if(npoint /= pc%nnod) call error('pf3/read/field, # number of values and points must agree')

  if(allocated(pc%fi)) deallocate(pc%fi)
  allocate(pc%fi(1))

  pc%fi(1)%name = word(trim(adjustl(line(index(lcase(line),'table of the values of')+len('table of the values of'):))),1)

  call info('Reading node field: '//trim(pc%fi(1)%name))

  if(allocated(pc%fi(1)%param)) deallocate(pc%fi(1)%param)
  allocate(pc%fi(1)%param(1))

  if(present(param)) then
    pc%fi(1)%param(1) = param
  else
    pc%fi(1)%param(1) = 0._real64
  endif

  if(allocated(pc%fi(1)%val)) deallocate(pc%fi(1)%val)
  allocate(pc%fi(1)%val(ncomp,npoint,1))

  do i = 1, pc%nnod
    read (unit=iu, fmt=*, iostat = ios)  (pc%fi(1)%val(j,i,1), j=1,ncomp)
    if (ios /= 0) call error('pf3/read/field, #'//trim(string(ios)))
  enddo

  if(.not. allocated(pc%fi)) call error('Piece without field.')

end subroutine



end module

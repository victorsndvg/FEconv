module module_utils_msh_fcnv

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
use module_alloc_common_bmod, only: search_multiple !basicmod does not directly provides this procedure
use module_pmh_fcnv, only: elgroup,piece,pmh_mesh
use module_fe_database_pmh_fcnv, only: check_fe,FEDB,FA_HEXA
implicit none

integer, parameter, private :: DEFAULT_ALLOC  = 1000 !initial size for allocation

!type msh_faces
!  integer              :: type
!  integer              :: zone
!  integer, allocatable :: mm(:)
!  integer, allocatable :: adcells(:)
!end type

type varlen
  integer, allocatable :: data(:)
end type

type msh_faces
  integer, allocatable :: type(:)
  integer, allocatable :: zone(:)
  type(varlen), allocatable :: mm(:)
  type(varlen), allocatable :: adcells(:)
end type

type msh_cells
  integer, allocatable :: type(:)
  integer, allocatable :: zone(:)
end type

type msh_zone
  integer, allocatable                :: id(:)
  character(len=MAXPATH), allocatable :: names(:)
  character(len=MAXPATH), allocatable :: types(:)
end type

contains



function is_blank_line(line) result(res)
  character(len=*), intent(in) :: line
  logical                      :: res

  res = .False.

  if(len_trim(line) == 0) res = .True.

end function


function get_section_index(line) result(res)
  character(len=*), intent(in) :: line
  character(len=len(line))     :: auxline
  integer                      :: res
  integer                      :: indx
  integer                      :: ios


  auxline(:) = line(:)
  indx = index(line,'(')
  if(indx == 0) then
    res = -1
  else
    call replace_char(auxline, '(', ' ')
    call replace_char(auxline, ')', ' ')
    read(auxline,*, iostat=ios) res
    if(ios /= 0) then
      res = -1
    endif
  endif

end function


subroutine replace_char(line, oldchar, newchar)
  character(len=*), intent(inout) :: line
  character, intent(in)           :: oldchar, newchar
  integer                         :: indx

  do
    indx = index(line,oldchar)
    if(indx == 0) exit
    line(indx:indx) = newchar
  enddo

end subroutine

function count_delimiter(line,del) result(res)
  character(len=*), intent(in) :: line
  character,        intent(in) :: del
  integer                      :: res
  integer                      :: i

  res = 0

  do i=1, len(line)
    if(line(i:i)==del) then
      res = res + 1
    endif
  enddo

end function


!-----------------------------------------------------------------------
! extend_elgroup(v, d): extends the element groups
!-----------------------------------------------------------------------
! v: array of elgroup
! d: new dimension of the array
!-----------------------------------------------------------------------

subroutine extend_elgroup(v, d)
  type(elgroup), allocatable          :: v(:), temp(:)
  integer, intent(in)           :: d !new dimension given by the user
  integer :: res, s, ns
  character(len=MAXPATH) :: cad

    if (.not. allocated(v)) then
      !DIMENSIONS
      ns = d
      !ALLOCATION
      allocate(v(ns), stat = res, errmsg = cad)
      if (res /= 0) call error('module_utils_pf3/extend_elgroup # unable to allocate variable v: '//trim(cad))
    else !v is already allocated
      s = size(v,1)
      if (d > s) then !reallocation is mandatory
        !DIMENSIONS
        ns = d
        !REALLOCATION
        allocate(temp(ns), stat = res, errmsg = cad)
        if (res /= 0) call error('module_utils_pf3/extend_elgroup # unable to allocate variable temp: '//trim(cad))
        temp(1:s)    = v
        call move_alloc(from=temp, to=v)
      end if
   end if
end subroutine


subroutine extend_varlen(v, d, fit)
  type(varlen), allocatable, intent(inout) :: v(:)
  type(varlen), allocatable                :: temp(:)
  integer, intent(in)           :: d !new dimension given by the user
  logical, optional, intent(in)           :: fit !new dimension given by the user
  integer :: res, s, ns
  character(len=MAXPATH) :: cad

    if (.not. allocated(v)) then
      !DIMENSIONS
      if(present(fit)) then
        if (fit) then; ns = d                        !we must fit to dimension given as argument
        else; ns = search_multiple(DEFAULT_ALLOC, d) !a multiple of DEFAULT_ALLOC must be taken as new dimension
        endif
      else
        ns=d
      end if
      !ALLOCATION
      allocate(v(ns), stat = res, errmsg = cad)
      if (res /= 0) call error('module_utils_pf3/extend_varlen # unable to allocate variable v: '//trim(cad))
    else !v is already allocated
      s = max(size(v,1),1)
      if (d > size(v,1)) then !reallocation is mandatory
        !DIMENSIONS
        if(present(fit)) then
          if (fit) then
            ns = d                        !we must fit to dimension given as argument
          else
            ns = search_multiple(s, d) !a multiple of DEFAULT_ALLOC must be taken as new dimension
          endif
        else
          ns=d
        end if
        !REALLOCATION
        allocate(temp(ns), stat = res, errmsg = cad)
        if (res /= 0) call error('module_utils_pf3/extend_varlen # unable to allocate variable temp: '//trim(cad))
        if(size(v,1)>0) temp(1:s) = v
        call move_alloc(from=temp, to=v)
      end if
   end if
   if(allocated(temp)) deallocate(temp)
end subroutine

!-----------------------------------------------------------------------
! reduce: reduce an array of type varlen
!-----------------------------------------------------------------------
subroutine reduce_varlen(v, d)
type(varlen), allocatable :: v(:), temp(:)
integer, intent(in)  :: d

if (size(v,1) == d) return !v has the right size
allocate(temp(d))
temp = v(1:d)
call move_alloc(from=temp, to=v)
end subroutine



function get_cell_pmh_type(tp) result(res)
  integer, intent(in) :: tp
  integer             :: res

  res = -1

  if(tp == 0) then        ! Mixed
    res = check_fe(.true., 0, 0, 0, 0)
  elseif(tp == 1) then    ! Triangle P1
    res = check_fe(.true., 3, 3, 3, 0)
  elseif(tp == 2) then    ! Tetrahedron P1
    res = check_fe(.true., 4, 4, 6, 4)
  elseif(tp == 3) then    ! Quadrangle P1
    res = check_fe(.true., 4, 4, 4, 0)
  elseif(tp == 4) then    ! Hexahedron P1
    res = check_fe(.true., 8, 8, 12, 6)
  elseif(tp == 5) then    ! Pyramid P1
    res = check_fe(.true., 5, 5, 8, 5)
  elseif(tp == 6) then    ! Wedge P1
    res = check_fe(.true., 6, 6, 9, 5)
  endif

end function


function get_face_pmh_type(tp) result(res)
  integer, intent(in) :: tp
  integer             :: res

  res = -1

  if(tp == 0) then        ! Mixed
    res = check_fe(.true., 0, 0, 0, 0)
  elseif(tp == 2) then    ! Edge P1
    res = check_fe(.true., 2, 2, 1, 0)
  elseif(tp == 3) then    ! Triangle P1
    res = check_fe(.true., 3, 3, 3, 0)
  elseif(tp == 4) then    ! Quadrangle P1
    res = check_fe(.true., 4, 4, 4, 0)
  endif

end function

function get_cell_msh_type(tp) result(res)
  integer, intent(in) :: tp
  integer             :: res

  res = -1

  if(tp == check_fe(.true., 0, 0, 0, 0)) then        ! Mixed
    res = 0
  elseif(tp == check_fe(.true., 3, 3, 3, 0)) then    ! Triangle P1
    res = 1
  elseif(tp == check_fe(.true., 4, 4, 6, 4)) then    ! Tetrahedron P1
    res = 2
  elseif(tp == check_fe(.true., 4, 4, 4, 0)) then    ! Quadrangle P1
    res = 3
  elseif(tp == check_fe(.true., 8, 8, 12, 6)) then   ! Hexahedron P1
    res = 4
  elseif(tp == check_fe(.true., 5, 5, 8, 5)) then    ! Pyramid P1
    res = 5
  elseif(tp == check_fe(.true., 6, 6, 9, 5)) then    ! Wedge P1
    res = 6
  endif

end function

subroutine check_interface_names(izones)
  type(msh_zone), intent(inout) :: izones
  integer                       :: i, j, count

  if(allocated(izones%names)) then
    call info('Checking mesh references')
    do i=1,size(izones%names,1)
      count = 0
      if(trim(izones%types(i))=='interior') then
      do j=1,size(izones%names,1)
        if(i == j .or. count>1) cycle
        if(index(trim(izones%names(i)),trim(izones%names(j))) > 0) count = count + 1
        if(count>1) then
          print*, 'Reference '//trim(string(izones%id(i)))//' ('//trim(izones%types(i))//&
            ') named "'//trim(izones%names(i))//'"'
          izones%id(i) = 0
        endif
      enddo
!      if(izones%id(i) /= 0) print*, 'Zone '//trim(string(izones%id(i)))//' named "'//&
!         trim(izones%names(i))//'" of type '//trim(izones%types(i))
      else
        print*, 'Reference '//trim(string(izones%id(i)))//' ('//trim(izones%types(i))//&
          ') named "'//trim(izones%names(i))//'"'
      endif
    enddo
  endif

end subroutine

!subroutine add_faces_to_pmh2(faces,izones,pmh)
!  type(msh_faces), allocatable, intent(inout) :: faces(:) ! msh faces
!  type(msh_zone),               intent(inout) :: izones
!  type(pmh_mesh),               intent(inout) :: pmh
!  integer                                     :: facetypes = 0
!  integer, allocatable                        :: groups(:) !groups(face type) = group number
!  integer                                     :: i, j, aux, numgroups, group, numface
!
!
!  if(.not. allocated(faces)) call error("MSH without faces, could be a binary file")
!  if(.not. allocated(faces(1)%mm)) call error("MSH without faces, could be a binary file")
!  if(.not. allocated(pmh%pc)) call error("PMH without pieces, could be a binary file")
!  if(.not. allocated(pmh%pc(1)%z)) call error("PMH without nodes, could be a binary file")
!
!  if(.not.allocated(groups)) allocate(groups(size(FEDB,1)))
!  groups = 0
!
!  numgroups = 0
!  if(allocated(pmh%pc(1)%el)) numgroups = size(pmh%pc(1)%el,1)
!
!  do i=1, size(faces,1)
!    if (groups(faces(i)%type) == 0) then
!      numgroups = numgroups + 1
!      groups(faces(i)%type) = numgroups
!    endif
!
!    group = groups(faces(i)%type)
!
!    if(.not. allocated(pmh%pc(1)%el)) allocate(pmh%pc(1)%el(group))
!    if(size(pmh%pc(1)%el,1)<group) call extend_elgroup(pmh%pc(1)%el,group)
!    if(pmh%pc(1)%el(group)%type == 0) pmh%pc(1)%el(group)%type = faces(i)%type
!    if(.not. allocated(pmh%pc(1)%el(group)%mm)) then
!      numface = 1
!      allocate(pmh%pc(1)%el(group)%mm(FEDB(faces(i)%type)%lnv,1))
!    else
!      numface = size(pmh%pc(1)%el(group)%mm,2)+1
!      call extend(pmh%pc(1)%el(group)%mm,FEDB(faces(i)%type)%lnv,numface)
!    endif
!      pmh%pc(1)%el(group)%mm(:,numface) = faces(i)%mm(:)
!    if(.not. allocated(pmh%pc(1)%el(group)%mm)) then
!      allocate(pmh%pc(1)%el(group)%mm(FEDB(faces(i)%type)%lnv,1))
!    else
!      call extend(pmh%pc(1)%el(group)%ref,numface)
!    endif
!      print*, group, numface
!      pmh%pc(1)%el(group)%ref(numface) = faces(i)%zone
!  enddo
!
!end subroutine

subroutine add_faces_to_pmh(faces,izones,pmh)
  type(msh_faces),              intent(inout) :: faces ! msh faces
  type(msh_zone),               intent(inout) :: izones
  type(pmh_mesh),               intent(inout) :: pmh
  integer, allocatable                        :: groups(:) !groups(face type) = group number
  integer                                     :: i, numgroups, group, numface,ft
  logical                                     :: tf(2) = [.true.,.false.]


  if(.not. allocated(faces%mm)) call error("MSH without faces, could be a binary file")
  if(.not. allocated(pmh%pc)) call error("PMH without pieces, could be a binary file")
  if(.not. allocated(pmh%pc(1)%z)) call error("PMH without nodes, could be a binary file")

  if(.not.allocated(groups)) allocate(groups(size(FEDB,1)))
  groups = 0

  numgroups = 0
  if(allocated(pmh%pc(1)%el)) numgroups = size(pmh%pc(1)%el,1)

  if (allocated(izones%types)) then
    do i=1,size(izones%id,1)
      if(trim(izones%types(i))/='interior') cycle
      where(faces%zone==izones%id(i))
        faces%zone = 0
      end where
    enddo
  end if  

  if (allocated(faces%type)) then
    do i=1, size(faces%type,1)
      if(faces%zone(i) == 0) cycle
      if (groups(faces%type(i)) == 0) then
        numgroups = numgroups + 1
        groups(faces%type(i)) = numgroups
        if(.not. allocated(pmh%pc(1)%el)) allocate(pmh%pc(1)%el(numgroups))
        if(size(pmh%pc(1)%el,1)<numgroups) then
          call extend_elgroup(pmh%pc(1)%el,numgroups)
          if(pmh%pc(1)%el(numgroups)%type == 0) pmh%pc(1)%el(numgroups)%type = faces%type(i)
          numface = 0
        endif
      endif
  
      ft = faces%type(i)
      group = groups(ft)
  
      numface = numface+1
  
      call set(2, pmh%pc(1)%el(group)%mm,faces%mm(i)%data,numface,fit=tf)
      call set(pmh%pc(1)%el(group)%ref,faces%zone(i), numface,.false.)
      pmh%pc(1)%el(group)%nel = pmh%pc(1)%el(group)%nel + 1
  
    enddo
  end if

  if(allocated(izones%id)) deallocate(izones%id)
  if(allocated(izones%names)) deallocate(izones%names)
  if(allocated(izones%types)) deallocate(izones%types)

end subroutine


subroutine build_cells(faces,cells, pmh)
  type(msh_faces),              intent(inout) :: faces ! msh faces
  type(msh_cells),              intent(inout) :: cells ! msh cells
  type(pmh_mesh),               intent(inout) :: pmh
  integer, allocatable                        :: nodesadded(:)
  integer, allocatable                        :: nels(:)
  integer, allocatable                        :: newnodes(:), posinsubcell(:), posincell(:), auxnodes(:)
  integer, allocatable                        :: groups(:) !groups(face type) = group number
  integer                                     :: i, j, k, l, aux1, aux2, aux3, aux4
  integer                                     :: wedgequadface(4),face(4)
  integer                                     :: numgroups,group,numcell,ft,ct,ac, snn, prevcell
  logical                                     :: tf(2) = [.true.,.false.]


  if(.not. allocated(faces%mm)) call error("MSH without faces, could be a binary file")
  if(.not. allocated(cells%type)) call error("MSH without cells, could be a binary file")
  if(.not. allocated(pmh%pc)) call error("PMH without pieces, could be a binary file")
  if(.not. allocated(pmh%pc(1)%z)) call error("PMH without nodes, could be a binary file")

  if(.not.allocated(groups)) allocate(groups(size(FEDB,1)))
  groups = 0

  if(.not.allocated(nodesadded)) allocate(nodesadded(size(cells%type,1)))
  nodesadded = 0

  if(.not.allocated(nels)) allocate(nels(size(cells%type,1)))
  nels = 0

! Renumber cells for PMH groups
  numcell = 1
  prevcell = 0

  do i=1, size(cells%type,1)
    if(prevcell /= cells%type(i)) numcell = 1
    nels(i) = numcell
    numcell = numcell + 1
    prevcell = cells%type(i)
  enddo

  numgroups = 0
  if(allocated(pmh%pc(1)%el)) numgroups = size(pmh%pc(1)%el,1)

  do i=1, size(faces%type,1)
    do j=1,size(faces%adcells(i)%data,1)
      ac = faces%adcells(i)%data(j)
      if(ac == 0) cycle
      ct = get_cell_pmh_type(cells%type(ac))
      if(nodesadded(ac) == FEDB(ct)%lnv) cycle
      ft = faces%type(i)


      if (groups(ct) == 0) then
        numgroups = numgroups + 1
        groups(ct) = numgroups
      endif

      group = groups(ct)

      if(.not. allocated(pmh%pc(1)%el)) allocate(pmh%pc(1)%el(group))
      if(size(pmh%pc(1)%el,1)<group) then
        call extend_elgroup(pmh%pc(1)%el,group)
        numcell = 0
        call info('Building high order elements group : '//trim(FEDB(ct)%desc))
      endif
      if(pmh%pc(1)%el(group)%type == 0) pmh%pc(1)%el(group)%type = ct

!
!      if(.not. allocated(pmh%pc(1)%el(group)%mm)) then
!        numcell = 1
!        nels(ac) = numcell
!      else
!        if(nels(ac)== 0) then
!          numcell = pmh%pc(1)%el(group)%nel+1
!          nels(ac) = numcell
!        endif
!      endif

      numcell = nels(ac)
      group = groups(ct)
      if(nodesadded(ac) == 0) then
        if(FEDB(pmh%pc(1)%el(group)%type)%tdim == 3) then
          face = FEDB(pmh%pc(1)%el(group)%type)%face(:,1)
        else
          face = [FEDB(pmh%pc(1)%el(group)%type)%edge(:,1),0,0]
        endif
        if(pmh%pc(1)%el(group)%type == check_fe(.true.,6,6,9,5)) then
          if(FEDB(ft)%lnv == 4) then
            if(j == 1 .and. .false.) then ! First adjacent cell: reverse orientation
              do k=1, FEDB(ft)%lnv
                call set(pmh%pc(1)%el(group)%mm, faces%mm(i)%data(k), FEDB(ft)%lnv-k+1, numcell, fit=tf)
              enddo
            else
              do k=1, FEDB(ft)%lnv
                call set(pmh%pc(1)%el(group)%mm, faces%mm(i)%data(k), k, numcell, fit=tf)
              enddo
            endif
            nodesadded(ac) = FEDB(ft)%lnv
          endif
        else
          if(j == 1) then ! First adjacent cell: reverse orientation
            do k=1, FEDB(ft)%lnv
              call set(pmh%pc(1)%el(group)%mm, faces%mm(i)%data(FEDB(ft)%lnv-k+1), face(k), numcell, fit=tf)
            enddo
          else
            do k=1, FEDB(ft)%lnv
              call set(pmh%pc(1)%el(group)%mm, faces%mm(i)%data(k), face(k), numcell, fit=tf)
            enddo
          endif
          nodesadded(ac) = FEDB(ft)%lnv
        endif
      elseif(nodesadded(ac) == FEDB(pmh%pc(1)%el(group)%type)%lnv) then
        cycle
      else
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! Build hexahedrons
        !!!!!!!!!!!!!!!!!!!!!!!!!
        if(pmh%pc(1)%el(group)%type == check_fe(.true.,8,8,12,6)) then
          call nodes_not_in_cell(pmh%pc(1)%el(group)%mm(1:4,numcell),faces%mm(i)%data(:), newnodes, posinsubcell, posincell)
          snn = size(newnodes,1)
          if(snn == 2) then
            call nodes_not_in_cell(pmh%pc(1)%el(group)%mm(5:,numcell),newnodes, auxnodes)
            if(size(auxnodes,1)/=snn) cycle
            call ssort(posinsubcell)
            if((posinsubcell(1) == 1 .and. posinsubcell(2) == 2) .or. &
               (posinsubcell(1) == 3 .and. posinsubcell(2) == 4)) then ! reverse newnodes
              aux1 = newnodes(1); newnodes(1) = newnodes(2); newnodes(2) = aux1
            endif
            do l=1,snn
              call set(pmh%pc(1)%el(group)%mm,newnodes(l),FA_HEXA(posincell(l),4),numcell)
            enddo
            nodesadded(ac) = nodesadded(ac) + snn
          endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! Build wedges
        !!!!!!!!!!!!!!!!!!!!!!!!!
        elseif(pmh%pc(1)%el(group)%type == check_fe(.true.,6,6,9,5)) then
          call nodes_not_in_cell(pmh%pc(1)%el(group)%mm(1:4,numcell),faces%mm(i)%data(:), newnodes, &
                               & posinsubcell, posincell,reverse=.true.)
          snn = size(newnodes,1)
          if(snn == 2 .and. FEDB(ft)%lnv ==4) then
            wedgequadface = pmh%pc(1)%el(group)%mm(1:4,numcell)
!!!!!!!!!!!!!!!!!!!!!! First and second triangle orientation
            if(all(sort(posincell)==posincell).eqv.all(sort(posinsubcell)==posinsubcell)) then
              aux1 = 2; aux2 = 3; aux3 = 3; aux4 = 2
            else
              aux1 = 3; aux2 = 2; aux3 = 2; aux4 = 3
            endif
!!!!!!!!!!!!!!!!!!!!!!! First triangle
            face = FEDB(pmh%pc(1)%el(group)%type)%face(:,1)
            call set(pmh%pc(1)%el(group)%mm,wedgequadface(posincell(1)),face(1),numcell)
            if(mod(posincell(1),4)+1/=posincell(2)) then
              call set(pmh%pc(1)%el(group)%mm,wedgequadface(mod(posincell(1),4)+1),face(aux1),numcell)
            else
              call set(pmh%pc(1)%el(group)%mm,wedgequadface(mod(posincell(1)+2,4)+1),face(aux1),numcell)
            endif
            if(mod(posinsubcell(1),4)+1/=posinsubcell(2)) then
              call set(pmh%pc(1)%el(group)%mm,faces%mm(i)%data(mod(posinsubcell(1),4)+1),face(aux2),numcell)
            else
              call set(pmh%pc(1)%el(group)%mm,faces%mm(i)%data(mod(posinsubcell(1)+2,4)+1),face(aux2),numcell)
            endif
!!!!!!!!!!!!!!!!!!!!!!! Second tringle
          face = FEDB(pmh%pc(1)%el(group)%type)%face(:,4)
            call set(pmh%pc(1)%el(group)%mm,wedgequadface(posincell(2)),4,numcell)
            if(mod(posincell(2),4)+1/=posincell(1)) then
              call set(pmh%pc(1)%el(group)%mm,wedgequadface(mod(posincell(2),4)+1),face(aux3),numcell)
            else
              call set(pmh%pc(1)%el(group)%mm,wedgequadface(mod(posincell(2)+2,4)+1),face(aux3),numcell)
            endif
            if(mod(posinsubcell(2),4)+1/=posinsubcell(1)) then
              call set(pmh%pc(1)%el(group)%mm,faces%mm(i)%data(mod(posinsubcell(2),4)+1),face(aux4),numcell)
            else
              call set(pmh%pc(1)%el(group)%mm,faces%mm(i)%data(mod(posinsubcell(2)+2,4)+1),face(aux4),numcell)
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checking wedges orientation
!print*, numcell,'-', trim(string(pmh%pc(1)%el(group)%mm(:,numcell)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            nodesadded(ac) = nodesadded(ac) + snn
          endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! Build quadrangles
        !!!!!!!!!!!!!!!!!!!!!!!!!
        elseif(pmh%pc(1)%el(group)%type == check_fe(.true.,4,4,4,0)) then

          call nodes_not_in_cell(pmh%pc(1)%el(group)%mm(1:2,numcell),faces%mm(i)%data(:), newnodes, posinsubcell, posincell)
          snn = size(newnodes,1)
          if(snn == 1) then
            call set(pmh%pc(1)%el(group)%mm,newnodes(1),4-posincell(1)+1,numcell)
            nodesadded(ac) = nodesadded(ac) + snn
          endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! Build tetrahedrons, triangles, quadrangles and wedges
        !!!!!!!!!!!!!!!!!!!!!!!!!
        else
          call nodes_not_in_cell(pmh%pc(1)%el(group)%mm(:,numcell),faces%mm(i)%data(:), newnodes)
          snn = size(newnodes,1)
          do k=1,snn
            call set(pmh%pc(1)%el(group)%mm, newnodes(k), FEDB(ft)%lnv+k, numcell, fit=tf)
          enddo
          nodesadded(ac) = nodesadded(ac) + snn
        endif
      endif

      call set(pmh%pc(1)%el(group)%ref,cells%zone(ac), numcell,.false.)
      if(pmh%pc(1)%el(group)%nel < numcell) pmh%pc(1)%el(group)%nel = numcell

    enddo

    if(allocated(faces%mm(i)%data)) deallocate(faces%mm(i)%data)
    if(allocated(faces%adcells(i)%data)) deallocate(faces%adcells(i)%data)

  enddo

end subroutine


subroutine nodes_not_in_cell(cell, subcell, nodes, posinsubcell, posincell, reverse)
  integer,              intent(in)    :: cell(:)
  integer,              intent(in)    :: subcell(:)
  integer, allocatable, intent(inout) :: nodes(:)
  integer                             :: i, j, numnodes, numold
  logical                             :: found
  integer, allocatable, optional      :: posinsubcell(:) ! Contains the position of the different nodes
  integer, allocatable, optional      :: posincell(:) ! Contains the position of coincident nodes
  logical, optional                   :: reverse ! posinsubcell contains the position of coincident nodes in subcell

  if(allocated(nodes)) deallocate(nodes)
  if(present(posinsubcell)) then; if(allocated(posinsubcell)) deallocate(posinsubcell); endif
  if(present(posincell)) then; if(allocated(posincell)) deallocate(posincell); endif
  numnodes = 0
  numold = 0
  do i=1, size(subcell,1)
    found = .false.
    do j=1, size(cell,1)
      if(subcell(i)==cell(j)) then
        found = .true.
        numold = numold + 1
        if(present(posincell)) call insert(posincell,j,numold)
        if(present(reverse)) then
          if(reverse) call insert(posinsubcell,i,numold)
        endif
        exit
      endif
    enddo
    if(.not. found) then

      numnodes = numnodes + 1
      call set(nodes,subcell(i),numnodes,.false.)
      if(.not. present(reverse)) then
        if(present(posinsubcell)) call set(posinsubcell,i,numnodes,.false.)
      elseif(.not. reverse) then
        if(present(posinsubcell)) call set(posinsubcell,i,numnodes,.false.)
      endif
    endif
  enddo
  call reduce(nodes,numnodes)
  if(present(posinsubcell) .and. allocated(posinsubcell)) then
    if(present(reverse)) then
      if(reverse) then
        call reduce(posinsubcell, numold)
      else
        call reduce(posinsubcell,numnodes)
      endif
    else
      call reduce(posinsubcell,numnodes)
    endif
  endif
  if(present(posincell) .and. allocated(posincell)) call reduce(posincell,numold)
!  if(present(oldpos) .and. allocated(oldpos)) print*, oldpos

end subroutine


subroutine face_positive_jacobian(mm,z)
  integer, intent(inout) :: mm(:)
  real(real64), allocatable,    intent(in) :: z(:,:)
  integer, allocatable :: pos(:)
  real(real64) ::  s, t
  logical :: QJac(4)


  if(.not. allocated(z)) call error("Positive jacobian: coordinates not allocated")

  if(size(mm,1) == 3) then ! Triangle P1
    if ( (z(1,mm(2)) - z(1,mm(1))) * (z(2,mm(3)) - z(2,mm(1))) & !only x and y coordinates are used
       - (z(2,mm(2)) - z(2,mm(1))) * (z(1,mm(3)) - z(1,mm(1))) < 0 ) &
    call swap(mm(2), mm(3))
  elseif(size(mm,1) == 4) then ! Quadrangle P1
    s = dot_product(z(:,mm(3))-z(:,mm(1)), z(:,mm(2))-z(:,mm(1))) !projection of a3-a1 in <{a2-a1}>
    t = dot_product(z(:,mm(3))-z(:,mm(1)), z(:,mm(4))-z(:,mm(1))) !projection of a3-a1 in <{a4-a1}>
    QJac = [QJ_pos([0._real64, 0._real64], [1._real64, 0._real64], [s, t], [0._real64, 1._real64], -1, -1), &
            QJ_pos([0._real64, 0._real64], [1._real64, 0._real64], [s, t], [0._real64, 1._real64],  1, -1), &
            QJ_pos([0._real64, 0._real64], [1._real64, 0._real64], [s, t], [0._real64, 1._real64],  1,  1), &
            QJ_pos([0._real64, 0._real64], [1._real64, 0._real64], [s, t], [0._real64, 1._real64], -1,  1)]
    call sfind(QJac, .false., pos) !check whether vertices in (s,t)-plane are well oriented
    if(size(pos,1) == 2) then
      call swap(mm(pos(1)), mm(pos(2)))
    elseif (size(pos,1) == 4) then !in theory, never achieved; it remains for security reasons
      call swap(mm(1), mm(2))
      call swap(mm(3), mm(4))
    end if
    !counterclockwise orientation is not ensured by a positive jacobian in the (s,t)-plane: instead, check if 3rd component
    !of the normal vector is positive
    if ((z(1,mm(2))-z(1,mm(1)))*(z(2,mm(4))-z(2,mm(1))) - &
      (z(2,mm(2))-z(2,mm(1)))*(z(1,mm(4))-z(1,mm(1))) < 0._real64) then
      call swap(mm(1), mm(2))
      call swap(mm(3), mm(4))
    end if
  endif

end subroutine

!-----------------------------------------------------------------------
! swap: swap the value of two variables
!-----------------------------------------------------------------------
subroutine swap(a,b)
integer :: a,b,c
c = a; a = b; b = c
end subroutine

!-----------------------------------------------------------------------
! QJ_pos: check whether the 4-node quadrilateral jacobian is positive
! Only (x,y) coordinates are used
!-----------------------------------------------------------------------
function QJ_pos(a1, a2, a3, a4, r, s) result(res)
real(real64), intent(in) :: a1(2), a2(2), a3(2), a4(2)
integer :: r, s
logical :: res

res = (a4(2)-a2(2))*(a3(1)-a1(1))-(a3(2)-a1(2))*(a4(1)-a2(1))+ &
     ((a3(2)-a4(2))*(a2(1)-a1(1))-(a2(2)-a1(2))*(a3(1)-a4(1)))*r+ &
     ((a4(2)-a1(2))*(a3(1)-a2(1))-(a3(2)-a2(2))*(a4(1)-a1(1)))*s > 0
end function


! Bubble sort code modification from:
! http://jean-pierre.moreau.pagesperso-orange.fr/f_sort.html
! return p,q in ascending order
subroutine order(p,q)
 integer, intent(inout) :: p,q
 integer                :: temp

  if (p>q) then
    temp=p
    p=q
    q=temp
  end if
  return

end subroutine

! Bubble sort code modification from:
! http://jean-pierre.moreau.pagesperso-orange.fr/f_sort.html
! Bubble sorting of integer array A
function bubble(A) result(R)
  integer, intent(in)  :: A(:)
  integer, allocatable :: R(:)
  integer              :: i, j, n

  n = size(A)
  if(.not. allocated(R)) allocate(R(n))
  R(:) = A(:)

  do i=1, n
    do j=n, i+1, -1
      call order(R(j-1), R(j))
    end do
  end do
  return
end function

!!-----------------------------------------------------------------------
!! PRIVATE PROCEDURES
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!! search_multiple: search the smallest value of 2 to the power of a that is bigger than b
!! 2^n*a > b  <=>  n > log2(b/a)
!!-----------------------------------------------------------------------
!integer function search_multiple(a,b)
!integer, intent(in) :: a, b
!
!if (b > a) then
!  search_multiple = int(2**real(ceiling(log(real(b)/a)/log(2.)))*a)
!else
!  search_multiple = a
!end if
!end function

end module

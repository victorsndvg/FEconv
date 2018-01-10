module module_utils_mphtxt_fcnv

!-----------------------------------------------------------------------
! Module to manage MPHTXT (Comsol) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 21/02/2014
!
! PUBLIC PROCEDURES:
! mphtxt_get_lnn:  return  the number of nodes for each MPHTXT FE type
! mphtxt_get_type: return the FE type from FEDB given the MPHTXT FE descriptor
! mphtxt_get_desc:  return the MPHTXT FE descriptor given the FE type from FEDB
! pmh_node_ordering: translates the node ordering from MPHTXT to PMH
! mphtxt_node_ordering: translates the node ordering from PMH to MPHTXT
! build_elements_baricenter: build the baricenter of the given elements
! build_faces_baricenter: build the baricenter of the faces of the given elements
!-----------------------------------------------------------------------

use basicmod, only: insert, reduce
use module_pmh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! mphtxt_get_lnn(num): return the number of nodes for each MPHTXT FE type
!-----------------------------------------------------------------------
! num: number describing the element type stored in FEDB
!-----------------------------------------------------------------------

function mphtxt_get_lnn(num) result(res)

  integer, intent(in) :: num
  integer             :: res

    res = 0
    if(num == check_fe(.false., 8, 4, 4, 0)) then        ! Quadrangle Lagrange P2
      res = 9                                            ! Quadrangle Lagrange P2 + baricentro
    elseif(num == check_fe(.false.,  20, 8,12, 6)) then  ! Hexahedron Lagrange P2
      res = 27                                           ! Hexahedron Lagrange P2 + baricentro + baricentros caras
    elseif(num == check_fe(.false., 15, 6,  9, 5)) then  ! Prism Lagrange P2
      res = 18                                           ! Prism Lagrange P2 + baricentros caras quads
    else
      res = FEDB(num)%lnn

  endif

end function


!-----------------------------------------------------------------------
! mphtxt_get_type(desc): return the FE type from FEDB given the MPHTXT FE descriptor
!-----------------------------------------------------------------------
! desc: MPHTXT FE descriptor
!-----------------------------------------------------------------------

function mphtxt_get_type(desc) result(res)

  character(len=*), intent(in) :: desc
  integer :: res
  integer :: nnod, nver, lnn, lnv, lne, lnf

    nnod=0; nver=0; lnn=0; lnv=0; lne=0; lnf=0

    if(trim(desc) == 'vtx') then        ! Node
      nnod=1; nver=1; lnn=1; lnv=1; lne=0; lnf=0
      call info('Element type: Node')
    elseif(trim(desc) == 'edg') then    ! Edge Lagrange P1
      nnod=2; nver=2; lnn=2; lnv=2; lne=1; lnf=0
      call info('Element type: Edge lagrange P1')
    elseif(trim(desc) == 'tri') then    ! Triangle Lagrange P1
      nnod=3; nver=3; lnn=3; lnv=3; lne=3; lnf=0
      call info('Element type: Triangle lagrange P1')
    elseif(trim(desc) == 'quad') then   ! Quadrangle Lagrange P1
      nnod=4; nver=4; lnn=4; lnv=4; lne=4; lnf=0
      call info('Element type: Quadrangle lagrange P1')
    elseif(trim(desc) == 'tet') then    ! Tetrahedron Lagrange P1
      nnod=4; nver=4; lnn=4; lnv=4; lne=6; lnf=4
      call info('Element type: Tetrahedron lagrange P1')
    elseif(trim(desc) == 'prism') then  ! Wedge Lagrange P1
      nnod=6; nver=6; lnn=6; lnv=6; lne=9; lnf=5
      call info('Element type: Wedge lagrange P1')
    elseif(trim(desc) == 'hex') then    ! Hexahedron Lagrange P1
      nnod=8; nver=8; lnn=8; lnv=8; lne=12; lnf=6
      call info('Element type: Hexahedron lagrange P1')
    elseif(trim(desc) == 'edg2') then   ! Edge Lagrange P2
      nnod=3; nver=2; lnn=3; lnv=2; lne=1; lnf=0
      call info('Element type: Edge lagrange P2')
    elseif(trim(desc) == 'tri2') then   ! Triangle Lagrange P2
      nnod=6; nver=3; lnn=6; lnv=3; lne=3; lnf=0
      call info('Element type: Triangle lagrange P2')
    elseif(trim(desc) == 'quad2') then  ! Quadrangle Lagrange P2
      nnod=8; nver=4; lnn=8; lnv=4; lne=4; lnf=0
      call info('Element type: Quadrangle lagrange P2')
    elseif(trim(desc) == 'tet2') then   ! Tetrahedron Lagrange P2
      nnod=10; nver=4; lnn=10; lnv=4; lne=6; lnf=4
      call info('Element type: Tetrahedron lagrange P2')
    elseif(trim(desc) == 'prism2') then ! Wedge Lagrange P2
      nnod=15; nver=6; lnn=15; lnv=6; lne=9; lnf=5
      call info('Element type: Wedge lagrange P2')
    elseif(trim(desc) == 'hex2') then   ! Hexahedron Lagrange P2
      nnod=20; nver=8; lnn=20; lnv=8; lne=12; lnf=6
      call info('Element type: Hexahedron lagrange P2')
  endif

res = check_fe(nnod==nver, lnn, lnv, lne, lnf)

end function


!-----------------------------------------------------------------------
! mphtxt_get_desc(num):  return the MPHTXT FE descriptor given the FE type from FEDB
!-----------------------------------------------------------------------
! num: number describing the element type stored in FEDB
!-----------------------------------------------------------------------

function mphtxt_get_desc(num) result(res)

  integer, intent(in) :: num
  character(len=MAXPATH) :: res

    res = ''


    if(num == check_fe(.true., 1, 1, 0, 0)) then      ! Node
      res = 'vtx'
      call info('Element type: Node')
    elseif(num == check_fe(.true., 2, 2, 1, 0)) then    ! Edge Lagrange P1
      res = 'edg'
      call info('Element type: Edge lagrange P1')
    elseif(num == check_fe(.true., 3, 3, 3, 0)) then    ! Triangle Lagrange P1
      res = 'tri'
      call info('Element type: Triangle lagrange P1')
    elseif(num == check_fe(.true., 4, 4, 4, 0)) then    ! Quadrangle Lagrange P1
      res = 'quad'
      call info('Element type: Quadrangle lagrange P1')
    elseif(num == check_fe(.true., 4, 4, 6, 4)) then    ! Tetrahedron Lagrange P1
      res = 'tet'
      call info('Element type: Tetrahedron lagrange P1')
    elseif(num == check_fe(.true., 6, 6, 9, 5)) then    ! Wedge Lagrange P1
      res = 'prism'
      call info('Element type: Prism lagrange P1')
    elseif(num == check_fe(.true., 15, 6, 9, 5)) then    ! Wedge Lagrange P2
      res = 'prism'
      call info('Element type: Prism lagrange P2')
    elseif(num == check_fe(.true., 8, 8, 12, 6)) then   ! Hexahedron Lagrange P1
      res = 'hex'
      call info('Element type: Hexahedron lagrange P1')
    elseif(num == check_fe(.false., 3, 2, 1, 0)) then   ! Edge Lagrange P2
      res = 'edg2'
      call info('Element type: Edge lagrange P2')
    elseif(num == check_fe(.false., 6, 3, 3, 0)) then   ! Triangle Lagrange P2
      res = 'tri2'
      call info('Element type: Triangle lagrange P2')
    elseif(num == check_fe(.false., 8, 4, 4, 0)) then   ! Quadrangle Lagrange P2
      res = 'quad2'
      call info('Element type: Quadrangle lagrange P2')
    elseif(num == check_fe(.false., 10, 4, 6, 4)) then  ! Tetrahedron Lagrange P2
      res = 'tet2'
      call info('Element type: Tetrahedron lagrange P2')
    elseif(num == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron Lagrange P2
      res = 'hex2'
      call info('Element type: Hexahedron lagrange P2')
    else
      call error('Finite element type not supported')
  endif

end function


!-----------------------------------------------------------------------
! pmh_node_ordering(el, tp):  translates the node ordering from MPHTXT to PMH
!-----------------------------------------------------------------------
! el: array of connectivities from an element stored in nn
! tp: FE type identifier from FEDB
!-----------------------------------------------------------------------

subroutine pmh_node_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer, intent(in) :: tp
  integer :: aux
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_read_mphtxt/node_ordering # Element type not supported')
    endif

    if (tp == check_fe(.true., 1, 1, 0, 0)) then     ! Nodes
        ! PMH and MPHTXT uses the same node ordering in nodes

    elseif (tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P1

    elseif (tp == check_fe(.true., 3, 3, 3, 0)) then ! Triangle Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in triangles lagrange P1

    elseif (tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle Lagrange P1
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = MPH[1,2,4,3]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux

    elseif (tp == check_fe(.true., 4, 4, 6, 4)) then ! Tetrahedron Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in tetrahedrons lagrange P1

    elseif (tp == check_fe(.true., 6, 6, 9, 5)) then ! Wedge Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in wedges lagrange P1

    elseif (tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron Lagrange P1
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,2,4,3,5,6,8,7]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(7); el(7) = aux

    elseif (tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge Lagrange P2
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P2

    elseif (tp == check_fe(.false., 6, 3, 3, 0)) then ! Triangle Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in triangles lagrange P2
        ! PMH[1,2,3,4,5,6] = MPH[1,2,3,4,6,5]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
          aux = el(6); el(6) = el(5); el(5) = aux

    elseif (tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadragle Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P2

        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(6); el(6) = aux
        aux = el(9); el(9) = el(7); el(7) = aux

    elseif (tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10] = MPH[1,2,4,3,5,7,6,8,9,10]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(7); el(7) = el(6); el(6) = aux

    elseif (tp == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] =
        ! MPH[1,2,4,3,5,6,8,7,9,12,21,10,11,13,23,14,22,27,25,16,26,15,17,20,24,18,19]



        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif

        if (allocated(auxel)) deallocate(auxel)
        allocate(auxel(size(el,1)))
        auxel(:) = el(:)

        el(1) = auxel(1); el(2) = auxel(2); el(3) = auxel(4); el(4) = auxel(3)
        el(5) = auxel(5); el(6) = auxel(6); el(7) = auxel(8); el(8) = auxel(7)
        el(9) = auxel(9); el(10) = auxel(12); el(11) = auxel(13); el(12) = auxel(10)
        el(13) = auxel(14); el(14) = auxel(16); el(15) = auxel(22); el(16) = auxel(20)
        el(17) = auxel(23); el(18) = auxel(26); el(19) = auxel(27); el(20) = auxel(24)
        el(21) = auxel(11); el(22) = auxel(17); el(23) = auxel(15); el(24) = auxel(25)
        el(25) = auxel(19); el(26) = auxel(21); el(27) = auxel(18)

        deallocate(auxel)

    elseif (tp == check_fe(.true., 15, 6, 9, 5)) then ! Wedge Lagrange P2
        ! PMH and MPHTXT uses the same node ordering in wedges lagrange P2


    endif


end subroutine


!-----------------------------------------------------------------------
! mphtxt_node_ordering(el, tp):  translates the node ordering from PMH to MPHTXT
!-----------------------------------------------------------------------
! el: array of connectivities from an element stored in nn
! tp: FE type identifier from FEDB
!-----------------------------------------------------------------------

subroutine mphtxt_node_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer, intent(in) :: tp
  integer :: aux
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_read_mphtxt/node_ordering # Element type not supported')
    endif

    if (tp == check_fe(.true., 1, 1, 0, 0)) then     ! Nodes
        ! PMH and MPHTXT uses the same node ordering in nodes

    elseif (tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P1

    elseif (tp == check_fe(.true., 3, 3, 3, 0)) then! Triangle Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in triangles lagrange P1

    elseif (tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle Lagrange P1
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = MPH[1,2,4,3]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux

    elseif (tp == check_fe(.true., 4, 4, 6, 4)) then ! Tetrahedron Lagrange P1
        ! PMH and MPHTXT uses the same node ordering in tetrahedrons lagrange P1

    elseif (tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron Lagrange P1
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,2,4,3,5,6,8,7]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(7); el(7) = aux

    elseif (tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge Lagrange P2
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P2

    elseif (tp == check_fe(.false., 6, 3, 3, 0)) then ! Triangle Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in triangles lagrange P2
        ! PMH[1,2,3,4,5,6] = MPH[1,2,3,4,6,5]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
          aux = el(6); el(6) = el(5); el(5) = aux

    elseif (tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadragle Lagrange P2
        ! PMH and MPHTXT uses the same node ordering in quadrangles lagrange P2
        if (size(el,1) /= 9) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(6); el(6) = aux
        aux = el(9); el(9) = el(7); el(7) = aux

    elseif (tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10] = MPH[1,2,4,3,5,7,6,8,9,10]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(7); el(7) = el(6); el(6) = aux

    elseif (tp == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron Lagrange P2
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] =
        ! MPH[1,2,4,3,5,6,8,7,9,12,21,10,11,13,23,14,22,27,25,16,26,15,17,20,24,18,19]

        if (size(el,1) /= 27) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif

        if (allocated(auxel)) deallocate(auxel)
        allocate(auxel(size(el,1)))
        auxel(:) = el(:)

        el(1) = auxel(1); el(2) = auxel(2); el(3) = auxel(4); el(4) = auxel(3)
        el(5) = auxel(5); el(6) = auxel(6); el(7) = auxel(8); el(8) = auxel(7)
        el(9) = auxel(9); el(10) = auxel(12); el(11) = auxel(21); el(12) = auxel(10)
        el(13) = auxel(11); el(14) = auxel(13); el(15) = auxel(23); el(16) = auxel(14)
        el(17) = auxel(22); el(18) = auxel(27); el(19) = auxel(25); el(20) = auxel(16)
        el(21) = auxel(26); el(22) = auxel(15); el(23) = auxel(17); el(24) = auxel(20)
        el(25) = auxel(24); el(26) = auxel(18); el(27) = auxel(19)

        deallocate(auxel)

    endif


end subroutine


!-----------------------------------------------------------------------
! build_elements_baricenter(pc, ip, ie, znod): build the baricenter of the given elements
!-----------------------------------------------------------------------
! pc:   piece of the mesh
! ip:   index of the piece
! ie:   index of the element group
! znod: array of node coordinates
!-----------------------------------------------------------------------

subroutine build_elements_baricenter(pc, ip, ie, znod)

  type(piece),               intent(inout) :: pc ! Piece
  integer,                   intent(in)    :: ip ! Piece number
  integer,                   intent(in)    :: ie ! Element group number
  real(real64),dimension(:,:),allocatable, intent(inout) :: znod(:,:)
  integer                                  :: i, j, tp, mphlnn, indx, maxindx
  integer, dimension(:,:), allocatable     :: nn
  real(real64),dimension(:), allocatable   :: val

    tp = pc%el(ie)%type
    mphlnn = mphtxt_get_lnn(tp)
    if(size(pc%el(ie)%nn,1) /= mphlnn) then
      if(allocated(nn)) deallocate(nn); allocate(nn(mphlnn,pc%el(ie)%nel)); nn = 0
      nn(1:size(pc%el(ie)%nn,1),:) = pc%el(ie)%nn(:,:)
      call move_alloc(from=nn,  to=pc%el(ie)%nn)
      if(allocated(nn)) deallocate(nn)
    endif

    if(allocated(val)) deallocate(val); allocate(val(pc%dim))

    do i = 1, pc%el(ie)%nel
      do j = 1, pc%dim
        val(j) = sum(pc%z(j, pc%el(ie)%mm(:,i))) / FEDB(tp)%lnv
      enddo
      if (pc%el(ie)%nn(mphlnn,i) == 0) pc%el(ie)%nn(mphlnn,i) = size(znod,1)+1
      indx = pc%el(ie)%nn(mphlnn,i)
      if(indx < size(znod,2)) then; call set(2, znod, val, indx)
      else; call insert(2, znod, val, indx); endif
      if(maxindx<max(indx,size(znod,2))) maxindx=max(indx,size(znod,2))
    end do

    if(size(znod,2)>maxindx) call reduce(znod, maxindx, pc%el(ie)%nel)
    call info('Added element baricenters in piece '//trim(string(ip))//' element group '//trim(string(ie)))

end subroutine


!-----------------------------------------------------------------------
! build_faces_baricenter(pc, ip, ie, znod): build the faces baricenter of the given elements
!-----------------------------------------------------------------------
! pc:   piece of the mesh
! ip:   index of the piece
! ie:   index of the element group
! znod: array of node coordinates
!-----------------------------------------------------------------------

subroutine build_faces_baricenter(pc, ip, ie, znod)

  type(piece),               intent(inout) :: pc ! Piece
  integer,                   intent(in)    :: ip ! Piece number
  integer,                   intent(in)    :: ie ! Element group number
  real(real64),dimension(:,:),allocatable, intent(inout) :: znod(:,:)
  integer, dimension(:,:), allocatable     :: nn
  real(real64),dimension(:), allocatable   :: val
  integer :: i, j, k, tp, indx, maxindx, mphlnn

    tp = pc%el(ie)%type
    mphlnn = mphtxt_get_lnn(tp)
    if(size(pc%el(ie)%nn,1) /= mphlnn) then
      if(allocated(nn)) deallocate(nn); allocate(nn(mphlnn,pc%el(ie)%nel)); nn = 0
      nn(1:size(pc%el(ie)%nn,1),:) = pc%el(ie)%nn(:,:)
      call move_alloc(from=nn,  to=pc%el(ie)%nn)
      deallocate(nn)
    endif

    if(allocated(val)) deallocate(val); allocate(val(pc%dim))
    val = 0

    do k = 1, pc%el(ie)%nel
      do i = 1, FEDB(tp)%lnf
        do j = 1, pc%dim
          indx = FEDB(tp)%lnv+FEDB(tp)%lne+i
          if (pc%el(ie)%nn(indx,k) == 0) pc%el(ie)%nn(indx,k) = size(znod,1)+1
          val(j) = sum(pc%z(j, pc%el(ie)%mm(FEDB(tp)%face(:,i),k))) / FEDB(FEDB(tp)%f_type)%lnv
          if(pc%el(ie)%nn(indx,k) < size(znod,2)) then; call set(2, znod, val, pc%el(ie)%nn(indx,k))
          else; call insert(2, znod, val, pc%el(ie)%nn(indx,k)); endif
          if(maxindx<max(pc%el(ie)%nn(indx,k),size(znod,2))) maxindx=max(pc%el(ie)%nn(indx,k),size(znod,2))
        end do
      end do
    end do


    if(size(znod,2)>maxindx) call reduce(znod, maxindx, pc%el(ie)%nel)
    call info('Added face baricenters in piece '//trim(string(ip))//' element group '//trim(string(ie)))

end subroutine



end module

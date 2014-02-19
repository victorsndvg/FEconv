module module_utils_mphtxt
use module_alloc, only:set_col,reduce
use module_pmh

contains

function mphtxt_get_lnn(num) result(res)

  integer, intent(in) :: num
  integer             :: res

    res = 0
    if((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (8==FEDB(num)%lnn) .and. &      ! Quadrangle Lagrange P2
        (4==FEDB(num)%lnv) .and. (4==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 9                                                                       ! Quadrangle Lagrange P2 + centroide
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (20==FEDB(num)%lnn) .and. & ! Hexahedron Lagrange P2
        (8==FEDB(num)%lnv) .and. (12==FEDB(num)%lne) .and. (6==FEDB(num)%lnf)) then
      res = 27                                                                      ! Hexahedron Lagrange P2 + centroide + baricentros caras
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (15==FEDB(num)%lnn) .and. & ! Prism Lagrange P2
        (6==FEDB(num)%lnv) .and. (9==FEDB(num)%lne) .and. (5==FEDB(num)%lnf)) then
      res = 18                                                                      ! Prism Lagrange P2 + centroide + baricentros caras
    else
      res = FEDB(num)%lnn

  endif

end function

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
    elseif(trim(desc) == 'prism') then  ! Prism Lagrange P1
      ! Prism FE not supported
      call error('Wedge lagrange P1 not supported')
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
    elseif(trim(desc) == 'prism2') then ! Prism Lagrange P2
      ! Quadratic prism FE not supported
      call error('Wedge lagrange P2 not supported')
    elseif(trim(desc) == 'hex2') then   ! Hexahedron Lagrange P2
      nnod=20; nver=8; lnn=20; lnv=8; lne=12; lnf=6
      call info('Element type: Hexahedron lagrange P2')
  endif

res = check_fe(nnod==nver, lnn, lnv, lne, lnf)

end function


function mphtxt_get_desc(num) result(res)

  integer, intent(in) :: num
  character(len=MAXPATH) :: res

    res = ''


    if((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (1==FEDB(num)%lnn) .and. &     ! Node
        (1==FEDB(num)%lnv) .and. (0==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'vtx'
      call info('Element type: Node')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (2==FEDB(num)%lnn) .and. & ! Edge Lagrange P1
        (2==FEDB(num)%lnv) .and. (1==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'edg'
      call info('Element type: Edge lagrange P1')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (3==FEDB(num)%lnn) .and. & ! Triangle Lagrange P1
        (3==FEDB(num)%lnv) .and. (3==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'tri'
      call info('Element type: Triangle lagrange P1')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (4==FEDB(num)%lnn) .and. & ! Quadrangle Lagrange P1
        (4==FEDB(num)%lnv) .and. (4==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'quad'
      call info('Element type: Quadrangle lagrange P1')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (4==FEDB(num)%lnn) .and. & ! Tetrahedron Lagrange P1
        (4==FEDB(num)%lnv) .and. (6==FEDB(num)%lne) .and. (4==FEDB(num)%lnf)) then
      res = 'tet'
      call info('Element type: Tetrahedron lagrange P1')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .true.) .and. (8==FEDB(num)%lnn) .and. & ! Hexahedron Lagrange P1
        (8==FEDB(num)%lnv) .and. (12==FEDB(num)%lne) .and. (6==FEDB(num)%lnf)) then
      res = 'hex'
      call info('Element type: Hexahedron lagrange P1')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (3==FEDB(num)%lnn) .and. & ! Edge Lagrange P2
        (2==FEDB(num)%lnv) .and. (1==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'edg2'
      call info('Element type: Edge lagrange P2')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (6==FEDB(num)%lnn) .and. & ! Triangle Lagrange P2
        (3==FEDB(num)%lnv) .and. (3==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'tri2'
      call info('Element type: Triangle lagrange P2')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (8==FEDB(num)%lnn) .and. & ! Quadrangle Lagrange P2
        (4==FEDB(num)%lnv) .and. (4==FEDB(num)%lne) .and. (0==FEDB(num)%lnf)) then
      res = 'quad2'
      call info('Element type: Quadrangle lagrange P2')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (10==FEDB(num)%lnn) .and. & ! Tetrahedron Lagrange P2
        (4==FEDB(num)%lnv) .and. (6==FEDB(num)%lne) .and. (4==FEDB(num)%lnf)) then
      res = 'tet2'
      call info('Element type: Tetrahedron lagrange P2')
    elseif((FEDB(num)%nver_eq_nnod .eqv. .false.) .and. (20==FEDB(num)%lnn) .and. & ! Hexahedron Lagrange P2
        (8==FEDB(num)%lnv) .and. (12==FEDB(num)%lne) .and. (6==FEDB(num)%lnf)) then
      res = 'hex2'
      call info('Element type: Hexahedron lagrange P2')
    else
      call error('Finite element type not supported')
  endif


end function


subroutine pmh_node_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer, intent(in) :: tp
  integer :: aux
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_read_mphtxt/node_ordering # Element type not supported')
    endif

    if ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. &     ! Nodes
        (FEDB(tp)%lnn == 1) .and. (FEDB(tp)%lnv == 1) .and. &
        (FEDB(tp)%lne == 0) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in nodes

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Edge Lagrange P1
            (FEDB(tp)%lnn == 2) .and. (FEDB(tp)%lnv == 2) .and. &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Triangle Lagrange P1
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 3) .and. &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in triangles lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Quadrangle Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = MPH[1,2,4,3]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Tetrahedron Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 6) .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT uses the same node ordering in tetrahedrons lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Hexahedron Lagrange P1
            (FEDB(tp)%lnn == 8) .and. (FEDB(tp)%lnv == 8) .and. &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,2,4,3,5,6,8,7]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(7); el(7) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Edge Lagrange P2
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 2) .and. &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P2

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Triangle Lagrange P2
            (FEDB(tp)%lnn == 6) .and. (FEDB(tp)%lnv == 3) .and. &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in triangles lagrange P2
        ! PMH[1,2,3,4,5,6] = MPH[1,2,3,4,6,5]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
          aux = el(6); el(6) = el(5); el(5) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Quadragle Lagrange P2
            (FEDB(tp)%lnn == 8) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P2

        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(6); el(6) = aux
        aux = el(9); el(9) = el(7); el(7) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Tetrahedron Lagrange P2
            (FEDB(tp)%lnn == 10) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 6) .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10] = MPH[1,2,4,3,5,7,6,8,9,10]
        if (size(el,1) /= mphtxt_get_lnn(tp)) then
          call error('module_read_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(7); el(7) = el(6); el(6) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Hexahedron Lagrange P2
            (FEDB(tp)%lnn == 20) .and. (FEDB(tp)%lnv == 8) .and. &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] =
        ! MPH[1,2,4,3,5,6,8,7,9,12,13,10,14,16,22,20,23,26,27,24,11,17,15,25,19,21,27]

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

    endif


end subroutine


subroutine mphtxt_node_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer, intent(in) :: tp
  integer :: aux
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_read_mphtxt/node_ordering # Element type not supported')
    endif

    if ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. &     ! Nodes
        (FEDB(tp)%lnn == 1) .and. (FEDB(tp)%lnv == 1) .and. &
        (FEDB(tp)%lne == 0) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in nodes

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Edge Lagrange P1
            (FEDB(tp)%lnn == 2) .and. (FEDB(tp)%lnv == 2) .and. &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Triangle Lagrange P1
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 3) .and. &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in triangles lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Quadrangle Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = MPH[1,2,4,3]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Tetrahedron Lagrange P1
            (FEDB(tp)%lnn == 4) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 6) .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT uses the same node ordering in tetrahedrons lagrange P1

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .true.) .and. & ! Hexahedron Lagrange P1
            (FEDB(tp)%lnn == 8) .and. (FEDB(tp)%lnv == 8) .and. &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,2,4,3,5,6,8,7]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(7); el(7) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Edge Lagrange P2
            (FEDB(tp)%lnn == 3) .and. (FEDB(tp)%lnv == 2) .and. &
            (FEDB(tp)%lne == 1) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in edges lagrange P2

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Triangle Lagrange P2
            (FEDB(tp)%lnn == 6) .and. (FEDB(tp)%lnv == 3) .and. &
            (FEDB(tp)%lne == 3) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT don't have the same node ordering in triangles lagrange P2
        ! PMH[1,2,3,4,5,6] = MPH[1,2,3,4,6,5]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
          aux = el(6); el(6) = el(5); el(5) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Quadragle Lagrange P2
            (FEDB(tp)%lnn == 8) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 4) .and. (FEDB(tp)%lnf == 0)) then
        ! PMH and MPHTXT uses the same node ordering in quadrangles lagrange P2
        if (size(el,1) /= 9) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(4); el(4) = el(3); el(3) = aux
        aux = el(8); el(8) = el(6); el(6) = aux
        aux = el(9); el(9) = el(7); el(7) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Tetrahedron Lagrange P2
            (FEDB(tp)%lnn == 10) .and. (FEDB(tp)%lnv == 4) .and. &
            (FEDB(tp)%lne == 6) .and. (FEDB(tp)%lnf == 4)) then
        ! PMH and MPHTXT don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10] = MPH[1,2,4,3,5,7,6,8,9,10]
        if (size(el,1) /= FEDB(tp)%lnn) then
          call error('module_write_mphtxt/node_ordering # Wrong element size' )
        endif
        aux = el(7); el(7) = el(6); el(6) = aux

    elseif ((FEDB(tp)%nver_eq_nnod .eqv. .false.) .and. & ! Hexahedron Lagrange P2
            (FEDB(tp)%lnn == 20) .and. (FEDB(tp)%lnv == 8) .and. &
            (FEDB(tp)%lne == 12) .and. (FEDB(tp)%lnf == 6)) then
        ! PMH and MPHTXT don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] =
        ! MPH[1,2,4,3,5,6,8,7,9,12,13,10,14,16,22,20,23,26,27,24,11,17,15,25,19,21,27]

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

    endif


end subroutine

subroutine build_elements_baricenter(pc, ip, ie, znod)
  type(piece),               intent(inout) :: pc ! Pieze
  integer,                   intent(in)    :: ip ! Pieze number
  integer,                   intent(in)    :: ie ! Element group number
  real(real64),dimension(:,:),allocatable, intent(inout) :: znod(:,:)
  integer                                  :: i, tp, mphlnn, indx, maxindx
  integer, dimension(:,:), allocatable     :: nn
  real(real64),dimension(:), allocatable   :: val

    tp = pc%el(ie)%type
    mphlnn = mphtxt_get_lnn(tp)
    if(size(pc%el(ie)%nn,1) /= mphlnn) then
      if(allocated(nn)) deallocate(nn); allocate(nn(mphlnn,pc%el(nelg)%nel)); nn = 0
      nn(1:size(pc%el(ie)%nn,1),:) = pc%el(ie)%nn(:,:)
      call move_alloc(from=nn,  to=pc%el(nelg)%nn)
      deallocate(nn)
    endif

    if(allocated(val)) deallocate(val); allocate(val(pc%dim))

    do i = 1, pc%el(ie)%nel
      val = sum(pc%z(1, pc%el(ie)%mm(:,i))) / FEDB(tp)%lnv
      if (pc%el(ie)%nn(mphlnn,i) == 0) pc%el(ie)%nn(mphlnn,i) = size(znod,1)+1
      indx = pc%el(ie)%nn(mphlnn,i)
      if(indx < size(znod,2)) then; call set_col(znod, val, indx)
      else; call set_col(znod, val, indx); endif
      if(maxindx<max(indx,size(znod,2))) maxindx=indx
    end do

    call reduce(znod, maxindx, pc%el(ie)%nel)
    call info('Added centroids in piece '//trim(string(ip))//' element group '//trim(string(ie)))

end subroutine

subroutine build_midface_coordinates(pc, ip, all_P1, znod)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! POR HACER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(piece),               intent(in)  :: pc
integer,                   intent(in)  :: ip
logical,                   intent(out) :: all_P1
real(real64), allocatable, intent(out) :: znod(:,:)
integer :: ig, i, j, k

!determine whether all elements are Lagrange P1
all_P1 = .true.
do ig = 1, size(pc%el,1)
  if (.not. FEDB(pc%el(ig)%type)%nver_eq_nnod) then
    all_P1 = .false.
    exit
  end if
end do
if (.not. all_P1) then
  !construct znod
  call alloc(znod, pc%dim, pc%nnod)
  do ig = 1, size(pc%el,1)
    associate(elg => pc%el(ig), tp => pc%el(ig)%type)
      if (FEDB(tp)%nver_eq_nnod) then
        !************************************* Lagrange P1 **************************************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnv
            znod(:, elg%mm(i,k)) = pc%z(:, elg%mm(i,k))
          end do
        end do
      elseif (FEDB(tp)%lnn == FEDB(tp)%lnv + FEDB(tp)%lne) then
        !************************************* Lagrange P2 **************************************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnv
            znod(:, elg%nn(i,k)) = pc%z(:, elg%mm(i,k))
          end do
          do i = 1, FEDB(tp)%lne
            znod(:, elg%nn(i+FEDB(tp)%lnv,k)) = (pc%z(:,elg%mm(FEDB(tp)%edge(1,i),k)) + pc%z(:,elg%mm(FEDB(tp)%edge(2,i),k)))/2
          end do
        end do
      elseif (tp == check_fe(.false.,  3, 3,  3, 0) .or. &
              tp == check_fe(.false.,  6, 4,  6, 4)) then
        !***** Triangle, Raviart-Thomas (edge) OR Tetrahedron, Nedelec (edge) ******************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lne
            znod(:, elg%nn(i,k)) = ( pc%z(:, elg%mm(FEDB(tp)%edge(1,i),k)) + pc%z(:, elg%mm(FEDB(tp)%edge(2,i),k)) )/2
          end do
        end do
      elseif (tp == check_fe(.false.,  4, 4,  6, 4)) then
        !************************************* Tetrahedron, Raviart-Thomas (face) ***************
        do k = 1, elg%nel
          do i = 1, FEDB(tp)%lnf
            do j = 1, pc%dim
              znod(j, elg%nn(i,k)) = sum(pc%z(1, elg%mm(FEDB(tp)%face(:,i),k))) / FEDB(FEDB(tp)%f_type)%lnv
            end do
          end do
        end do
      else
        call info('(module_pmh/build_node_coordinates) build node coordinates for element type '//trim(string(FEDB(tp)%desc))//&
        &' is not implemented: piece '//trim(string(ip))//', group '//trim(string(ig))//'; variable znod is not created')
      end if
    end associate
  end do
end if
end subroutine



end module

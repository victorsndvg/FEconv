module module_utils_pf3_fcnv

!-----------------------------------------------------------------------
! Module to manage PF3 (Flux) files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Victor Sande, victor(dot)sande(at)usc(dot)es
! Last update: 19/03/2014
!
! PUBLIC PROCEDURES:
! pf3_assign_element_type: returns element type given the PF3 element descriptors
! pf32pmh_ordering:  translates the node ordering from PF3 to PMH
! pmh2pf3_ordering(el, tp, prevnod):  translates the node ordering from PMH to PF3
! pf3_get_element_desc1(tp): returns the first PF3 element descriptor
! pf3_get_element_desc2(tp): returns the second PF3 element descriptor
! pf3_get_element_desc3(tp): returns the third PF3 element descriptor
! extend_elgroup(v, d): extends the element groups
! linear_search(n, x, u): implements a linear search
! write_header_line(iu,line,ch,comm): write a line in the PF3 file
!-----------------------------------------------------------------------

use basicmod, only: insert, reduce
use module_pmh_fcnv

implicit none


contains


!-----------------------------------------------------------------------
! pf3_assign_element_type(desc1, desc2, desc3): returns element type given the PF3 element descriptors
!-----------------------------------------------------------------------
! desc1: integer,first element descriptor in PF3 files
! desc2: integer,second element descriptor in PF3 files
! desc3: integer,third element descriptor in PF3 files
!-----------------------------------------------------------------------

function pf3_assign_element_type(desc1, desc2, desc3) result(res)
  integer, intent(in)  :: desc1, desc2, desc3
  integer              :: res

! Not supported yet:
! -----------------
! line first order shell elements --> (2,202,23)
! line second order shell elements --> (2,303,24)
! triangle first order shell elements --> (3,207,19)
! triangle second order shell elements --> (3,307,20)
! rectangle first order shell elements --> (4,2202,21)
! rectangle second order shell elements --> (4,3302,22)

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
  elseif(desc1 == 4 .and. desc2 == 202 .and. desc3 == 7) then ! Quadrangle P1
    res = check_fe(.true., 4, 4, 4, 0)
    call info('Element type: Triangle Lagrange P1')
  elseif(desc1 == 4 .and. desc2 == 303 .and. desc3 == 8) then ! Quadrangle P2
    res = check_fe(.false., 8, 4, 4, 0)
    call info('Element type: Triangle Lagrange P2')
  elseif(desc1 == 5 .and. desc2 == 4 .and. desc3 == 10) then ! Tetrahedron  P1
    res = check_fe(.true., 4, 4, 6, 4)
    call info('Element type: Tetrahedron Lagrange P1')
  elseif(desc1 == 5 .and. desc2 == 15 .and. desc3 == 11) then ! Tetrahedron  P2
    res = check_fe(.false., 10, 4, 6, 4)
    call info('Element type: Tetrahedron Lagrange P2')
  elseif(desc1 == 6 .and. desc2 == 207 .and. desc3 == 12) then ! Wedge P1
    res = check_fe(.true., 6, 6, 9, 5)
    call info('Element type: Wedge Lagrange P1')
  elseif(desc1 == 6 .and. desc2 == 307 .and. desc3 == 13) then ! Wedge P2
    res = check_fe(.false., 15, 6, 9, 5)
    call info('Element type: Wedge Lagrange P2')
  elseif(desc1 == 7 .and. desc2 == 2202 .and. desc3 == 15) then ! Hexahedron  P1
    res = check_fe(.true., 8, 8, 12, 6)
    call info('Element type: Hexahedron Lagrange P1')
  elseif(desc1 == 7 .and. desc2 == 3303 .and. desc3 == 16) then ! Hexahedron  P2
    res = check_fe(.false., 20, 8, 12, 6)
    call info('Element type: Hexahedron Lagrange P2')
  elseif(desc1 == 8 .and. desc2 == 4202 .and. desc3 == 17) then ! Pyramid  P1
    res = check_fe(.true., 5, 5, 8, 5)
    call info('Element type: Pyramid Lagrange P1')
  elseif(desc1 == 8 .and. desc2 == 4203 .and. desc3 == 18) then ! Pyramid  P2
    res = check_fe(.false., 13, 8, 12, 6)
    call info('Element type: Pyramid Lagrange P2')
  else
    call info('Element type: Unknown element type # '&
      //string(desc1)//' '//string(desc2)//' '//string(desc3))
    res = 0
  endif

end function


!-----------------------------------------------------------------------
! pf32pmh_ordering(el, tp):  translates the node ordering from PF3 to PMH
!-----------------------------------------------------------------------
! el: array of connectivities from an element stored in nn
! tp: FE type identifier from FEDB
!-----------------------------------------------------------------------

subroutine pf32pmh_ordering(el, tp)

  integer, dimension(:), intent(inout) :: el
  integer, intent(in) :: tp
  integer :: aux
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_utils_pf3/node_ordering # Element type not supported')
    endif

    if (tp == check_fe(.true., 1, 1, 0, 0)) then     ! Nodes
        ! PMH and PF3 uses the same node ordering in nodes

    elseif (tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge Lagrange P1
        ! PMH and PF3 uses the same node ordering in edges lagrange P1

    elseif (tp == check_fe(.true., 3, 3, 3, 0)) then ! Triangle Lagrange P1
        ! PMH and PF3 uses the same node ordering in triangles lagrange P1

    elseif (tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle Lagrange P1
        ! PMH and PF3 don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = PF3[1,4,3,2]
        aux = el(2); el(2) = el(4); el(4) = aux

    elseif (tp == check_fe(.true., 4, 4, 6, 4)) then ! Tetrahedron Lagrange P1
        ! PMH and PF3 don't have the same node ordering in tetrahedrons lagrange P1
        ! PMH[1,2,3,4] = PF3[1,3,2,4]
        aux = el(2); el(2) = el(3); el(3) = aux

    elseif (tp == check_fe(.true., 6, 6, 9, 5)) then ! Wedge Lagrange P1
        ! PMH and PF3 don't have the same node ordering in wedges lagrange P1
        ! PMH[1,2,3,4,5,6] = PF3[1,3,2,4,6,5]
        aux = el(2); el(2) = el(3); el(3) = aux
        aux = el(5); el(5) = el(6); el(6) = aux

    elseif (tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron Lagrange P1
        ! PMH and PF3 don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = MPH[1,4,3,2,5,8,7,6]
        aux = el(2); el(2) = el(4); el(4) = aux
        aux = el(6); el(6) = el(8); el(8) = aux

    elseif (tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge Lagrange P2
        ! PMH and PF3 uses the same node ordering in edges lagrange P2

    elseif (tp == check_fe(.false., 6, 3, 3, 0)) then ! Triangle Lagrange P2
        ! PMH and PF3 have the same node ordering in triangles lagrange P2

    elseif (tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadragle Lagrange P2
        ! PMH and PF3 don't have the same node ordering in quadrangles lagrange P2
        ! PMH[1,2,3,4,5,6,7,8] = PF3[1,4,3,2,8,7,6,5]
        if (allocated(auxel)) deallocate(auxel)
        allocate(auxel(size(el,1)))
        auxel(:) = el(:)

        el(1) = auxel(1); el(2) = auxel(4); el(3) = auxel(3); el(4) = auxel(2)
        el(5) = auxel(8); el(6) = auxel(7); el(7) = auxel(6); el(8) = auxel(5)

        deallocate(auxel)

    elseif (tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron Lagrange P2
        ! PMH and PF3 don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4, 5,6,7,8, 9,10] = MPH[1,3,2,4, 6,8,5,7, 9,10]
        aux = el(2); el(2) = el(3); el(3) = aux
        aux = el(5); el(5) = el(6); el(6) = el(8)
        el(8) = el(7); el(7) = aux

    elseif (tp == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron Lagrange P2
        ! PMH and PF3 don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] =
        ! PF3[1,4,3,2,5,8,7,6,12,11,10,9,17,20,19,18,16,15,14,13]
        if (allocated(auxel)) deallocate(auxel)
        allocate(auxel(size(el,1)))
        auxel(:) = el(:)

        el(1) = auxel(1); el(2) = auxel(4); el(3) = auxel(3); el(4) = auxel(2)
        el(5) = auxel(5); el(6) = auxel(8); el(7) = auxel(7); el(8) = auxel(6)
        el(9) = auxel(12); el(10) = auxel(11); el(11) = auxel(10); el(12) = auxel(9)
        el(13) = auxel(17); el(14) = auxel(20); el(15) = auxel(19); el(16) = auxel(18)
        el(17) = auxel(16); el(18) = auxel(15); el(19) = auxel(14); el(20) = auxel(13)

        deallocate(auxel)


    endif

end subroutine



!-----------------------------------------------------------------------
! pf3_get_element_desc1(tp): returns the first PF3 element descriptor
!-----------------------------------------------------------------------
! tp: type of PMH element
!-----------------------------------------------------------------------

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
  elseif(tp == check_fe(.true., 4, 4, 4, 0) .or. tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadrangle P1-P2
    res = 4
  elseif(tp == check_fe(.true., 4, 4, 6, 4) .or. tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron P1-P2
    res = 5
  elseif(tp == check_fe(.true., 6, 6, 9, 5) .or. tp == check_fe(.false., 15, 6, 9, 5)) then ! Wedge P1-P2
    res = 6
  elseif(tp == check_fe(.true., 8, 8, 12, 6) .or. tp==check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron  P1-P2
    res = 7
  elseif(tp == check_fe(.true., 5, 5, 8, 5) .or. tp == check_fe(.false., 13, 5, 8, 5)) then ! Pyramid  P1-P2
    res = 8
  else
    call error('module_utils_pf3/pf3_get_element_desc1 # Unknown element type #'//trim(string(tp)))
  endif

end function


!-----------------------------------------------------------------------
! pf3_get_element_desc2(tp): returns the second PF3 element descriptor
!-----------------------------------------------------------------------
! tp: type of PMH element
!-----------------------------------------------------------------------

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
  elseif(tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle P1
    res = 202
  elseif(tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadrangle P2
    res = 303
  elseif(tp == check_fe(.true., 4, 4, 6, 4)) then  ! Tetrahedron P1
    res = 4
  elseif(tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron P2
    res = 15
  elseif(tp == check_fe(.true., 6, 6, 9, 5)) then ! Wedge P1
    res = 207
  elseif(tp == check_fe(.false., 15, 6, 9, 5)) then ! Wedge P2
    res = 307
  elseif(tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron  P1
    res = 2202
  elseif(tp==check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron  P2
    res = 3303
  elseif(tp == check_fe(.true., 5, 5, 8, 5)) then !Pyramid  P1
    res = 4202
  elseif(tp == check_fe(.false., 13, 5, 8, 5)) then ! Pyramid  P2
    res = 4203
  else
    call error('module_utils_pf3/pf3_get_element_desc2 # Unknown element type #'//trim(string(tp)))
  endif


end function


!-----------------------------------------------------------------------
! pf3_get_element_desc3(tp): returns the third PF3 element descriptor
!-----------------------------------------------------------------------
! tp: type of PMH element
!-----------------------------------------------------------------------

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
  elseif(tp == check_fe(.true., 3, 3, 3, 0)) then  ! Triangle P1
    res = 5
  elseif(tp == check_fe(.false., 6, 3, 3, 0)) then  ! Triangle P2
    res = 6
  elseif(tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle P1
    res = 7
  elseif(tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadrangle P2
    res = 8
  elseif(tp == check_fe(.true., 4, 4, 6, 4)) then  ! Tetrahedron P1
    res = 10
  elseif(tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron P2
    res = 11
  elseif(tp == check_fe(.true., 6, 6, 9, 5)) then ! Wedge P1
    res = 12
  elseif(tp == check_fe(.false., 15, 6, 9, 5)) then ! Wedge P2
    res = 13
  elseif(tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron  P1
    res = 15
  elseif(tp==check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron  P2
    res = 16
  elseif(tp == check_fe(.true., 5, 5, 8, 5)) then !Pyramid  P1
    res = 17
  elseif(tp == check_fe(.false., 13, 5, 8, 5)) then ! Pyramid  P2
    res = 8
  else
    call error('module_utils_pf3/pf3_get_element_desc3 # Unknown element type #'//trim(string(tp)))
  endif


end function


!-----------------------------------------------------------------------
! pmh2pf3_ordering(el, tp, prevnod):  translates the node ordering from PMH to PF3
!-----------------------------------------------------------------------
! el: array of connectivities from an element stored in nn
! tp: FE type identifier from FEDB
! prevnnod: number of nodes of the previous pieces
!-----------------------------------------------------------------------


function pmh2pf3_ordering(el, tp, prevnnod) result(auxel)

  integer, dimension(:), intent(in) :: el
  integer, intent(in) :: tp, prevnnod
  integer, dimension(:), allocatable :: auxel

    if (tp <= 0) then
      call error('module_utils_pf3/node_ordering # Element type not supported')
    endif

    if (allocated(auxel)) deallocate(auxel)
    allocate(auxel(size(el,1)))
    auxel(:) = el(:)

    if (tp == check_fe(.true., 1, 1, 0, 0)) then     ! Nodes
        ! PMH and PF3 uses the same node ordering in nodes

    elseif (tp == check_fe(.true., 2, 2, 1, 0)) then ! Edge Lagrange P1
        ! PMH and PF3 uses the same node ordering in edges lagrange P1

    elseif (tp == check_fe(.true., 3, 3, 3, 0)) then ! Triangle Lagrange P1
        ! PMH and PF3 uses the same node ordering in triangles lagrange P1

    elseif (tp == check_fe(.true., 4, 4, 4, 0)) then ! Quadrangle Lagrange P1
        ! PMH and PF3 don't have the same node ordering in quadrangles lagrange P1
        ! PMH[1,2,3,4] = PF3[1,4,3,2]
        auxel(2) = el(4); auxel(4) = el(2)

    elseif (tp == check_fe(.true., 4, 4, 6, 4)) then ! Tetrahedron Lagrange P1
        ! PMH and PF3 don't have the same node ordering in tetrahedrons lagrange P1
        ! PMH[1,2,3,4] = PF3[1,3,2,4]
        auxel(2) = el(3); auxel(3) = el(2)

    elseif (tp == check_fe(.true., 6, 6, 9, 5)) then ! Wedge Lagrange P1
        ! PMH and PF3 don't have the same node ordering in wedges lagrange P1
        ! PMH[1,2,3,4,5,6] = PF3[1,3,2,4,6,5]
        auxel(2) = el(3); auxel(3) = el(2)
        auxel(5) = el(6); auxel(6) = el(5)

    elseif (tp == check_fe(.true., 8, 8, 12, 6)) then ! Hexahedron Lagrange P1
        ! PMH and PF3 don't have the same node ordering in hexahedrons lagrange P1
        ! PMH[1,2,3,4,5,6,7,8] = PF3[1,4,3,2,5,8,7,6]
        auxel(2) = el(4); auxel(4) = el(2)
        auxel(6) = el(8); auxel(8) = el(6)

    elseif (tp == check_fe(.false., 3, 2, 1, 0)) then ! Edge Lagrange P2
        ! PMH and PF3 uses the same node ordering in edges lagrange P2

    elseif (tp == check_fe(.false., 6, 3, 3, 0)) then ! Triangle Lagrange P2
        ! PMH and PF3 have the same node ordering in triangles lagrange P2

    elseif (tp == check_fe(.false., 8, 4, 4, 0)) then ! Quadragle Lagrange P2
        ! PMH and PF3 don't have the same node ordering in quadrangles lagrange P2
        ! PMH[1,2,3,4,5,6,7,8] = PF3[1,4,3,2,8,7,6,5]
        auxel(2) = el(4); auxel(4) = el(2)
        auxel(5) = el(8); auxel(6) = el(7); auxel(7) = el(6); auxel(8) = el(5)


    elseif (tp == check_fe(.false., 10, 4, 6, 4)) then ! Tetrahedron Lagrange P2
        ! PMH and PF3 don't have the same node ordering in tetrahedrons lagrange P2
        ! PMH[1,2,3,4, 5,6,7,8, 9,10] = PF3[1,3,2,4, 6,8,5,7, 9,10]
        auxel(2) = el(3); auxel(3) = el(2)
        auxel(5) = el(7); auxel(6) = el(5)
        auxel(7) = el(8); auxel(8) = el(6)

    elseif (tp == check_fe(.false., 20, 8, 12, 6)) then ! Hexahedron Lagrange P2
        ! PMH and PF3 don't have the same node ordering in hexahedrons lagrange P2
        ! PMH[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] =
        ! PF3[1,4,3,2,5,8,7,6,12,11,10,9,17,20,19,18,16,15,14,13]
        auxel(1) = el(1); auxel(2) = el(4); auxel(3) = el(3); auxel(4) = el(2)
        auxel(5) = el(5); auxel(6) = el(8); auxel(7) = el(7); auxel(8) = el(6)
        auxel(9) = el(12); auxel(10) = el(11); auxel(11) = el(10); auxel(12) = el(9)
        auxel(13) = el(20); auxel(14) = el(19); auxel(15) = el(18); auxel(16) = el(17)
        auxel(17) = el(13); auxel(18) = el(16); auxel(19) = el(15); auxel(20) = el(14)

    endif

    auxel(:) = auxel(:) + prevnnod

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
  character(maxpath) :: cad

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


!-----------------------------------------------------------------------
! linear_search(n, x, u): If 'u' is in 'x' then its index is returned.
! Otherwise −ind is returned. Negative values on 'x' will be considered as empty positions.
!-----------------------------------------------------------------------
! n: size of the sequence
! x: integer array where to search
! u: value to search
!-----------------------------------------------------------------------

function linear_search(n, x, u) result(pos)
integer,               intent(in) :: n !components of seq
integer, dimension(:), intent(in) :: x !array, where to search
integer,               intent(in) :: u !value to search
integer :: pos, i, j

pos = -1
if (size(x,1) <= 0 .or. n <= 0) return

i=1; j=n
do i = 1, size(x,1)
  if(x(i) < 0) then
    pos = -i
    return
  elseif(x(i) == u) then
    pos = i
    return
  end if
end do

end function


!-----------------------------------------------------------------------
! write_header_line(iu,line,ch,comm): write a line in the PF3 file
!-----------------------------------------------------------------------
! iu:   unit number of the PF3 file
! line: text included in one line
! ch:   comments character (Optional)
! comm: commentary (Optional)
!-----------------------------------------------------------------------

subroutine write_header_line(iu,line,num)
  integer, intent(in) :: iu ! File unit number
  character(len=*), intent(in) :: line ! String
  integer, optional, intent(in) :: num ! String: Comment character
  integer :: ios

  if(present(num)) then
    write(unit=iu, fmt='(1I8,a)', iostat = ios) num, '           '//trim(line)
    if (ios /= 0) call error('write_utils_pf3/header_line, #'//trim(string(ios)))
  else
    write(unit=iu, fmt='(a)', iostat = ios)  ' '//trim(line)
    if (ios /= 0) call error('write_utils_pf3/header_line, #'//trim(string(ios)))
  endif



end subroutine

end module

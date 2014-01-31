module module_fe_database_pmh
!-----------------------------------------------------------------------
! Module to manage the finite elements database for PMH format
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   check_fe: given the finite element parameters, return the index of its position in the database
!
! REMARKS:
!   To check the type of a mesh written in another format:
!      it = check_fe(...)
!      if (it == 0) call error('...: element type not implemented for PMH format')
!      pmesh%pc(...)%el(...)%type = fedb%type(it)
!      ...
!
!   To check whether a PMH mesh can be written in another format:
!      use module_alloc
!      (The integer array valid_fe is created using check_fe_param)
!      if (find_first(valid_fe, pmesh%pc(...)%el(...)%type) == 0) call error(...: element type ... cannot be written')
!-----------------------------------------------------------------------
implicit none

!Constants
!EDGE_TRIA(i,j), vertex #i of edge #j of a triangle
integer, parameter :: EDGE_TRIA(2,12)  = reshape([1,2, 2,3, 3,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
!EDGE_QUAD(i,j), vertex #i of edge #j of a quadrangle
integer, parameter :: EDGE_QUAD(2,12)  = reshape([1,2, 2,3, 3,4, 4,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
!EDGE_TETRA(i,j), vertex #i of edge #j of a tetrahedron
integer, parameter :: EDGE_TETRA(2,12) = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
!EDGE_HEXA(i,j), vertex #i of edge #j of a hexahedron
integer, parameter :: EDGE_HEXA(2,12)  = reshape([1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5], [2,12])
!FACE_TETRA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_TETRA(4,6) = reshape([1,3,2,0, 1,4,3,0, 1,2,4,0, 2,3,4,0, 0,0,0,0, 0,0,0,0], [4,6])
!FACE_HEXA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_HEXA(4,6)  = reshape([1,4,3,2, 1,5,8,4, 1,2,6,5, 5,6,7,8, 2,3,7,6, 3,4,8,7], [4,6])

!Types
type :: fe_db_pmh
  character(len=46) :: desc       = ' '       !description
  logical           :: nver_eq_nnod = .false. !is the total number of vertices equal to the total number of nodes?
  integer           :: tdim       = 0         !topological dimension (0 for nodes, 1 for edges, etc.)
  integer           :: lnn        = 0         !local number of nodes
  integer           :: lnv        = 0         !local number of vertices
  integer           :: lne        = 0         !local number of edges
  integer           :: lnf        = 0         !local number of faces
  integer           :: edge(2,12) = 0         !local numbering of vertices in each edge (lne <= 12)
  integer           :: lnv_f      = 0         !local number of vertices per face
  integer           :: face(4,6)  = 0         !local numbering of vertices in each edge (lnf <= 4, lnv_f <= 6)   
end type

!Constants                                     char(46)   td   n  v   e  f         ed  vf          fa
type(fe_db_pmh), parameter :: FEDB(17) = [ &
  fe_db_pmh('Node                              ', .true.,  0,  1, 1,  0, 0,         0,  0,          0), &
  fe_db_pmh('Edge, Lagrange P1                 ', .true.,  1,  2, 2,  1, 0,         0,  0,          0), &
  fe_db_pmh('Edge, Lagrange P2                 ', .false., 1,  3, 2,  1, 0,         0,  0,          0), &
  fe_db_pmh('Triangle, Lagrange P1             ', .true.,  2,  3, 3,  3, 0, EDGE_TRIA,  0,          0), &
  fe_db_pmh('Triangle, Lagrange P2             ', .false., 2,  6, 3,  3, 0, EDGE_TRIA,  0,          0), &
  fe_db_pmh('Triangle, Raviart-Thomas (edge)   ', .false., 2,  3, 3,  3, 0, EDGE_TRIA,  0,          0), &
  fe_db_pmh('Quadrangle, Lagrange P1           ', .true.,  2,  4, 4,  4, 0, EDGE_QUAD,  0,          0), &
  fe_db_pmh('Quadrangle, Lagrange P2           ', .false., 2,  8, 4,  4, 0, EDGE_QUAD,  0,          0), &
  fe_db_pmh('Quadrangle, Raviart-Thomas (edge) ', .false., 2,  4, 4,  4, 0, EDGE_QUAD,  0,          0), &
  fe_db_pmh('Tetrahedron, Lagrange P1          ', .true.,  3,  4, 4,  6, 4, EDGE_TETRA, 3, FACE_TETRA), &
  fe_db_pmh('Tetrahedron, Lagrange P2          ', .false., 3, 10, 4,  6, 4, EDGE_TETRA, 3, FACE_TETRA), &
  fe_db_pmh('Tetrahedron, Raviart-Thomas (face)', .false., 3,  4, 4,  6, 4, EDGE_TETRA, 3, FACE_TETRA), &
  fe_db_pmh('Tetrahedron, Nedelec (edge)       ', .false., 3,  6, 4,  6, 4, EDGE_TETRA, 3, FACE_TETRA), &
  fe_db_pmh('Hexahedron, Lagrange P1           ', .true.,  3,  8, 8, 12, 6, EDGE_HEXA,  3, FACE_HEXA),  &
  fe_db_pmh('Hexahedron, Lagrange P2           ', .false., 3, 16, 8, 12, 6, EDGE_HEXA,  4, FACE_HEXA),  &
  fe_db_pmh('Hexahedron, Raviart-Thomas (face) ', .false., 3,  6, 8, 12, 6, EDGE_HEXA,  4, FACE_HEXA),  &
  fe_db_pmh('Hexahedron, Nedelec (edge)        ', .false., 3, 12, 8, 12, 6, EDGE_HEXA,  4, FACE_HEXA)]

contains

!-----------------------------------------------------------------------
! check_fe: given the finite element parameters, return the index of its position in the database
!-----------------------------------------------------------------------
function check_fe(nnod, nver, lnn, lnv, lne, lnf) result(res)
integer, intent(in) :: nnod, nver, lnn, lnv, lne, lnf
integer :: res, i

res = 0
do i = 1, size(FEDB,1)
  if (((nver==nnod).eqv.FEDB(i)%nver_eq_nnod) .and. lnn==FEDB(i)%lnn .and. lnv==FEDB(i)%lnv .and. &
                                                    lne==FEDB(i)%lne .and. lnf==FEDB(i)%lnf) then
    res = i
    return
  end if
enddo
end function

end module

!Innecesary/old constants:
!Element types
!integer, parameter :: NODE2 = 1, NODE3 = 2, &                                   !nodes
!ED2_P1 =  3, ED3_P1 =  4, ED2_P2 =  5, ED3_P2 =  6, &                           !edges 
!TR2_P1 =  7, TR3_P1 =  8, TR2_P2 =  8, TR3_P2 = 10, TR2_RT = 11, TR3_RT = 12, & !triangles
!QU2_P1 = 13, QU3_P1 = 14, QU2_P2 = 15, QU3_P2 = 16, QU2_RT = 17, QU3_RT = 18, & !quadrangles
!TET_P1 = 19, TET_P2 = 20, TET_RT = 21, TET_ND = 22, &                           !tetrahedra
!HEX_P1 = 23, HEX_P2 = 24, HEX_RT = 25, HEX_ND = 26                              !hexahedra
!EDGE_WEDGE(i,j), vertex #i of edge #j of a wedge
!integer, parameter :: EDGE_WEDGE(2,9) = reshape([1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4], [2,9])
!FACE_WEDGE(i,j), vertex #i of face #j of a wedge; note that faces 1 and 4 have only three valid vertices
!integer, parameter :: FACE_WEDGE(4,5) = reshape([1,3,2,0, 1,4,6,3, 1,2,5,4, 4,5,6,0, 2,3,6,5], [4,5]) 


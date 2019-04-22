module module_fe_database_pmh_fcnv
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
!ED_????(i,j), vertex #i of edge #j of element ????
integer, parameter :: ED_EDGE(2,12) = reshape([1,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
integer, parameter :: ED_TRIA(2,12) = reshape([1,2, 2,3, 3,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
integer, parameter :: ED_QUAD(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
integer, parameter :: ED_TETR(2,12) = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0], [2,12])
integer, parameter :: ED_HEXA(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5], [2,12])
integer, parameter :: ED_WEDG(2,12) = reshape([1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4, 0,0, 0,0, 0,0], [2,12])
!ED_????_P2(i,j), node #i of edge #j of element Lagrange P2 ????
integer, parameter :: ED_TRIA_P2(3,12) = reshape([1,2,4,  2,3,5,  3,1,6,  0,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0, &
 0,0,0,  0,0,0], [3,12])
integer, parameter :: ED_QUAD_P2(3,12) = reshape([1,2,5,  2,3,6,  3,4,7,  4,1,8,  0,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0, &
 0,0,0,  0,0,0], [3,12])
integer, parameter :: ED_TETR_P2(3,12) = reshape([1,2,5,  2,3,6,  3,1,7,  1,4,8,  2,4,9, 3,4,10,  0,0,0,  0,0,0,  0,0,0,  0,0,0, &
 0,0,0,  0,0,0], [3,12])
integer, parameter :: ED_HEXA_P2(3,12) = reshape([1,2,9, 2,3,10, 3,4,11, 4,1,12, 1,5,13, 2,6,14, 3,7,15, 4,8,16, 5,6,17, 6,7,18, &
7,8,19, 8,5,20], [3,12])
!FA_????(i,j), vertex #i of face #j of element ????
integer, parameter :: FA_TETR(4,6) = reshape([1,3,2,0, 1,4,3,0, 1,2,4,0, 2,3,4,0, 0,0,0,0, 0,0,0,0], [4,6])
integer, parameter :: FA_HEXA(4,6) = reshape([1,4,3,2, 1,5,8,4, 1,2,6,5, 5,6,7,8, 2,3,7,6, 3,4,8,7], [4,6])
integer, parameter :: FA_WEDG(4,6) = reshape([1,3,2,0, 1,4,6,3, 1,2,5,4, 4,5,6,0, 2,3,6,5, 0,0,0,0], [4,6]) 
integer, parameter :: VF_WEDG(5)   = [3, 4, 4, 3, 4] !vertices per face in a wedge
!FA_????_P2(i,j), node #i of face #j of element Lagrange P2 ????
integer, parameter :: FA_TETR_P2(8,6) = reshape([1,3,2,7,6,5,0,0, 1,4,3,8,10,7,0,0, 1,2,4,5,9,8,0,0, 2,3,4,6,10,9,0,0, &
0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0], [8,6]) 
integer, parameter :: FA_HEXA_P2(8,6) = reshape([1,4,3,2,12,11,10,9, 1,5,8,4,13,20,16,12, 1,2,6,5,9,14,17,13, 5,6,7,8,17,18,19,20, &
2,3,7,6,10,15,18,14, 3,4,8,7,11,16,19,15], [8,6]) 
 
!Types
type :: fe_db_pmh
  character(len=46) :: desc         = ' '     !description
  logical           :: nver_eq_nnod = .false. !is the total number of vertices equal to the total number of nodes?
  integer           :: tdim         = 0       !topological dimension (0 for nodes, 1 for edges, etc.)
  integer           :: lnn          = 0       !local number of nodes
  integer           :: lnv          = 0       !local number of vertices
  integer           :: lne          = 0       !local number of edges
  integer           :: lnf          = 0       !local number of faces
  integer           :: v_type       = 0       !f.e. type induced in vertices
  integer           :: e_type       = 0       !f.e. type induced in edges (for P1, RT, ND, edges P1 are induced) 
  integer           :: edge( 2,12)  = 0       !local numbering of vertices in each edge (lnv/edge  = 2, lne <= 12)
  integer           :: nedge(3,12)  = 0       !local numbering of nodes    in each edge (lnn/edge  = 3, lne <= 12)
  integer           :: f_type       = 0       !f.e. type induced in faces (for P1, RT, ND, faces P1 are induced)
  integer           :: face( 4,6)   = 0       !local numbering of vertices in each edge (lnv/face <= 4, lnf <= 6)
  integer           :: nface(8,6)   = 0       !local numbering of nodes    in each edge (lnn/face <= 8, lnf <= 6)
end type

!Constants             char(46)                         td   n  v   e  f vt et       ev          en ft       fv
type(fe_db_pmh), parameter :: FEDB(16) = [ &
fe_db_pmh('Vertex                            ', .true.,  0,  1, 1,  0, 0, 0, 0,       0,          0, 0,       0,          0), & ! 1
fe_db_pmh('Edge, Lagrange P1                 ', .true.,  1,  2, 2,  1, 0, 1, 0, ED_EDGE,          0, 0,       0,          0), & ! 2 
fe_db_pmh('Edge, Lagrange P2                 ', .false., 1,  3, 2,  1, 0, 1, 0, ED_EDGE,          0, 0,       0,          0), & ! 3
fe_db_pmh('Triangle, Lagrange P1             ', .true.,  2,  3, 3,  3, 0, 1, 2, ED_TRIA,          0, 0,       0,          0), & ! 4
fe_db_pmh('Triangle, Lagrange P2             ', .false., 2,  6, 3,  3, 0, 1, 3, ED_TRIA, ED_TRIA_P2, 0,       0,          0), & ! 5
fe_db_pmh('Triangle, Raviart-Thomas (edge)   ', .false., 2,  3, 3,  3, 0, 1, 2, ED_TRIA,          0, 0,       0,          0), & ! 6
fe_db_pmh('Quadrangle, Lagrange P1           ', .true.,  2,  4, 4,  4, 0, 1, 2, ED_QUAD,          0, 0,       0,          0), & ! 7
fe_db_pmh('Quadrangle, Lagrange P2           ', .false., 2,  8, 4,  4, 0, 1, 3, ED_QUAD, ED_QUAD_P2, 0,       0,          0), & ! 8
fe_db_pmh('Tetrahedron, Lagrange P1          ', .true.,  3,  4, 4,  6, 4, 1, 2, ED_TETR,          0, 4, FA_TETR,          0), & ! 9
fe_db_pmh('Tetrahedron, Lagrange P2          ', .false., 3, 10, 4,  6, 4, 1, 3, ED_TETR, ED_TETR_P2, 5, FA_TETR, FA_TETR_P2), & !10
fe_db_pmh('Tetrahedron, Raviart-Thomas (face)', .false., 3,  4, 4,  6, 4, 1, 2, ED_TETR,          0, 4, FA_TETR,          0), & !11
fe_db_pmh('Tetrahedron, Nedelec (edge)       ', .false., 3,  6, 4,  6, 4, 1, 2, ED_TETR,          0, 4, FA_TETR,          0), & !12
fe_db_pmh('Tetrahedron, Nedelec 2 (edge)     ', .false., 3, 20, 4,  6, 4, 1, 2, ED_TETR,          0, 4, FA_TETR,          0), & !13
fe_db_pmh('Hexahedron, Lagrange P1           ', .true.,  3,  8, 8, 12, 6, 1, 2, ED_HEXA,          0, 7, FA_HEXA,          0), & !14
fe_db_pmh('Hexahedron, Lagrange P2           ', .false., 3, 20, 8, 12, 6, 1, 3, ED_HEXA, ED_HEXA_P2, 8, FA_HEXA, FA_HEXA_P2), & !15
fe_db_pmh('Wedge, Lagrange P1                ', .true.,  3,  6, 6,  9, 5, 1, 2, ED_WEDG,          0, 0, FA_WEDG, 0)]!vf variable 16

contains

!-----------------------------------------------------------------------
! check_fe: given the finite element parameters, return the index of its position in the database
!-----------------------------------------------------------------------
function check_fe(nver_eq_nnod, lnn, lnv, lne, lnf) result(res)
logical, intent(in) :: nver_eq_nnod
integer, intent(in) :: lnn, lnv, lne, lnf
integer :: res, i

res = 0
do i = 1, size(FEDB,1)
  if ((nver_eq_nnod.eqv.FEDB(i)%nver_eq_nnod) .and. lnn==FEDB(i)%lnn .and. lnv==FEDB(i)%lnv .and. &
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


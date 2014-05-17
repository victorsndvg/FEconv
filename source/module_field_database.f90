module module_field_database
!-----------------------------------------------------------------------
! Module to manage the field database
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 16/05/2014
!
! PUBLIC PROCEDURES:
! id_mesh_ext: given the mesh extension, return the index of its position in the database
!-----------------------------------------------------------------------
implicit none
 
!Types
type :: field_db
  character(len=10) :: mesh_ext         = ' '     !mesh  extension
  character(len=10) :: field_ext        = ' '     !field extension
  logical           :: is_field_outside = .false. !is field outside mesh?
end type

type(field_db), parameter :: FLDB(2) = [       &
field_db('mfm       ', 'mff       ', .true.), &
field_db('mum       ', 'muf       ', .true.)]

contains

!-----------------------------------------------------------------------
! id_mesh_ext: given the mesh extension, return the index of its position in the database
!-----------------------------------------------------------------------
function id_mesh_ext(mesh_ext) result(res)
character(*), intent(in) :: mesh_ext
integer :: res, i

res = 0
do i = 1, size(FLDB,1)
  if (mesh_ext == FLDB(i)%mesh_ext) then
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


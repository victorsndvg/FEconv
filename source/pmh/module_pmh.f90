module module_pmh
!-----------------------------------------------------------------------
! Module to manage piecewise meshes (PMH)
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: 19/01/2014
!
! PUBLIC PROCEDURES:
!   none
! REMARKS:
!   A mesh is divided into pieces 
!   Each piece has common vertices/nodes and it can contain one or several groups of elements
!   Each group of elements belong to a type of element defined in module_eltype
!   If either znod or z is known, the other variable can be constructed using el()%type, el()%nn, , el()%mm
!   For P1 elements, el()%nn and el()%mm are redundant
!
! PASOS A DAR:
!   Definición de pmh
!   Creacion de funciones para procesar las coordenadas, el nn, etc
!   Opcionalmente, creación de funciones para procesar refs, orientación, etc.
!   Conversor a/de mfm
!   Conversor a vtu
!   Uso en comsol
!-----------------------------------------------------------------------
use module_fe_database_pmh
use module_ALLOC
implicit none

!Constants
!EDGE_TRIA(i,j), vertex #i of edge #j of a triangle
integer, parameter :: EDGE_TRIA(2,3) = reshape([1,2, 2,3, 3,1], [2,3])

!EDGE_QUAD(i,j), vertex #i of edge #j of a quadrangle
integer, parameter :: EDGE_QUAD(2,4) = reshape([1,2, 2,3, 3,4, 4,1], [2,4])

!EDGE_TETRA(i,j), vertex #i of edge #j of a tetrahedron
integer, parameter :: EDGE_TETRA(2,6) = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2,6])

!FACE_TETRA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_TETRA(3,4) = reshape([1,3,2, 1,4,3, 1,2,4, 2,3,4], [3,4])

!EDGE_HEXA(i,j), vertex #i of edge #j of a hexahedron
integer, parameter :: EDGE_HEXA(2,12) = reshape([1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5], [2,12])

!FACE_HEXA(i,j), vertex #i of face #j of a tetrahedron
integer, parameter :: FACE_HEXA(4,6) = reshape([1,4,3,2, 1,5,8,4, 1,2,6,5, 5,6,7,8, 2,3,7,6, 3,4,8,7], [4,6])

!EDGE_WEDGE(i,j), vertex #i of edge #j of a wedge
integer, parameter :: EDGE_WEDGE(2,9) = reshape([1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4], [2,9])

!FACE_WEDGE(i,j), vertex #i of face #j of a wedge
integer, parameter :: FACE_WEDGE(4,5) = reshape([1,3,2,0, 1,4,6,3, 1,2,5,4, 4,5,6,0, 2,3,6,5], [4,5]) !note that faces 1 and 4 have only three valid vertices

!Types
type elgroup
  integer              :: type = 0 !element type (one of those defined in module_eltype)
  integer              :: nel  = 0 !total number of elements
  integer              :: lnn  = 0 !local number of nodes per element
  integer              :: lnv  = 0 !local number of vertices per element
  integer, allocatable :: nn(:,:)  !global numbering of nodes
  integer, allocatable :: mm(:,:)  !global numbering of vertices
  integer, allocatable :: ref(:)   !reference numbering
end type

type piece
  integer                    :: dim  = 0  ! space dimension of the node/vertex coordinates
  integer                    :: nnod = 0  ! total number of nodes
  integer                    :: nver = 0  ! total number of vertices
  real(real64),  allocatable :: znod(:,:) ! node coordinates
  real(real64),  allocatable :: z(:,:)    ! vertex coordinates
  type(elgroup), allocatable :: el(:)     ! element groups  
end type

type pmh_mesh
  type(piece), allocatable :: pc(:) ! pieces that compose the mesh
end type  

end module

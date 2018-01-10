module module_FE_DB_fcnv
!-----------------------------------------------------------------------
! Module to manage the fe descriptor id's database
! Last update: 04/04/2010
!
! MODULE CONSTANTS (public):
!   FE_DB: fe database
!
! OUTPUT PROCEDURES:
!   find_descriptorID: checks whether the given descriptor is in the database
!-----------------------------------------------------------------------
use basicmod, only: maxpath
implicit none

!Constants
!ND_????_P2(i,j), Salome node order corresponding to a PMH nn from a Lagrange P2 element
integer, parameter :: ND_EDGE_P2(20) = [1,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
integer, parameter :: ND_TRIA_P2(20) = [1,4,2,5,3,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
integer, parameter :: ND_QUAD_P2(20) = [1,5,2,6,3,7,4,8,0,0,0,0,0,0,0,0,0,0,0,0]
integer, parameter :: ND_TETR_P2(20) = [1,5,2,6,3,7,8,9,10,4,0,0,0,0,0,0,0,0,0,0]
integer, parameter :: ND_HEXA_P2(20) = [1,9,2,10,3,11,4,12,13,14,15,16,5,17,6,18,7,19,8,20]


!Types
type fe_db_type
  character(len=MAXPATH) :: desc = ' '       !description
  integer                :: DIM  = 0         !space dimension
  integer                :: LNN  = 0         !local number of nodes
  integer                :: LNV  = 0         !local number of vertices
  integer                :: LNE  = 0         !local number of edges
  integer                :: LNF  = 0         !local number of faces
  logical                :: Beam = .false.   !ef is beam
  integer                :: nn_order(20) = 0 !Salome P2 node order
end type

!Constants
type(fe_db_type), public, dimension(232), save :: FE_DB

contains

!-----------------------------------------------------------------------
! init: give initial values to FE_DB
!-----------------------------------------------------------------------
subroutine FE_DB_init()

FE_DB( 11) = fe_db_type('Rod',                                 1, 2,2, 1,0,.true.,          0)
FE_DB( 21) = fe_db_type('Linear beam',                         1, 2,2, 1,0,.true.,          0)
FE_DB( 22) = fe_db_type('Tapered beam',                        1, 3,2, 1,0,.true., ND_EDGE_P2)

FE_DB( 24) = fe_db_type('Parabolic beam',                      1, 3,2, 1,0,.true., ND_EDGE_P2)

FE_DB( 41) = fe_db_type('Plane Stress Linear Triangle',        2, 3,3, 3,0,.false.,         0)
FE_DB( 42) = fe_db_type('Plane Stress Parabolic Triangle',     2, 6,3, 3,0,.false.,ND_TRIA_P2)
FE_DB( 44) = fe_db_type('Plane Stress Linear Quadrilateral',   2, 4,4, 4,0,.false.,         0)
FE_DB( 45) = fe_db_type('Plane Stress Parabolic Quadrilateral',2, 8,4, 4,0,.false.,ND_QUAD_P2)

FE_DB( 81) = fe_db_type('Axisymetric Solid Linear Triangle',   2, 3,3, 3,0,.false.,         0)
FE_DB( 82) = fe_db_type('Axisymetric Solid Parabolic Triangle',3, 6,3, 3,0,.false.,ND_TRIA_P2)
FE_DB( 84) = fe_db_type('Axisymetric Solid Linear Quad',       3, 4,4, 4,0,.false.,         0)
FE_DB( 85) = fe_db_type('Axisymetric Solid Parabolic Quad',    3, 8,4, 4,0,.false.,ND_QUAD_P2)

FE_DB( 91) = fe_db_type('Thin Shell Linear Triangle',          3, 3,3, 3,0,.false.,         0)
FE_DB( 92) = fe_db_type('Thin Shell Parabolic Triangle',       3, 6,3, 3,0,.false.,ND_TRIA_P2)
FE_DB( 94) = fe_db_type('Thin Shell Linear Quadrilateral',     3, 4,4, 4,0,.false.,         0)
FE_DB( 95) = fe_db_type('Thin Shell Parabolic Quadrilateral',  3, 8,4, 4,0,.false.,ND_QUAD_P2)

FE_DB(111) = fe_db_type('Solid Linear Tetrahedron',            3, 4,4, 6,4,.false.,         0)
FE_DB(112) = fe_db_type('Solid Linear Wedge',                  3, 6,6, 9,5,.false.,         0)
FE_DB(115) = fe_db_type('Solid Linear Brick',                  3, 8,8,12,6,.false.,         0)
FE_DB(116) = fe_db_type('Solid Parabolic Brick',               3,20,8,12,6,.false.,ND_HEXA_P2)
FE_DB(118) = fe_db_type('Solid Parabolic Tetrahedron',         3,10,4, 6,4,.false.,ND_TETR_P2)

FE_DB(122) = fe_db_type('Rigid Element',                       2, 4,4, 4,0,.false.,         0)
end subroutine

!-----------------------------------------------------------------------
! check_unv_fe: given the finite element parameters, return the index of its position in the database
!-----------------------------------------------------------------------
function check_unv_fe(tdim, lnn, lnv, lne, lnf) result(res)
integer, intent(in) :: tdim, lnn, lnv, lne, lnf
integer :: res, i

res = 0
do i = 1, size(FE_DB,1)
  if (tdim == 1) then
     if (lnn == FE_DB(i)%LNN .and. lnv == FE_DB(i)%LNV) then
       res = i
       return
     end if
   elseif (tdim == 2) then
     if (tdim <= FE_DB(i)%DIM .and. lnn == FE_DB(i)%LNN .and. lnv == FE_DB(i)%LNV .and. lne == FE_DB(i)%LNE) then
       res = i
       return
     end if
   elseif (tdim == 3) then
     if (tdim == FE_DB(i)%DIM .and. lnn == FE_DB(i)%LNN .and. lnv == FE_DB(i)%LNV .and. lne == FE_DB(i)%LNE .and. &
     lnf == FE_DB(i)%LNF) then
       res = i
       return
     end if
   end if
enddo
end function

end module

!-----------------------------
! FE Descriptor Id definitions
!-----------------------------
!
!   11  Rod
!   21  Linear beam
!   22  Tapered beam
!   23  Curved beam
!   24  Parabolic beam
!   31  Straight pipe
!   32  Curved pipe
!   41  Plane Stress Linear Triangle
!   42  Plane Stress Parabolic Triangle
!   43  Plane Stress Cubic Triangle
!   44  Plane Stress Linear Quadrilateral
!   45  Plane Stress Parabolic Quadrilateral
!   46  Plane Strain Cubic Quadrilateral
!   51  Plane Strain Linear Triangle
!   52  Plane Strain Parabolic Triangle
!   53  Plane Strain Cubic Triangle
!   54  Plane Strain Linear Quadrilateral
!   55  Plane Strain Parabolic Quadrilateral
!   56  Plane Strain Cubic Quadrilateral
!   61  Plate Linear Triangle
!   62  Plate Parabolic Triangle
!   63  Plate Cubic Triangle
!   64  Plate Linear Quadrilateral
!   65  Plate Parabolic Quadrilateral
!   66  Plate Cubic Quadrilateral
!   71  Membrane Linear Quadrilateral
!   72  Membrane Parabolic Triangle
!   73  Membrane Cubic Triangle
!   74  Membrane Linear Triangle
!   75  Membrane Parabolic Quadrilateral
!   76  Membrane Cubic Quadrilateral
!   81  Axisymetric Solid Linear Triangle
!   82  Axisymetric Solid Parabolic Triangle
!   84  Axisymetric Solid Linear Quadrilateral
!   85  Axisymetric Solid Parabolic Quadrilateral
!   91  Thin Shell Linear Triangle
!   92  Thin Shell Parabolic Triangle
!   93  Thin Shell Cubic Triangle
!   94  Thin Shell Linear Quadrilateral
!   95  Thin Shell Parabolic Quadrilateral
!   96  Thin Shell Cubic Quadrilateral
!   101 Thick Shell Linear Wedge
!   102 Thick Shell Parabolic Wedge
!   103 Thick Shell Cubic Wedge
!   104 Thick Shell Linear Brick
!   105 Thick Shell Parabolic Brick
!   106 Thick Shell Cubic Brick
!   111 Solid Linear Tetrahedron
!   112 Solid Linear Wedge
!   113 Solid Parabolic Wedge
!   114 Solid Cubic Wedge
!   115 Solid Linear Brick
!   116 Solid Parabolic Brick
!   117 Solid Cubic Brick
!   118 Solid Parabolic Tetrahedron
!   121 Rigid Bar
!   122 Rigid Element
!   136 Node To Node Translational Spring
!   137 Node To Node Rotational Spring
!   138 Node To Ground Translational Spring
!   139 Node To Ground Rotational Spring
!   141 Node To Node Damper
!   142 Node To Gound Damper
!   151 Node To Node Gap
!   152 Node To Ground Gap
!   161 Lumped Mass
!   171 Axisymetric Linear Shell
!   172 Axisymetric Parabolic Shell
!   181 Constraint
!   191 Plastic Cold Runner
!   192 Plastic Hot Runner
!   193 Plastic Water Line
!   194 Plastic Fountain
!   195 Plastic Baffle
!   196 Plastic Rod Heater
!   201 Linear node-to-node interface
!   202 Linear edge-to-edge interface
!   203 Parabolic edge-to-edge interface
!   204 Linear face-to-face interface
!   208 Parabolic face-to-face interface
!   212 Linear axisymmetric interface
!   213 Parabolic axisymmetric interface
!   221 Linear rigid surface
!   222 Parabolic rigin surface
!   231 Axisymetric linear rigid surface
!   232 Axisymentric parabolic rigid surface

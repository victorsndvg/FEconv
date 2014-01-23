module module_FE_DB
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
use module_OS_DEPENDANT, only: MAXPATH
implicit none

!Types
type fe_db_type
  character(len=MAXPATH) :: desc = ' '     !description
  integer                :: DIM  = 0       !space dimension
  integer                :: LNN  = 0       !local number of nodes
  integer                :: LNV  = 0       !local number of vertices
  integer                :: LNE  = 0       !local number of edges
  integer                :: LNF  = 0       !local number of faces
  logical                :: Beam = .false. !ef is beam
end type

!Constants
type(fe_db_type), public, dimension(232), save :: FE_DB

contains

!-----------------------------------------------------------------------
! init: give initial values to FE_DB
!-----------------------------------------------------------------------
subroutine FE_DB_init()

FE_DB( 11) = fe_db_type('Rod',                            1, 2,2,0,0,.true. )
FE_DB( 22) = fe_db_type('Tapered beam',                   1, 3,2,0,0,.true. )
FE_DB( 41) = fe_db_type('Plane Stress Linear Triangle',   2, 3,3,3,0,.false.)
FE_DB( 42) = fe_db_type('Plane Stress Parabolic Triangle',2, 6,3,3,0,.false.)
FE_DB( 94) = fe_db_type('Thin Shell Linear Quadrilateral',3, 4,4,4,0,.false.)
FE_DB(111) = fe_db_type('Solid Linear Tetrahedron',       3, 4,4,6,4,.false.)
FE_DB(118) = fe_db_type('Solid Parabolic Tetrahedron',    3,10,4,6,4,.false.)

end subroutine

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


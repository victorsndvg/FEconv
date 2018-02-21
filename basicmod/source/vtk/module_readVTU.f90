module module_readVTU_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es and Andres Prieto, andres.prieto(at)usc.es
!! Date: 31/12/2017
!!
!! Module for reading VTU binary files.  
!!
!! The steps to read a VTU file are:  
!! 1. [[module_readvtu_bmod(module):VTU_read_open(subroutine)]]: open a VTU file.  
!! 2. [[module_readvtu_bmod(module):VTU_read_mesh(subroutine)]]: read mesh.  
!! 3. [[module_readvtu_bmod(module):VTU_read_data(function)]]: read node/cell field when exists.  
!! 4. [[module_readvtu_bmod(module):VTU_read_close(subroutine)]]: close the VTU file.  
!!
!! Here we show a complete example:  
!!
!! `program test_readvtu`  
!! `use basicmod, only: real64, VTU_read_open, VTU_read_mesh, VTU_read_data, VTU_read_close`  
!! `implicit none`  
!! `integer                   :: nel      ! Global number of elements`  
!! `integer                   :: nnod     ! Global number of nodes`  
!! `integer,      allocatable :: nn(:,:)  ! Mesh connectivities`  
!! `real(real64), allocatable :: z(:,:)   ! Node coordinates`  
!! `character(20)             :: celldesc ! Cell description (see `[[module_readvtu_bmod(module):VTU_read_mesh(subroutine)]] 
!! `documentation)`  
!! `real(real64), allocatable :: nfield   ! Node field`  
!! `real(real64), allocatable :: cfield   ! Cell field`  
!! `character(20)             :: typef    ! Field type (see `[[module_readvtu_bmod(module):VTU_read_data(function)]] 
!! `documentation)` 
!! `!  `
!! `call VTU_read_open('test.vtu')`  
!! `call VTU_read_mesh(nel, nnod, nn, z, celldesc)`  
!! `if  (VTU_read_data(nfield, 'pd1', 'node', typef) /= 0) print*, 'Field pd1 not found.'`  
!! `if  (VTU_read_data(cfield, 'pd5', 'cell', typef) /= 0) print*, 'Field pd5 not found.'`  
!! `call VTU_read_close()`  
!! `end program`  
use module_compiler_dependant_bmod, only: real64
use module_report_bmod, only: info, error
use module_alloc_bmod, only: alloc, dealloc
use module_convers_bmod, only: string
use lib_vtk_io_bmod
use lib_vtk_io_read_bmod
implicit none

!Private procedures
private :: cell_type

contains

!-----------------------------------------------------------------------
! VTU_read_open
!-----------------------------------------------------------------------
subroutine VTU_read_open(ficho)
!! Open a VTU file for reading.
character(*), intent(in) :: ficho !! Input VTU file.
integer(I4P):: np

if (vtk_ini_xml_read('BINARY', ficho, 'UnstructuredGrid',np) /= 0) &
  call error('(module_readvtu/VTU_read_open) vtk_ini_xml_read error.')
call info('(module_readvtu/VTU_read_open) reading a VTU binary file '//trim(ficho)//'...')
end subroutine

!-----------------------------------------------------------------------
! VTU_read_mesh
!-----------------------------------------------------------------------
subroutine VTU_read_mesh(nel, nnod, nn, z, celldesc)
!! Read mesh data.  
!!
!! @note The valid cell descriptions are identified by a number that corresponds with those exposed in Figures 2. and 3. of 
!! [File Formats for vtk](www.vtk.org/VTK/img/file-formats.pdf). In particular, the cells accepted by this module are: 
!! `3_I1P` (line), `5_I1P` (triangle), 22_I1P` (triangle2), `9_I1P`(quad), `10_I1P`(tetra), `24_I1P`(tetra2), 
!! `12_I1P` (hexahedron), `13_I1P` (wedge), `14_I1P` (pyramid).  
integer(I4P), intent(out)              :: nel      !! Global number of elements.
integer(I4P), intent(out)              :: nnod     !! Global number of nodes.
integer(I4P), intent(out), allocatable :: nn(:,:)  !! Mesh conectivities.
real(R8P),    intent(out), allocatable :: z(:,:)   !! Node coordinates.
character(*), intent(out)              :: celldesc !! Cell description (see Note below).
real(R8P), allocatable :: z1(:), z2(:), z3(:)
integer(I1P), dimension(:), allocatable :: aux_CTYPE(:)
integer(I4P), dimension(:), allocatable :: aux_LNN(:), nn_aux(:)
! for the cell_type function
integer:: DIM   !! Space dimension.
integer:: LNV   !! Local number of vertices.
integer:: LNN   !! Local number of nodes.
integer:: LNE   !! Local number of edges.
integer:: LNF   !! Local number of faces.

! read vertices mesh
if (vtk_geo_xml_read(nnod, nel, z1, z2, z3) /=0 ) &
  call error('(module_readvtu/VTU_read_mesh) vtk_geo_xml_read error on coordinates.')
call alloc(z, 3, nnod)
z(1,:) = z1; z(2,:) = z2; z(3,:) = z3
! read connectivity arrays
if (vtk_con_xml_read(nel, nn_aux, aux_LNN, aux_CTYPE)/=0) &
  call error('(module_readvtu/VTU_read_mesh) vtk_geo_xml_read error on connectivities.')
! check the finite element type (assuming that all the elements have the same ctype)
if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, aux_CTYPE(1)) /= 0) &
  call error('(module_readvtu/VTU_read_mesh) cell_type error.')
DIM = size(z,1) ! change DIM for triangles in 3D
! To read arrays with rank bigger than 1 can crash with old ifort
call alloc(nn, LNN, nel)
nn = reshape(nn_aux+1, (/LNN,nel/))
call dealloc(nn_aux); call dealloc(aux_LNN); deallocate(aux_CTYPE);
call dealloc(z1);     call dealloc(z2);      call dealloc(z3)
end subroutine

!-----------------------------------------------------------------------
! VTU_read_data
!-----------------------------------------------------------------------
function VTU_read_data(field, namef, dataf, typef) result(res)
!! Check the existence of a node/cell field in the VTU file and read when it exists.  
!!
!! @note The valid field types are identified with the `NumberOfComponents` tag value (see specifications for XML file format 
!! in [File Formats for vtk](www.vtk.org/VTK/img/file-formats.pdf). In particular, the field types accepted by this module are: 
!! `1_I4P` (scalar), `3_I4P` (vector), 6_I4P` (symmetric_tensor).  
real(R8P), allocatable, intent(out) :: field(:) !! Field.
character(*),           intent(in)  :: namef    !! Field name.
character(*),           intent(in)  :: dataf    !! Data type (`node` or `cell`).
character(*),           intent(out) :: typef    !! Field type (see Note below).
integer                             :: res      !! zero when read was successful; non-zero otherwise.
integer(I4P) :: ncomp,nnod

call info('(module_readvtu/VTU_read_data) reading pointdata field '//trim(namef)//'...')
if (dataf /= 'node' .and. dataf /= 'cell') &
  call error('(module_readvtu/VTU_read_data) data type must be node or cell: '//trim(string(dataf)))
res = vtk_var_xml_read(dataf, nnod, ncomp, trim(namef), field)
if (res /= 0) return
select case(ncomp)
case (1_I4P)
  typef='scalar'
case (3_I4P)
  typef='vector'
case (6_I4P)
  typef='symmetric_tensor'
case default
  call error('(module_readvtu/VTU_read_data) NumberOfComponents value not implemented: '//trim(string(ncomp)))
end select
end function

!-----------------------------------------------------------------------
! VTU_read_close
!-----------------------------------------------------------------------
subroutine VTU_read_close()
!! Close the VTU file.

if (vtk_end_xml_read()/=0) call error('(module_readvtu/VTU_read_close) vtk_ini_xml_read error.')
call info('(module_readvtu/VTU_read_close) Done!')
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! cell_type
!-----------------------------------------------------------------------
function cell_type(name, DIM, LNV, LNN, LNE, LNF, CTYPE) result(res)
!! Give the associated dimensions of a FE.
character(*), intent(out) :: name  !! FE name.
integer,      intent(out) :: DIM   !! Space dimension.
integer,      intent(out) :: LNV   !! Local number of vertices.
integer,      intent(out) :: LNN   !1 Llocal number of nodes.
integer,      intent(out) :: LNE   !! local number of edges.
integer,      intent(out) :: LNF   !! Local number of faces.
integer(I1P), intent(in)  :: CTYPE !! Cell type.
integer                   :: res   !! 0 if the FE is known, -1 otherwise.

res = 0
select case(CTYPE)
 case(3_I1P)
  name='line'
  DIM = 1; LNV = 2; LNN = 2; LNE = 1;  LNF = 0
 case(5_I1P)
  name='triangle'
  DIM = 2; LNV = 3; LNN = 3; LNE = 3;  LNF = 1
 case(22_I1P)
  name='triangle2'
  DIM = 2; LNV = 3; LNN = 6; LNE = 3;  LNF = 1
 case(9_I1P)
  name='quad'
  DIM = 2; LNV = 4; LNN = 4; LNE = 4;  LNF = 1
 case(10_I1P)
  name='tetra'
  DIM = 3; LNV = 4; LNN = 4; LNE = 6;  LNF = 4
 case(24_I1P)
  name='tetra2'
  DIM = 3; LNV = 4; LNN = 10; LNE = 6;  LNF = 4
 case(12_I1P)
  name='hexahedron'
  DIM = 3; LNV = 8; LNN = 8; LNE = 12; LNF = 6
 case(13_I1P)
  name='wedge'
  DIM = 3; LNV = 6; LNN = 6; LNE = 9;  LNF = 5
 case(14_I1P)
  name='pyramid'
  DIM = 3; LNV = 5; LNN = 5; LNE = 8;  LNF = 5
case default
  res = -1
end select
end function

end module

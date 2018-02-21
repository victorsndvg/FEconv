module module_writeVTU_bmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module for writing VTU binary files.  
!!
!! The steps to write a vtu are ([] means optional):  
!! 1.  [[module_writeVTU_bmod(module):VTU_open(subroutine)]]: opens a VTU file.  
!! 2.  [[module_writeVTU_bmod(module):VTU_write_mesh(subroutine)]]: writes the mesh.  
!! [3.] [[module_writeVTU_bmod(module):VTU_begin_pointdata(subroutine)]]: begins the poindata record.  
!! [4.] [[module_writeVTU_bmod(module):VTU_write_pointdata(subroutine)]]: writes a pointdata field.  
!! [5.] [[module_writeVTU_bmod(module):VTU_end_pointdata(subroutine)]]: ends the poindata record.  
!! [6.] [[module_writeVTU_bmod(module):VTU_begin_celldata(subroutine)]]: begins the celldata record.  
!! [7.] [[module_writeVTU_bmod(module):VTU_write_celldata(subroutine)]]: writes a celldata field.  
!! [8.] [[module_writeVTU_bmod(module):VTU_end_celldata(subroutine)]]: ends celldata record.  
!! 9.  [[module_writeVTU_bmod(module):VTU_close(subroutine)]]: closes the file.  
!!
!! The previous steps can be summarized calling [[module_writeVTU_bmod(module):writeVTU(subroutine)]] if you only attach 
!! a single field in the mesh file.  
!!
!! Here we show a complete example using steps 1.-9. For a simpler example, check 
!! [[module_writeVTU_bmod(module):writeVTU(subroutine)]] documentation.  
!!
!! `program test_writeVTU`  
!! `use basicmod, only: real64, VTU_open, VTU_write_mesh, VTU_begin_pointdata, VTU_write_pointdata, VTU_end_pointdata, &  
!!                                                        VTU_begin_celldata,  VTU_write_celldata,  VTU_end_celldata, VTU_close`  
!! `implicit none`  
!! `integer :: i`  
!! `integer                   :: nel = 2  ! Number of cells`  
!! `integer                   :: nnod = 5 ! Number of points`  
!! `integer, allocatable      :: nn(:,:)  ! Tetrahedron connectivities, 4 x nel`  
!! `real(real64), allocatable :: z(:,:)   ! Point coordinates, 3 x nnod`  
!! `real(real64), allocatable :: field(:) ! field array`  
!!
!! `allocate(nn(4,2)) ! Tetrahedron connectivities, 4 x nel`  
!! `nn(:,1) = [1,2,3,4]; nn(:,2) = [1,3,2,5]`  
!! `allocate(z(3,5)) ! Point coordinates, 3 x nnod`  
!! `z(:,1) = [ 0.,  0.,  0.]`  
!! `z(:,2) = [ 1.,  0.,  0.]`  
!! `z(:,3) = [ 0.,  1.,  0.]`  
!! `z(:,4) = [0.3, 0.3,  1.]`  
!! `z(:,5) = [0.3, 0.3, -1.]`  
!! `call VTU_open('test.vtu')                              ! Open file`  
!! `call VTU_write_mesh(nel, nnod, nn, z, 'tetra')         ! Write mesh`  
!! `call VTU_begin_pointdata(Scalars='pd2', Vectors='pd4') ! Begin pointdata record`  
!! `allocate(field(5))`  
!! `field = [(real(1,8)*i, i = 1,5)]`  
!! `call VTU_write_pointdata(field, 'pd1', 'scalar')       ! Write pointdata scalar field`  
!! `field = [(real(1,8)*i, i = 5,1,-1)]`  
!! `call VTU_write_pointdata(field, 'pd2', 'scalar')       ! Write pointdata scalar field`  
!! `deallocate(field); allocate(field(15))`  
!! `field = [(real(1,8)*i, i = 1,15)]`  
!! `call VTU_write_pointdata(field, 'pd3', 'vector')       ! Write pointdata vector field`  
!! `field = [(real(1,8)*i, i = 15,1,-1)]`  
!! `call VTU_write_pointdata(field, 'pd4', 'vector')       ! Write pointdata vector field`  
!! `call VTU_end_pointdata()                               ! End pointdata record`  
!! `call VTU_begin_celldata(Scalars='pd5', Vectors='pd7')  ! Begin celldata record`  
!! `deallocate(field); allocate(field(2))`  
!! `field = [(real(1,8)*i, i = 1,2)]`  
!! `call VTU_write_celldata(field, 'pd5', 'scalar')        ! Write celldata scalar field`  
!! `field = [(real(1,8)*i, i = 2,1,-1)]`  
!! `call VTU_write_celldata(field, 'pd6', 'scalar')        ! Write celldata scalar field`  
!! `deallocate(field); allocate(field(6))`  
!! `field = [(real(1,8)*i, i = 1,6)]`  
!! `call VTU_write_celldata(field, 'pd7', 'vector')        ! Write celldata vector field`  
!! `field = [(real(1,8)*i, i = -1,-6,-1)]`  
!! `call VTU_write_celldata(field, 'pd8', 'vector')        ! Write celldata vector field`  
!! `call VTU_end_celldata()                                ! End celldata record`  
!! `call VTU_close()                                       ! Close file`  
!! `end program`  
use module_report_bmod, only: info, error
use module_convers_bmod, only: lcase
use lib_vtk_io_bmod
implicit none

!Private variables
integer,      private :: DIM   !space dimension
integer,      private :: LNV   !local number of vertices
integer,      private :: LNN   !local number of nodes
integer,      private :: LNE   !local number of edges
integer,      private :: LNF   !local number of faces
integer(I1P), private :: CTYPE !cell type
integer(I4P), private :: nl    !global number of cells
integer(I4P), private :: nd    !global number of nodes

!Private procedures
private :: cell_type

contains

!-----------------------------------------------------------------------
! VTU_open: open a vtu file
!-----------------------------------------------------------------------
subroutine VTU_open(ficho)
!! Open a VTU file.
character(*), intent(in) :: ficho !! Output VTU file.
if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) call error('(module_writeVTU/VTU_open) vtk_ini_xml error.')
call info('(module_writeVTU/VTU_open) writing a VTU binary file '//trim(ficho)//'...')
end subroutine

!-----------------------------------------------------------------------
! VTU_write_mesh: writes mesh info
!-----------------------------------------------------------------------
subroutine VTU_write_mesh(nel, nnod, nn, z, celldesc)
!! Write the mesh.
!!
!! @note The valid cell descriptions are: `line`, `triangle`, `triangle2`, `quad`, `tetra`, `tetra2`, 
!! `hexahedron`, `wedge` and `pyramid`.  
integer(I4P), intent(in) :: nel      !! Global number of elements.
integer(I4P), intent(in) :: nnod     !! Global number of nodes.
integer(I4P)             :: nn(:,:)  !! Nodes index array.
real(R8P)                :: z(:,:)   !! Vertices coordinates array.
character(*), intent(in) :: celldesc !! Cell description; see the Note below.
real(R8P), allocatable :: z1(:), z2(:), z3(:)
integer(I4P), allocatable :: aux_LNN(:)
integer(I1P), allocatable :: aux_CTYPE(:)
integer(I4P), allocatable :: nn_aux(:)
integer :: i

nl = nel; nd = nnod
if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) call error('(module_writeVTU/VTU_write_mesh) cell_type error.')
DIM = size(z,1) !change DIM for triangles in 3D
call info('(module_writeVTU/VTU_write_mesh) saving mesh data...')
!saving coords in the current piece
allocate(z1(nnod),z2(nnod),z3(nnod))
if (DIM < 3) then
  z1=z(1,:);z2=z(2,:);z3=(/(real(0,R8P),i=1,nnod)/)
  if (vtk_geo_xml(nnod, nel, z1, z2, z3) /= 0) call error('(module_writeVTU/VTU_write_mesh) vtk_geo_xml error 1.')
else
  z1=z(1,:);z2=z(2,:);z3=z(3,:)
  if (vtk_geo_xml(nnod, nel, z1, z2, z3) /= 0) call error('(module_writeVTU/VTU_write_mesh) vtk_ini_xml error 2.')
end if
!saving connectivity
! To insert reshape as argument can crash with old ifort (without heap arrays 0)
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) &
!   call error('(module_writeVTU/VTU_write_mesh) vtk_con_xml error.')
allocate(nn_aux(LNN*nel),aux_LNN(nel),aux_CTYPE(nel))
nn_aux = reshape(nn-1, (/LNN*nel/))
aux_LNN=(/(i*LNN,i=1,nel)/);aux_CTYPE=(/(CTYPE,i=1,nel)/)
if (vtk_con_xml(nel, nn_aux, aux_LNN, aux_CTYPE) /= 0) &
  call error('(module_writeVTU/VTU_write_mesh) vtk_con_xml error.')
deallocate(nn_aux,aux_LNN,aux_CTYPE,z1,z2,z3)
end subroutine

!-----------------------------------------------------------------------
! VTU_begin_pointdata: permits the writing of poindata fields
!-----------------------------------------------------------------------
subroutine VTU_begin_pointdata(Scalars, Vectors, Normals, Tensors, TCoords)
!! Begin the poindata record.
character(*), intent(in), optional :: Scalars !! Name of the active scalars array, if any.
character(*), intent(in), optional :: Vectors !! Name of the active vectors array, if any.
character(*), intent(in), optional :: Normals !! Name of the active normals array, if any.
character(*), intent(in), optional :: Tensors !! Name of the active tensors array, if any.
character(*), intent(in), optional :: TCoords !! Name of the active texture coordinates array, if any.

if (vtk_dat_xml('node', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) &
  call error('(module_writeVTU/VTU_begin_pointdata) vtk_dat_xml error.')
end subroutine

!-----------------------------------------------------------------------
! VTU_write_pointdata: writes a pointdata field
!-----------------------------------------------------------------------
subroutine VTU_write_pointdata(field, namef, typef)
!! Write a pointdata real64 field.  
real(R8P), intent(in)    :: field(:) !! Pointdata real64 field array.
character(*), intent(in) :: namef    !! Field name (e.g. `Temperature`).
character(*), intent(in) :: typef    !! Field type, `scalar` or `vector`.
real(R8P), allocatable :: z1(:), z2(:), z3(:)
integer :: i

call info('(module_writeVTU/VTU_write_pointdata) saving pointdata field '//trim(namef)//'...')
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nd, trim(namef), field) /= 0) call error('(module_writeVTU/VTU_write_pointdata) vtk_var_xml error.')
case('vector')
  allocate(z1(nd),z2(nd),z3(nd))
  select case(DIM)
  case (2)
    z1=field(1:nd*DIM:DIM);z2=field(2:nd*DIM:DIM);z3=(/(real(0,R8P),i=1,nd)/)
    if (vtk_var_xml(nd, trim(namef), z1, z2, z3) /= 0) call error('(module_writeVTU/VTU_write_pointdata) vtk_var_xml error.')
  case (3)
    z1=field(1:nd*DIM:DIM);z2=field(2:nd*DIM:DIM);z3=field(3:nd*DIM:DIM)
    if (vtk_var_xml(nd, trim(namef), z1, z2, z3) /= 0) call error('(module_writeVTU/VTU_write_pointdata) vtk_var_xml error.')
  case default; call error('(module_writeVTU/VTU_write_pointdata) Not implemented.')
  end select
  deallocate(z1,z2,z3)
case default; call error('(module_writeVTU/VTU_write_pointdata) Not implemented.')
end select
end subroutine

!-----------------------------------------------------------------------
! VTU_write_pointdata_R4P: writes a pointdata field
!-----------------------------------------------------------------------
subroutine VTU_write_pointdata_R4P(field, namef, typef)
!! Write a pointdata real field.
real(R4P)                :: field(:) !! Pointdata real64 field array.
character(*), intent(in) :: namef    !! Field name (e.g. `Temperature`).
character(*), intent(in) :: typef    !! Field type, `scalar` or `vector`.
integer :: i

call info('(module_writeVTU/VTU_write_pointdata_R4P) saving pointdata field '//trim(namef)//'...')
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nd, trim(namef), field) /= 0) call error('(module_writeVTU/VTU_write_pointdata_R4P) vtk_var_xml error.')
case('vector')
  select case(DIM)
  case (2)
    if (vtk_var_xml(nd, trim(namef), field(1:nd*DIM:DIM), field(2:nd*DIM:DIM), &
    (/(real(0,R4P),i=1,nd)/)) /= 0) call error('(module_writeVTU/VTU_write_pointdata_R4P) vtk_var_xml error.')
  case (3)
    if (vtk_var_xml(nd, trim(namef), field(1:nd*DIM:DIM), field(2:nd*DIM:DIM), &
    field(3:nd*DIM:DIM)) /= 0) call error('(module_writeVTU/VTU_write_pointdata_R4P) vtk_var_xml error.')
  case default; call error('(module_writeVTU/VTU_write_pointdata_R4P) Not implemented.')
  end select
case default; call error('(module_writeVTU/VTU_write_pointdata_R4P) Not implemented.')
end select
end subroutine

!-----------------------------------------------------------------------
! VTU_end_pointdata: stops the writing of poindata fields
!-----------------------------------------------------------------------
subroutine VTU_end_pointdata()
!! End the poindata record.

if (vtk_dat_xml('node', 'CLOSE')/=0) call error('(module_writeVTU/VTU_end_pointdata) vtk_dat_xml error.')
end subroutine

!-----------------------------------------------------------------------
! VTU_begin_celldata: permits the writing of celldata fields
!-----------------------------------------------------------------------
subroutine VTU_begin_celldata(Scalars, Vectors, Normals, Tensors, TCoords)
!! Begin the celldata record.
character(*), intent(in), optional :: Scalars !! Name of the active scalars array, if any.
character(*), intent(in), optional :: Vectors !! Name of the active vectors array, if any.
character(*), intent(in), optional :: Normals !! Name of the active normals array, if any.
character(*), intent(in), optional :: Tensors !! Name of the active tensors array, if any.
character(*), intent(in), optional :: TCoords !! Name of the active texture coordinates array, if any.

if (vtk_dat_xml('cell', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) &
  call error('(module_writeVTU/VTU_begin_celldata) vtk_dat_xml error.')
end subroutine

!-----------------------------------------------------------------------
! VTU_write_celldata: writes a celldata field
!-----------------------------------------------------------------------
subroutine VTU_write_celldata(field, namef, typef)
!! Write a celldata real64 field.
real(R8P)                :: field(:) !! Celldata real64 field array.
character(*), intent(in) :: namef    !! Field name (e.g. `Heat flux`).
character(*), intent(in) :: typef    !! Field type, `scalar` or `vector`.
integer :: i

call info('(module_writeVTU/VTU_write_celldata) saving celldata field '//trim(namef)//'...')
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nl, trim(namef), field) /= 0) call error('(module_writeVTU/VTU_write_celldata) vtk_var_xml error.')
case('vector')
  select case(DIM)
  case (2)
    if (vtk_var_xml(nl, trim(namef), field(1:nl*DIM:DIM), field(2:nl*DIM:DIM), &
    (/(real(0,R8P),i=1,nl)/)) /= 0) call error('(module_writeVTU/VTU_write_celldata) vtk_var_xml error.')
  case (3)
    if (vtk_var_xml(nl, trim(namef), field(1:nl*DIM:DIM), field(2:nl*DIM:DIM), &
    field(3:nl*DIM:DIM)) /= 0) call error('(module_writeVTU/VTU_write_celldata) vtk_var_xml error.')
  case default; call error('(module_writeVTU/VTU_write_celldata) Not implemented.')
  end select
case default; call error('(module_writeVTU/VTU_write_celldata) Not implemented.')
end select
end subroutine

!-----------------------------------------------------------------------
! VTU_end_celldata: stops the writing of cellata fields
!-----------------------------------------------------------------------
subroutine VTU_end_celldata()
!! End celldata record.

if (vtk_dat_xml('cell', 'CLOSE')/=0) call error('(module_writeVTU/VTU_end_celldata) vtk_dat_xml error.')
end subroutine

!-----------------------------------------------------------------------
! VTU_close: closes a vtu file
!-----------------------------------------------------------------------
subroutine VTU_close(geo)
!! Closes the file.
logical, optional :: geo

if(present(geo)) then
  if(geo) then
    if (vtk_geo_xml()/=0) call error('(module_writeVTU/VTU_close) vtk_geo_xml error.')
  endif
else
  if (vtk_geo_xml()/=0) call error('(module_writeVTU/VTU_close) vtk_geo_xml error.')
endif
call info('(module_writeVTU/VTU_close) saving appended data...')
if (vtk_end_xml()/=0) call error('(module_writeVTU/VTU_close) vtk_geo_xml error.')
end subroutine

!-----------------------------------------------------------------------
! Write a RECONVXX mesh (and a B field) into a VTU binary file
! Writes reals in double precision (8 bytes) - R8P
!-----------------------------------------------------------------------
subroutine writeVTU(nel, nnod, nn, z, celldesc, field, namef, typef, DOFType, ficho, &
Scalars, Vectors, Normals, Tensors, TCoords)
!! Save a mesh and a single real64 field into a VTU file.  
!! __Example:__  
!! `program test_writeVTU`  
!! `use basicmod, only: real64, writeVTU`  
!! `implicit none`  
!! `integer      :: nel = 2  ! Number of cells`  
!! `integer      :: nnod = 5 ! Number of points`  
!! `integer      :: nn(4,2)  ! Tetrahedron connectivities, 4 x nel`  
!! `real(real64) :: z(3,5)   ! Point coordinates, 3 x nnod`  
!! `real(real64) :: field(5) ! Pointdata field array`  
!! `integer :: i`  
!! `nn(:,1) = [1,2,3,4]; nn(:,2) = [1,3,2,5] ! Tetrahedron connectivities`  
!! `z(:,1) = [ 0.,  0.,  0.] ! Point coordinates`  
!! `z(:,2) = [ 1.,  0.,  0.]`  
!! `z(:,3) = [ 0.,  1.,  0.]`  
!! `z(:,4) = [0.3, 0.3,  1.]`  
!! `z(:,5) = [0.3, 0.3, -1.]`  
!! `call writeVTU(nel, nnod, nn, z, 'tetra', field, 'Potential (V)', 'scalar', 'node', 'poten.vtu')`  
!! `end program`  
!!
!! @note The valid cell descriptions are: `line`, `triangle`, `triangle2`, `quad`, `tetra`, `tetra2`, 
!! `hexahedron`, `wedge` and `pyramid`.  
integer(I4P), intent(in) :: nel      !! Global number of elements.
integer(I4P), intent(in) :: nnod     !! Global number of nodes.
integer(I4P)             :: nn(:,:)  !! Nodes index array.
real(R8P)                :: z(:,:)   !! Nodes coordinates array.
character(*), intent(in) :: celldesc !! Cell description; see the Note below.
real(R8P)                :: field(:) !! Real64 field array.
character(*), intent(in) :: namef    !! Field name (e.g. `Temperature`).
character(*), intent(in) :: typef    !! Field type, `scalar` or `vector`.
character(*), intent(in) :: DOFtype  !! Field type, `node` or `cell`.
character(*), intent(in) :: ficho    !! Output VTU file.
character(*), intent(in), optional :: Scalars !! Name of the active scalars array, if any.
character(*), intent(in), optional :: Vectors !! Name of the active vectors array, if any.
character(*), intent(in), optional :: Normals !! Name of the active normals array, if any.
character(*), intent(in), optional :: Tensors !! Name of the active tensors array, if any.
character(*), intent(in), optional :: TCoords !! Name of the active texture coordinates array, if any.
integer      :: DIM   !! Space dimension.
integer      :: LNV   !! Local number of vertices.
integer      :: LNE   !! Local number of edges.
integer      :: LNF   !! Local number of faces.
integer      :: LNN   !! Local number of nodes.
integer(I1P) :: CTYPE !! Cell linear type, according to VTK format.
integer(I4P), dimension(:), allocatable  :: nn_aux
integer :: i

if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) &
  call error('(module_writeVTU/writeVTU) cell_type error.')
call info('(module_writeVTU/writeVTU) writing a VTU binary file...')
call info('(module_writeVTU/writeVTU) saving mesh data...')
! Create output
if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) call error('(module_writeVTU/writeVTU) vtk_ini_xml error.')
! Saving coords in the current piece
if (DIM < 3) then
  if (vtk_geo_xml(nnod, nel, z(1,:), z(2,:), (/(real(0,R8P),i=1,nnod)/)) /= 0) &
    call error('(module_writeVTU/writeVTU) vtk_geo_xml error 1.')
else
  if (vtk_geo_xml(nnod, nel, z(1,:), z(2,:), z(3,:)) /= 0) &
    call error('(module_writeVTU/writeVTU) vtk_geo_xml error 2.')
end if
! Saving connectivity
! To insert reshape as argument can crash with old ifort (without heap arrays 0)
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) &
!   call error('(module_writeVTU/writeVTU) vtk_con_xml error.')
allocate(nn_aux(LNN*nel))
nn_aux = reshape(nn-1, (/LNN*nel/))
if (vtk_con_xml(nel, nn_aux, (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) &
  call error('(module_writeVTU/writeVTU) vtk_con_xml error.')
deallocate(nn_aux) ! fix ifort
call info('(module_writeVTU/writeVTU) saving field data...')
select case(trim(lcase(DOFType)))
case('node')
  ! Node fields
  if (vtk_dat_xml('node', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) &
    call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nnod, trim(namef), field) /= 0) call error('(module_writeVTU/writeVTU) vtk_var_xml error.')
  case('vector')
    select case(DIM)
    case (2)
      if (vtk_var_xml(nnod, trim(namef), field(1:nnod*DIM:DIM), field(2:nnod*DIM:DIM), &
      (/(real(0,R8P),i=1,nnod)/)) /= 0) call error('(module_writeVTU/writeVTU) vtk_var_xml error.')
    case (3)
      if (vtk_var_xml(nnod, trim(namef), field(1:nnod*DIM:DIM), field(2:nnod*DIM:DIM), &
      field(3:nnod*DIM:DIM)) /= 0) call error('(module_writeVTU/writeVTU) vtk_var_xml error.')
    case default; call error('(module_writeVTU/writeVTU) not implemented (node).')
    end select
  case default; call error('(module_writeVTU/writeVTU) not implemented (node).')
  end select
  if (vtk_dat_xml('node', 'CLOSE')/=0) call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
case('cell')
  ! Cell Fields
  if (vtk_dat_xml('cell', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) &
    call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nel, trim(namef), field) /= 0) call error('(module_writeVTU/writeVTU) vtk_var_xml error.')
  case('vector')
    select case(DIM)
    case (2)
      if (vtk_var_xml(nel, trim(namef), field(1:nel*DIM:DIM), field(2:nel*DIM:DIM), &
      (/(real(0,R8P),i=1,nel)/)) /= 0) call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
    case (3)
      if (vtk_var_xml(nel, trim(namef), field(1:nel*DIM:DIM), field(2:nel*DIM:DIM), &
      field(3:nel*DIM:DIM)) /= 0) call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
    case default; call error('(module_writeVTU/writeVTU) not implemented (cell).')
    end select
  case default;  call error('(module_writeVTU/writeVTU) not implemented (cell).')
  end select
  if (vtk_dat_xml('cell', 'CLOSE')/=0) call error('(module_writeVTU/writeVTU) vtk_dat_xml error.')
case default; call error('(module_writeVTU/writeVTU) not implemented (close).')
end select
if (vtk_geo_xml()/=0) call error('(module_writeVTU/writeVTU) vtk_geo_xml error.')
call info('(module_writeVTU/writeVTU) saving appended data...')
if (vtk_end_xml()/=0) call error('(module_writeVTU/writeVTU) vtk_end_xml error.')
call info('(module_writeVTU/writeVTU) done!')
end subroutine

!-----------------------------------------------------------------------
! Write a RECONVXX mesh (and a B field) into a VTU binary file
! Modified to write reals in single precision (4 bytes) - R4P
!-----------------------------------------------------------------------
subroutine writeVTU4(nel, nnod, nn, z, celldesc, field, namef, typef, DOFType, ficho)
!! Save a mesh and a single real64 field into a VTU file with R4P precision (4 bytes).  
!! __Example:__  
!! `program test_writeVTU`  
!! `use basicmod, only: real64, writeVTU4`  
!! `implicit none`  
!! `integer      :: nel = 2  ! Number of cells`  
!! `integer      :: nnod = 5 ! Number of points`  
!! `integer      :: nn(4,2)  ! Tetrahedron connectivities, 4 x nel`  
!! `real(real64) :: z(3,5)   ! Point coordinates, 3 x nnod`  
!! `real(real64) :: field(5) ! Pointdata field array`  
!! `integer :: i`  
!! `nn(:,1) = [1,2,3,4]; nn(:,2) = [1,3,2,5] ! Tetrahedron connectivities`  
!! `z(:,1) = [ 0.,  0.,  0.] ! Point coordinates`  
!! `z(:,2) = [ 1.,  0.,  0.]`  
!! `z(:,3) = [ 0.,  1.,  0.]`  
!! `z(:,4) = [0.3, 0.3,  1.]`  
!! `z(:,5) = [0.3, 0.3, -1.]`  
!! `call writeVTU4(nel, nnod, nn, z, 'tetra', field, 'Potential (V)', 'scalar', 'node', 'poten.vtu')`  
!! `end program`  
!!
!! @note The valid cell descriptions are: `line`, `triangle`, `triangle2`, `quad`, `tetra`, `tetra2`,  
!! `hexahedron`, `wedge` and `pyramid`.  
integer(I4P), intent(in) :: nel      !! Global number of elements.
integer(I4P), intent(in) :: nnod     !! Global number of nodes.
integer(I4P)             :: nn(:,:)  !! Nodes index array.
real(R8P)                :: z(:,:)   !! Nodes coordinates array.
character(*), intent(in) :: celldesc !! Cell description; see the Note below.
real(R8P)                :: field(:) !! Real64 field array.
character(*), intent(in) :: namef    !! Field name (e.g. `Temperature`).
character(*), intent(in) :: typef    !! Field type, `scalar` or `vector`.
character(*), intent(in) :: DOFtype  !! Field type, `node` or `cell`.
character(*), intent(in) :: ficho    !! Output VTU file.
integer      :: DIM   !! Space dimension.
integer      :: LNV   !! Local number of vertices.
integer      :: LNE   !! Local number of edges.
integer      :: LNF   !! Local number of faces.
integer      :: LNN   !! Local number of nodes.
integer(I1P) :: CTYPE !! Cell linear type, according to VTK format.
integer(I4P), allocatable  :: nn_aux(:)
integer :: i

if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) &
  call error('(module_writeVTU/writeVTU4) cell_type error.')
call info('(module_writeVTU/writeVTU4) writing a VTU binary file...')
call info('(module_writeVTU/writeVTU4) saving mesh data...')
! Create output
if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) call error('(module_writeVTU/writeVTU4) vtk_ini_xml error.')
! Saving coords in the current piece
if (DIM < 3) then
  if (vtk_geo_xml(nnod, nel, real(z(1,:),R4P), real(z(2,:),R4P), (/(real(0,R4P),i=1,nnod)/)) /= 0) &
    call error('(module_writeVTU/writeVTU4) vtk_geo_xml error 1.')
else
  if (vtk_geo_xml(nnod, nel, real(z(1,:),R4P), real(z(2,:),R4P), real(z(3,:),R4P)) /= 0) &
    call error('(module_writeVTU/writeVTU4) vtk_geo_xml error 2.')
end if
! Saving connectivity
! To insert reshape as argument can crash with old ifort (without heap arrays 0)
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0)  &
!   call error('(module_writeVTU/writeVTU4) vtk_con_xml error.')
allocate(nn_aux(LNN*nel))
nn_aux = reshape(nn-1, (/LNN*nel/))
if (vtk_con_xml(nel, nn_aux, (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) &
  call error('(module_writeVTU/writeVTU4) vtk_con_xml error.')
deallocate(nn_aux)
call info('(module_writeVTU/writeVTU4) saving field data...')
select case(trim(lcase(DOFType)))
case('node')
  ! Node Fields
  if (vtk_dat_xml('node', 'OPEN')/=0) call error('(module_writeVTU/writeVTU4) vtk_dat_xml error.')
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nnod, trim(namef), real(field,R4P)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
  case('vector')
    select case(DIM)
    case (2)
      if (vtk_var_xml(nnod, trim(namef),&
        &real(field(1:nnod*DIM:DIM),R4P),&
        &real(field(2:nnod*DIM:DIM),R4P),&
        &(/(real(0,R4P),i=1,nnod)/)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
    case (3)
      if (vtk_var_xml(nnod, trim(namef),&
        &real(field(1:nnod*DIM:DIM),R4P),&
        &real(field(2:nnod*DIM:DIM),R4P),&
        &real(field(3:nnod*DIM:DIM),R4P)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
    case default; call error('(module_writeVTU/writeVTU4) not implemented.')
    end select
  case default; call error('(module_writeVTU/writeVTU) not implemented.')
  end select
  if (vtk_dat_xml('node', 'CLOSE')/=0) call error('(module_writeVTU/writeVTU4) vtk_dat_xml error.')
case('cell')
  ! Cell Fields
  if (vtk_dat_xml('cell', 'OPEN')/=0) call error('(module_writeVTU/writeVTU4) vtk_dat_xml error.')
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nel, trim(namef), real(field,R4P)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
  case('vector')
    select case(DIM)
    case (2)
      if (vtk_var_xml(nel, trim(namef),&
        &real(field(1:nel*DIM:DIM),R4P),&
        &real(field(2:nel*DIM:DIM),R4P),&
        &(/(real(0,R4P),i=1,nel)/)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
    case (3)
      if (vtk_var_xml(nel, trim(namef),&
        &real(field(1:nel*DIM:DIM),R4P),&
        &real(field(2:nel*DIM:DIM),R4P),&
        &real(field(3:nel*DIM:DIM),R4P)) /= 0) call error('(module_writeVTU/writeVTU4) vtk_var_xml error.')
    case default; call error('(module_writeVTU/writeVTU4) not implemented.')
    end select
  case default; call error('(module_writeVTU/writeVTU4) not implemented.')
  end select
  if (vtk_dat_xml('cell', 'CLOSE')/=0) call error('(module_writeVTU/writeVTU4) vtk_dat_xml error.')
case default; call error('(module_writeVTU/writeVTU4) not implemented.')
end select
if (vtk_geo_xml()/=0) call error('(module_writeVTU/writeVTU4) vtk_geo_xml error.')
call info('(module_writeVTU/writeVTU4) saving appended data...')
if (vtk_end_xml()/=0) call error('(module_writeVTU/writeVTU4) vtk_end_xml error.')
call info('(module_writeVTU/writeVTU4) done!')
end subroutine

!-----------------------------------------------------------------------
! PRIVATE PROCEDURES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! cell_type: give the associated dimensions of a FE
! RESULT: 0 if the FE is known, -1 otherwise
!-----------------------------------------------------------------------
function cell_type(name, DIM, LNV, LNN, LNE, LNF, CTYPE) result(res)

character(*), intent(in)  :: name  !FE name
integer,          intent(out) :: DIM   !space dimension
integer,          intent(out) :: LNV   !local number of vertices
integer,          intent(out) :: LNN   !local number of nodes
integer,          intent(out) :: LNE   !local number of edges
integer,          intent(out) :: LNF   !local number of faces
integer(I1P),     intent(out) :: CTYPE !cell type
integer :: res

res = 0
select case(trim(adjustl(lcase(name))))
case('line')
  DIM = 1; LNV = 2; LNN = 2; LNE = 1;  LNF = 0; CTYPE = 3
case('triangle')
  DIM = 2; LNV = 3; LNN = 3; LNE = 3;  LNF = 1; CTYPE = 5
case('triangle2')
  DIM = 2; LNV = 3; LNN = 6; LNE = 3;  LNF = 1; CTYPE = 22
case('quad')
  DIM = 2; LNV = 4; LNN = 4; LNE = 4;  LNF = 1; CTYPE = 9
case('tetra')
  DIM = 3; LNV = 4; LNN = 4; LNE = 6;  LNF = 4; CTYPE = 10
case('tetra2')
  DIM = 3; LNV = 4; LNN = 10; LNE = 6;  LNF = 4; CTYPE = 24
case('hexahedron')
  DIM = 3; LNV = 8; LNN = 8; LNE = 12; LNF = 6; CTYPE = 12
case('wedge')
  DIM = 3; LNV = 6; LNN = 6; LNE = 9;  LNF = 5; CTYPE = 13
case('pyramid')
  DIM = 3; LNV = 5; LNN = 5; LNE = 8;  LNF = 5; CTYPE = 14
case default
  res = -1
end select
end function

end module

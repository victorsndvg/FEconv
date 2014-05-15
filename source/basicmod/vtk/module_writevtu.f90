module module_writevtu
!-----------------------------------------------------------------------
! module for writing VTU binary files
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran.pena(at)usc.es
!         LiB_VTK_IO is a creation of stefano.zaghi(at)gmail.com
! Last update: 27/01/2012
!
! PUBLIC PROCEDURES: the steps to write a vtu are ([] mean optional steps):
!     1  VTU_open: opens a vtu file 
!     2  VTU_write_mesh: writes mesh info 
!    [3] VTU_begin_pointdata: permits the writing of poindata fields
!    [4] VTU_write_pointdata: writes a pointdata field  
!    [5] VTU_end_pointdata: stops the writing of poindata fields
!    [6] VTU_begin_celldata: permits the writing of celldata fields
!    [7] VTU_write_celldata: writes a celldata field  
!    [8] VTU_end_celldata: stops the writing of cellata fields
!     9  VTU_close: closes a vtu file
!
! To write a single field you can use 'writeVTU' instead
!-----------------------------------------------------------------------
use lib_vtk_io
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
private :: cell_type, lcase

contains

!-----------------------------------------------------------------------
! VTU_open: open a vtu file 
!-----------------------------------------------------------------------
subroutine VTU_open(ficho)
character(len=*), intent(in) :: ficho              !output file

if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) stop 'writeVTU: vtk_ini_xml error'
print'(a)', 'Writing a VTU binary file '//trim(ficho)//'...'
end subroutine

!-----------------------------------------------------------------------
! VTU_write_mesh: writes mesh info 
!-----------------------------------------------------------------------
subroutine VTU_write_mesh(nel, nnod, nn, z, celldesc)

integer(I4P), intent(in) :: nel                    !global number of elements
integer(I4P), intent(in) :: nnod                   !global number of nodes
integer(I4P), dimension(:,:) :: nn    !nodes index array
real(R8P),    dimension(:,:) :: z     !vertices coordinates array
real(R8P), dimension(:), allocatable :: z1,z2,z3  
integer(I4P), dimension(:), allocatable :: aux_LNN
integer(I1P), dimension(:), allocatable :: aux_CTYPE
character(len=*), intent(in) :: celldesc           !cell description
integer :: i

integer(I4P), dimension(:), allocatable  :: nn_aux ! fix ifort

nl = nel; nd = nnod
!print*, '----------', 1
if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) stop 'writeVTU: cell_type error'
DIM = size(z,1) !change DIM for triangles in 3D
print'(a)', 'Saving mesh data...'

!saving coords in the current piece
allocate(z1(nnod),z2(nnod),z3(nnod))
if (DIM < 3) then
  z1=z(1,:);z2=z(2,:);z3=(/(real(0,R8P),i=1,nnod)/)
  if (vtk_geo_xml(nnod, nel, z1, z2, z3) /= 0) stop 'writeVTU: vtk_geo_xml error 1'
else
  z1=z(1,:);z2=z(2,:);z3=z(3,:)
  if (vtk_geo_xml(nnod, nel, z1, z2, z3) /= 0) stop 'writeVTU: vtk_ini_xml error 2'
end if

!saving connectivity
  allocate(nn_aux(LNN*nel),aux_LNN(nel),aux_CTYPE(nel)) ! fix ifort
  nn_aux = reshape(nn-1, (/LNN*nel/)) ! fix ifort
! esto poderia cascar con ifort e mallachip, en reshape: [sen heap arrays 0]
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) & 
!    stop 'writeVTU: vtk_con_xml error'
  aux_LNN=(/(i*LNN,i=1,nel)/);aux_CTYPE=(/(CTYPE,i=1,nel)/)
  if (vtk_con_xml(nel, nn_aux, aux_LNN, aux_CTYPE) /= 0) & 
     stop 'writeVTU: vtk_con_xml error' ! fix ifort
  deallocate(nn_aux,aux_LNN,aux_CTYPE,z1,z2,z3) ! fix ifort

end subroutine

!-----------------------------------------------------------------------
! VTU_begin_pointdata: permits the writing of poindata fields
!-----------------------------------------------------------------------
subroutine VTU_begin_pointdata(Scalars, Vectors, Normals, Tensors, TCoords)
character(len=*), intent(in), optional :: Scalars ! The name of the active scalars array, if any.
character(len=*), intent(in), optional :: Vectors ! The name of the active vectors array, if any.
character(len=*), intent(in), optional :: Normals ! The name of the active normals array, if any.
character(len=*), intent(in), optional :: Tensors ! The name of the active tensors array, if any.
character(len=*), intent(in), optional :: TCoords ! The name of the active texture coordinates array, if any.

if (vtk_dat_xml('node', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) stop

end subroutine

!-----------------------------------------------------------------------
! VTU_write_pointdata: writes a pointdata field  
!-----------------------------------------------------------------------
subroutine VTU_write_pointdata(field, namef, typef)

real(R8P), intent(in), dimension(:)   :: field !field
character(len=*), intent(in) :: namef              !field name (ex., temperature)
character(len=*), intent(in) :: typef              !field type (scalar/vector)
real(R8P), allocatable, dimension(:)   :: z1,z2,z3
integer :: i

print'(a)', 'Saving pointdata field '//trim(namef)//'...'
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nd, trim(namef), field) /= 0) stop
case('vector')
  allocate(z1(nd),z2(nd),z3(nd))
  select case(DIM)
  case (2)      
    z1=field(1:nd*DIM:DIM);z2=field(2:nd*DIM:DIM);z3=(/(real(0,R8P),i=1,nd)/)
    if (vtk_var_xml(nd, trim(namef), z1, z2, z3) /= 0) stop
  case (3)     
    z1=field(1:nd*DIM:DIM);z2=field(2:nd*DIM:DIM);z3=field(3:nd*DIM:DIM)
    if (vtk_var_xml(nd, trim(namef), z1, z2, z3) /= 0) stop
  case default; stop 'Not implemented'
  end select
  deallocate(z1,z2,z3)
case default; stop 'Not implemented'
end select

end subroutine

!-----------------------------------------------------------------------
! VTU_write_pointdata_R4P: writes a pointdata field  
!-----------------------------------------------------------------------
subroutine VTU_write_pointdata_R4P(field, namef, typef)

real(R4P),    dimension(:)   :: field !field
character(len=*), intent(in) :: namef              !field name (ex., temperature)
character(len=*), intent(in) :: typef              !field type (scalar/vector)
integer :: i

print'(a)', 'Saving pointdata field '//trim(namef)//'...'
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nd, trim(namef), field) /= 0) stop
case('vector')
  select case(DIM)
  case (2)      
    if (vtk_var_xml(nd, trim(namef), field(1:nd*DIM:DIM), field(2:nd*DIM:DIM), &
    (/(real(0,R4P),i=1,nd)/)) /= 0) stop
  case (3)      
    if (vtk_var_xml(nd, trim(namef), field(1:nd*DIM:DIM), field(2:nd*DIM:DIM), &
    field(3:nd*DIM:DIM)) /= 0) stop
  case default; stop 'Not implemented'
  end select
case default; stop 'Not implemented'
end select

end subroutine

!-----------------------------------------------------------------------
! VTU_end_pointdata: stops the writing of poindata fields
!-----------------------------------------------------------------------
subroutine VTU_end_pointdata()

if (vtk_dat_xml('node', 'CLOSE')/=0) stop

end subroutine

!-----------------------------------------------------------------------
! VTU_begin_celldata: permits the writing of celldata fields
!-----------------------------------------------------------------------
subroutine VTU_begin_celldata(Scalars, Vectors, Normals, Tensors, TCoords)
character(len=*), intent(in), optional :: Scalars ! The name of the active scalars array, if any.
character(len=*), intent(in), optional :: Vectors ! The name of the active vectors array, if any.
character(len=*), intent(in), optional :: Normals ! The name of the active normals array, if any.
character(len=*), intent(in), optional :: Tensors ! The name of the active tensors array, if any.
character(len=*), intent(in), optional :: TCoords ! The name of the active texture coordinates array, if any.

if (vtk_dat_xml('cell', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) stop

end subroutine

!-----------------------------------------------------------------------
! VTU_write_celldata: writes a celldata field  
!-----------------------------------------------------------------------
subroutine VTU_write_celldata(field, namef, typef)

real(R8P),    dimension(:)   :: field !field
character(len=*), intent(in) :: namef              !field name (ex., temperature)
character(len=*), intent(in) :: typef              !field type (scalar/vector)
integer :: i

print'(a)', 'Saving celldata field '//trim(namef)//'...'
select case(trim(lcase(typef)))
case('scalar')
  if (vtk_var_xml(nl, trim(namef), field) /= 0) stop
case('vector')
  select case(DIM)
  case (2)      
    if (vtk_var_xml(nl, trim(namef), field(1:nl*DIM:DIM), field(2:nl*DIM:DIM), &
    (/(real(0,R8P),i=1,nl)/)) /= 0) stop
  case (3)      
    if (vtk_var_xml(nl, trim(namef), field(1:nl*DIM:DIM), field(2:nl*DIM:DIM), &
    field(3:nl*DIM:DIM)) /= 0) stop
  case default; stop 'Not implemented'
  end select
case default; stop 'Not implemented'
end select

end subroutine

!-----------------------------------------------------------------------
! VTU_end_celldata: stops the writing of cellata fields
!-----------------------------------------------------------------------
subroutine VTU_end_celldata()

if (vtk_dat_xml('cell', 'CLOSE')/=0) stop

end subroutine

!-----------------------------------------------------------------------
! VTU_close: closes a vtu file
!-----------------------------------------------------------------------
subroutine VTU_close(geo)
  logical,optional ::geo

if(present(geo)) then
  if(geo) then; if (vtk_geo_xml()/=0) stop; endif
else
  if (vtk_geo_xml()/=0) stop
endif
print'(a)', 'Saving appended data...'
if (vtk_end_xml()/=0) stop

end subroutine

!-----------------------------------------------------------------------
! Write a RECONVXX mesh (and a B field) into a VTU binary file
! Writes reals in double precision (8 bytes) - R8P
!-----------------------------------------------------------------------
subroutine writeVTU(nel, nnod, nn, z, celldesc, field, namef, typef, DOFType, ficho, &
Scalars, Vectors, Normals, Tensors, TCoords)
integer(I4P), intent(in) :: nel                    !global number of elements
integer(I4P), intent(in) :: nnod                   !global number of nodes
integer(I4P), dimension(:,:) :: nn    !nodes index array
real(R8P),    dimension(:,:) :: z     !vertices coordinates array
character(len=*), intent(in) :: celldesc           !cell description
real(R8P),    dimension(:)   :: field !field
character(len=*), intent(in) :: namef              !field name (ex., temperature)
character(len=*), intent(in) :: typef              !field type (scalar/vector)
character(len=*), intent(in) :: DOFtype            !field DOF type (node/cell)
character(len=*), intent(in) :: ficho              !output file
character(len=*), intent(in), optional :: Scalars ! The name of the active scalars array, if any.
character(len=*), intent(in), optional :: Vectors ! The name of the active vectors array, if any.
character(len=*), intent(in), optional :: Normals ! The name of the active normals array, if any.
character(len=*), intent(in), optional :: Tensors ! The name of the active tensors array, if any.
character(len=*), intent(in), optional :: TCoords ! The name of the active texture coordinates array, if any.

integer      :: DIM   !space dimension
integer      :: LNV   !local number of vertices
integer      :: LNE   !local number of edges
integer      :: LNF   !local number of faces
integer      :: LNN   !local number of nodes
integer(I1P) :: CTYPE !cell linear type (according to VTK format)
integer :: i

integer(I4P), dimension(:), allocatable  :: nn_aux ! fix ifort
!print*, '----------', 2
if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) stop 'writeVTU: cell_type error'
print'(a)', 'Writing a VTU binary file...'
print'(a)', 'Saving mesh data...'

!Create output
if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) stop 'writeVTU: vtk_ini_xml error'
!saving coords in the current piece
if (DIM < 3) then
  if (vtk_geo_xml(nnod, nel, z(1,:), z(2,:), (/(real(0,R8P),i=1,nnod)/)) /= 0) & 
   stop 'writeVTU: vtk_geo_xml error 1'
else
  if (vtk_geo_xml(nnod, nel, z(1,:), z(2,:), z(3,:)) /= 0) & 
   stop 'writeVTU: vtk_geo_xml error 2'
end if

!saving connectivity
  allocate(nn_aux(LNN*nel)) ! fix ifort
  nn_aux = reshape(nn-1, (/LNN*nel/)) ! fix ifort
! esto cascaba con ifort e mallachip, en reshape: [sen heap arrays 0]
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) & 
!   stop 'writeVTU: vtk_con_xml error'
  if (vtk_con_xml(nel, nn_aux, (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) &
    stop 'writeVTU: vtk_con_xml error' ! fix ifort
  deallocate(nn_aux) ! fix ifort

print'(a)', 'Saving field data...'
!Node Fields 
select case(trim(lcase(DOFType)))
case('node')
  if (vtk_dat_xml('node', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) stop
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nnod, trim(namef), field) /= 0) stop
  case('vector')
    select case(DIM)
    case (2)      
      if (vtk_var_xml(nnod, trim(namef), field(1:nnod*DIM:DIM), field(2:nnod*DIM:DIM), &
      (/(real(0,R8P),i=1,nnod)/)) /= 0) stop
    case (3)      
      if (vtk_var_xml(nnod, trim(namef), field(1:nnod*DIM:DIM), field(2:nnod*DIM:DIM), &
      field(3:nnod*DIM:DIM)) /= 0) stop
    case default; stop 'Not implemented (node)'
    end select
  case default; stop 'Not implemented (node)'
  end select
  if (vtk_dat_xml('node', 'CLOSE')/=0) stop
case('cell')
!Cell Fields 
  if (vtk_dat_xml('cell', 'OPEN', Scalars, Vectors, Normals, Tensors, TCoords)/=0) stop
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nel, trim(namef), field) /= 0) stop
  case('vector')
    select case(DIM)
    case (2)      
      if (vtk_var_xml(nel, trim(namef), field(1:nel*DIM:DIM), field(2:nel*DIM:DIM), &
      (/(real(0,R8P),i=1,nel)/)) /= 0) stop
    case (3)      
      if (vtk_var_xml(nel, trim(namef), field(1:nel*DIM:DIM), field(2:nel*DIM:DIM), &
      field(3:nel*DIM:DIM)) /= 0) stop
    case default; stop 'Not implemented (cell)'
    end select
  case default;  stop 'Not implemented (cell)'
  end select
  if (vtk_dat_xml('cell', 'CLOSE')/=0) stop
case default; stop 'Not implemented (close)'
end select
if (vtk_geo_xml()/=0) stop
print'(a)', 'Saving appended data...'
if (vtk_end_xml()/=0) stop
print'(a)', 'Done!'

end subroutine

!-----------------------------------------------------------------------
! Write a RECONVXX mesh (and a B field) into a VTU binary file
! Modified to write reals in single precision (4 bytes) - R4P
!-----------------------------------------------------------------------
subroutine writeVTU4(nel, nnod, nn, z, celldesc, field, namef, typef, DOFType, ficho)

integer(I4P), intent(in) :: nel                    !global number of elements
integer(I4P), intent(in) :: nnod                   !global number of nodes
integer(I4P), dimension(:,:) :: nn    !nodes index array
real(R8P),    dimension(:,:) :: z     !vertices coordinates array
character(len=*), intent(in) :: celldesc           !cell description
real(R8P),    dimension(:)   :: field !field
character(len=*), intent(in) :: namef              !field name (ex., temperature)
character(len=*), intent(in) :: typef              !field type (scalar/vector)
character(len=*), intent(in) :: DOFtype            !field DOF type (node/cell)
character(len=*), intent(in) :: ficho              !output file

integer      :: DIM   !space dimension
integer      :: LNV   !local number of vertices
integer      :: LNE   !local number of edges
integer      :: LNF   !local number of faces
integer      :: LNN   !local number of nodes
integer(I1P) :: CTYPE !cell linear type (according to VTK format)
integer :: i

integer(I4P), dimension(:), allocatable  :: nn_aux ! fix ifort
!print*, '----------', 3
if (cell_type(celldesc, DIM, LNV, LNN, LNE, LNF, CTYPE) /= 0) stop 'writeVTU: cell_type error'
print'(a)', 'Writing a VTU binary file...'
print'(a)', 'Saving mesh data...'

!Create output
if (vtk_ini_xml('BINARY', ficho, 'UnstructuredGrid') /= 0) stop 'writeVTU: vtk_ini_xml error'
!saving coords in the current piece
if (DIM < 3) then
  if (vtk_geo_xml(nnod, nel, real(z(1,:),R4P), real(z(2,:),R4P), (/(real(0,R4P),i=1,nnod)/)) /= 0)&
     & stop 'writeVTU: vtk_geo_xml error 1'
else
  if (vtk_geo_xml(nnod, nel, real(z(1,:),R4P), real(z(2,:),R4P), real(z(3,:),R4P)) /= 0) & 
   stop 'writeVTU: vtk_geo_xml error 2'
end if

!saving connectivity
  allocate(nn_aux(LNN*nel)) ! fix ifort
  nn_aux = reshape(nn-1, (/LNN*nel/)) ! fix ifort
! esto poderia cascar con ifort e mallachip, en reshape: [sen heap arrays 0]
! if (vtk_con_xml(nel, reshape(nn-1, (/LNN*nel/)), (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0)  &
!   stop 'writeVTU: vtk_con_xml error'
  if (vtk_con_xml(nel, nn_aux, (/(i*LNN,i=1,nel)/), (/(CTYPE,i=1,nel)/)) /= 0) & 
    stop 'writeVTU: vtk_con_xml error' ! fix ifort
  deallocate(nn_aux) ! fix ifort

print'(a)', 'Saving field data...'
!Node Fields 
select case(trim(lcase(DOFType)))
case('node')
  if (vtk_dat_xml('node', 'OPEN')/=0) stop
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nnod, trim(namef), real(field,R4P)) /= 0) stop
  case('vector')
    select case(DIM)
    case (2)      
      if (vtk_var_xml(nnod, trim(namef),&
        &real(field(1:nnod*DIM:DIM),R4P),&
        &real(field(2:nnod*DIM:DIM),R4P),&
        &(/(real(0,R4P),i=1,nnod)/)) /= 0) stop
    case (3)      
      if (vtk_var_xml(nnod, trim(namef),&
        &real(field(1:nnod*DIM:DIM),R4P),&
        &real(field(2:nnod*DIM:DIM),R4P),&
        &real(field(3:nnod*DIM:DIM),R4P)) /= 0) stop
    case default; stop 'Not implemented'
    end select
  case default; stop 'Not implemented'
  end select
  if (vtk_dat_xml('node', 'CLOSE')/=0) stop
case('cell')
!Cell Fields 
  if (vtk_dat_xml('cell', 'OPEN')/=0) stop
  select case(trim(lcase(typef)))
  case('scalar')
    if (vtk_var_xml(nel, trim(namef), real(field,R4P)) /= 0) stop
  case('vector')
    select case(DIM)
    case (2)      
      if (vtk_var_xml(nel, trim(namef),&
        &real(field(1:nel*DIM:DIM),R4P),&
        &real(field(2:nel*DIM:DIM),R4P),&
        &(/(real(0,R4P),i=1,nel)/)) /= 0) stop
    case (3)      
      if (vtk_var_xml(nel, trim(namef),&
        &real(field(1:nel*DIM:DIM),R4P),&
        &real(field(2:nel*DIM:DIM),R4P),&
        &real(field(3:nel*DIM:DIM),R4P)) /= 0) stop
    case default; stop 'Not implemented'
    end select
  case default; stop 'Not implemented'
  end select
  if (vtk_dat_xml('cell', 'CLOSE')/=0) stop
case default; stop 'Not implemented'
end select
if (vtk_geo_xml()/=0) stop
print'(a)', 'Saving appended data...'
if (vtk_end_xml()/=0) stop
print'(a)', 'Done!'

end subroutine

!-----------------------------------------------------------------------
! cell_type: give the associated dimensions of a FE
! RESULT: 0 if the FE is known, -1 otherwise
!-----------------------------------------------------------------------
function cell_type(name, DIM, LNV, LNN, LNE, LNF, CTYPE) result(res)

character(len=*), intent(in)  :: name  !FE name
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

!print*, '-------------------------------------------------'
!print*, name, DIM, LNV, LNN, LNE, LNF, CTYPE
!print*, '-------------------------------------------------'

end function

!-----------------------------------------------------------------------
! lcase: converts a string to lower case
!-----------------------------------------------------------------------
function lcase(str) result(res)
 
  character(len=*), intent(in) :: str
  character(len=len(str)) :: res
  integer :: diff, i

  res = str
  diff = ichar('A') - ichar('a')
  do i = 1, len_trim(str)
    if (str(i:i) < 'A' .or. str(i:i) > 'Z') cycle
    res(i:i) = char(ichar(str(i:i)) - diff)
  end do
 
end function

end module

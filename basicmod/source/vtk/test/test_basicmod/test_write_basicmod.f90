program test_writevtu
!use iso_fortran_env, only: real64
use basicmod, only: real64, VTU_open, VTU_write_mesh, VTU_begin_pointdata, VTU_write_pointdata, VTU_end_pointdata, & 
VTU_begin_celldata, VTU_write_celldata, VTU_end_celldata, VTU_close
implicit none
integer :: i
integer :: nel = 2 ! Number of cells
integer :: nver = 5 ! Number of points
integer, allocatable :: mm(:,:) ! Tetrahedron connectivities, 4 x nel
real(real64), allocatable :: z(:,:) ! Point coordinates, 3 x nver
real(real64), allocatable :: field(:) ! field array

allocate(mm(4,2)) ! Tetrahedron connectivities, 4 x nel
mm(:,1) = [1,2,3,4]; mm(:,2) = [1,3,2,5]
allocate(z(3,5)) ! Point coordinates, 3 x nver
z(:,1) = [ 0., 0., 0.]
z(:,2) = [ 1., 0., 0.]
z(:,3) = [ 0., 1., 0.]
z(:,4) = [0.3, 0.3, 1.]
z(:,5) = [0.3, 0.3, -1.]
call VTU_open('test.vtu') ! Open file
call VTU_write_mesh(nel, nver, mm, z, 'tetra') ! Write mesh
call VTU_begin_pointdata(Scalars='pd2', Vectors='pd4') ! Begin pointdata record
allocate(field(5))
field = [(real(1,8)*i, i = 1,5)]
call VTU_write_pointdata(field, 'pd1', 'scalar') ! Write pointdata scalar field
field = [(real(1,8)*i, i = 5,1,-1)]
call VTU_write_pointdata(field, 'pd2', 'scalar') ! Write pointdata scalar field
deallocate(field); allocate(field(15))
field = [(real(1,8)*i, i = 1,15)]
call VTU_write_pointdata(field, 'pd3', 'vector') ! Write pointdata vector field
field = [(real(1,8)*i, i = 15,1,-1)]
call VTU_write_pointdata(field, 'pd4', 'vector') ! Write pointdata vector field
call VTU_end_pointdata() ! End pointdata record
call VTU_begin_celldata(Scalars='pd5', Vectors='pd7') ! Begin celldata record
deallocate(field); allocate(field(2))
field = [(real(1,8)*i, i = 1,2)]
call VTU_write_celldata(field, 'pd5', 'scalar') ! Write celldata scalar field
field = [(real(1,8)*i, i = 2,1,-1)]
call VTU_write_celldata(field, 'pd6', 'scalar') ! Write celldata scalar field
deallocate(field); allocate(field(6))
field = [(real(1,8)*i, i = 1,6)]
call VTU_write_celldata(field, 'pd7', 'vector') ! Write celldata vector field
field = [(real(1,8)*i, i = -1,-6,-1)]
call VTU_write_celldata(field, 'pd8', 'vector') ! Write celldata vector field
call VTU_end_celldata() ! End celldata record
call VTU_close() ! Close file
end program 
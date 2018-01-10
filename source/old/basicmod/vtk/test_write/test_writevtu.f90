program test_writeVTU
!Procedimiento para crear desde fortran un fichero vtu con varios campos
!1) Apertura del fichero
!call VTU_open('proba.vtu')
!2) Escritura de la malla (en este ejemplo, de tetraedros)
!call VTU_write_mesh(nel, nver, mm, z, 'tetra')
!3) En caso de que haya campos por vértices, se abre el bloque pointdata
!   Se puede indicar en los argumentos opcionales Scalars y Vectors el nombre de los campos
!   que un visualizador como Paraview abrirá por defecto
!call VTU_begin_pointdata(Scalars='pd2', Vectors='pd4')
!4) Escritura de los campos pointdata, uno detrás de otro (el orden no importa)
!call VTU_write_pointdata(field1, 'pd1', 'scalar')
!call VTU_write_pointdata(field2, 'pd2', 'scalar')
!call VTU_write_pointdata(field3, 'pd3', 'vector')
!call VTU_write_pointdata(field4, 'pd4', 'vector')
!6) En caso de que se hayan escrito campos por vértices, se cierra el bloque pointdata
!call VTU_end_pointdata()
!7) En caso de que haya campos por elementos, se abre el bloque celldata
!   Se puede indicar en los argumentos opcionales Scalars y Vectors el nombre de los campos
!   que un visualizador como Paraview abrirá por defecto
!call VTU_begin_celldata(Scalars='pd5', Vectors='pd7')
!8) Escritura de los campos celldata, uno detrás de otro (el orden no importa)
!call VTU_write_celldata(field5, 'pd5', 'scalar')
!call VTU_write_celldata(field6, 'pd6', 'scalar')
!call VTU_write_celldata(field7, 'pd7', 'vector')
!call VTU_write_celldata(field8, 'pd8', 'vector')
!9) En caso de que se hayan escrito campos por vértices, se cierra el bloque pointdata
!call VTU_end_celldata()
!10) Finalmente, se cierra el fichero 
!call VTU_close()

use module_writeVTU
implicit none

integer :: nel = 2
integer :: nver = 5
integer, allocatable, dimension(:,:) :: mm 
real(8), allocatable, dimension(:,:) :: z
real(8), allocatable, dimension(:)   :: field 
integer :: i

allocate(mm(4,2))
mm(:,1) = [1,2,3,4]; mm(:,2) = [1,3,2,5]

allocate(z(3,5))
z(:,1) = [ 0.,  0.,  0.]
z(:,2) = [ 1.,  0.,  0.] 
z(:,3) = [ 0.,  1.,  0.] 
z(:,4) = [0.3, 0.3,  1.] 
z(:,5) = [0.3, 0.3, -1.] 
!mesh
call VTU_open('proba.vtu')
call VTU_write_mesh(nel, nver, mm, z, 'tetra')

!pointdata
call VTU_begin_pointdata(Scalars='pd2', Vectors='pd4')
allocate(field(5))
field = [(real(1,8)*i, i = 1,5)];    call VTU_write_pointdata(field, 'pd1', 'scalar')
field = [(real(1,8)*i, i = 5,1,-1)]; call VTU_write_pointdata(field, 'pd2', 'scalar')

deallocate(field); allocate(field(15))
field = [(real(1,8)*i, i = 1,15)];    call VTU_write_pointdata(field, 'pd3', 'vector')
field = [(real(1,8)*i, i = 15,1,-1)]; call VTU_write_pointdata(field, 'pd4', 'vector')
call VTU_end_pointdata()

!celldata
call VTU_begin_celldata(Scalars='pd5', Vectors='pd7')
deallocate(field); allocate(field(2))
field = [(real(1,8)*i, i = 1,2)];    call VTU_write_celldata(field, 'pd5', 'scalar')
field = [(real(1,8)*i, i = 2,1,-1)]; call VTU_write_celldata(field, 'pd6', 'scalar')

deallocate(field); allocate(field(6))
field = [(real(1,8)*i, i = 1,6)];    call VTU_write_celldata(field, 'pd7', 'vector')
field = [(real(1,8)*i, i = -1,-6,-1)]; call VTU_write_celldata(field, 'pd8', 'vector')
call VTU_end_celldata()

call VTU_close()

end program      

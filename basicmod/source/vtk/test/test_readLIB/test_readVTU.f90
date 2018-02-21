program prueba
use LIB_VTK_IO
use LIB_VTK_IO_READ
implicit none
integer(I4P) :: np = 5, nn, nc, nco
real(R8P), allocatable :: X(:), Y(:), Z(:), var(:)
integer(I4P), allocatable :: connect(:), offset(:)
integer(I1P), allocatable :: cell_type(:)


if (vtk_ini_xml_read('Binary','proba.vtu','UnstructuredGrid', np)/=0) stop 'Error'
print*, np

if (vtk_geo_xml_read(nn, nc, X, Y, Z, 1_I4P)/=0) stop 'Error'
print*, nn, nc,X,Y,Z

if (vtk_con_xml_read(nc,connect,offset,cell_type)/=0) stop 'Error'
print*, nc,connect,offset,cell_type

if(vtk_var_xml_read('node',nn,nco,'pd1',var)/=0) stop 'Error'
print*, nn,nco,var

if(vtk_var_xml_read('cell',nn,nco,'pd8',var)/=0) stop 'Error'
print*, nn,nco,var

if (vtk_end_xml_read()/=0) stop 'Error'

end program


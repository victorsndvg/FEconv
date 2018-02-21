program test_readvtu
use module_readvtu, only: real64, VTU_read_open, VTU_read_mesh, VTU_read_data, VTU_read_close
implicit none
integer                   :: nel       ! Global number of elements
integer                   :: nnod      ! Global number of nodes
integer,      allocatable :: nn(:,:)   ! Mesh connectivities
real(real64), allocatable :: z(:,:)    ! Node coordinates
character(20)             :: celldesc  ! Cell description (see VTU_read_mesh documentation)
real(real64), allocatable :: nfield(:) ! Node field
real(real64), allocatable :: cfield(:) ! Cell field
character(20)             :: typef     ! Field type (see VTU_read_pointdata_documentation)

call VTU_read_open('test.vtu')
call VTU_read_mesh(nel, nnod, nn, z, celldesc)
print*, 'Cell description: '//trim(celldesc)
if  (VTU_read_data(nfield, 'pd4', 'node', typef) /= 0) print*, 'Field pd4 not found.'
print*, 'Field pd4: '//trim(typef)
if  (VTU_read_data(nfield, 'XXX', 'node', typef) /= 0) print*, 'Field XXX not found.'
if  (VTU_read_data(cfield, 'pd7', 'cell', typef) /= 0) print*, 'Field pd7 not found.'
print*, 'Field pd7: '//trim(typef)
if  (VTU_read_data(cfield, 'XXX', 'cell', typef) /= 0) print*, 'Field XXX not found.'
call VTU_read_close()
end program
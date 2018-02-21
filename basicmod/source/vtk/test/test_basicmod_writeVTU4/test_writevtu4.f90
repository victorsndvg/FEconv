program test_writevtu
use basicmod, only: real64, writeVTU4
implicit none
integer :: nel = 2 ! Number of cells
integer :: nver = 5 ! Number of points
integer :: mm(4,2) ! Tetrahedron connectivities, 4 x nel
real(real64) :: z(3,5) ! Point coordinates, 3 x nver
real(real64) :: field(5) ! Pointdata field array
integer :: i
mm(:,1) = [1,2,3,4]; mm(:,2) = [1,3,2,5] ! Tetrahedron connectivities
z(:,1) = [ 0., 0., 0.] ! Point coordinates
z(:,2) = [ 1., 0., 0.]
z(:,3) = [ 0., 1., 0.]
z(:,4) = [0.3, 0.3, 1.]
z(:,5) = [0.3, 0.3, -1.]
call writeVTU4(nel, nver, mm, z, 'tetra', field, 'Potential (V)', 'scalar', 'node', 'poten.vtu')
end program
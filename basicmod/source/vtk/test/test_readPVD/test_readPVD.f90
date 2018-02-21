program test_readPVD
use module_os_dependant, only: maxpath
use module_readPVD, only: real64, readPVD
implicit none
integer :: i
real(real64), allocatable :: time(:)
character(maxpath), allocatable :: vtu(:)
call readPVD('temp_evolution.pvd', time, vtu)
print*, time
print*, (trim(vtu(i))//' ', i = 1, size(vtu,1)) 
end program

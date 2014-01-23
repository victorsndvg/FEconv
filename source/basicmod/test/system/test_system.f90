program test_system
use module_os_dependant, only: set_os, get_os
use module_system
implicit none
integer, allocatable :: v(:)
integer :: res, i
character(maxpath), allocatable :: list(:)
character(maxpath) :: o

print*,' '
print*, 'Set the operating system:'
print*, 'call set_os()'
         call set_os()

print*,' '
print*, 'Get the operating system:'
print*, 'o = get_os()'
         o = get_os()
print*, '(RESULT) o = ',trim(o)         
         
print*,' '
print*, 'List the files of folder:'
print*, 'res = ls(''.'', list)'
         res = ls('.', list)
print*, '(RESULT) res = ',res
do i = 1, size(list,1)
  print*, '(RESULT) list(',i,') = ', trim(list(i))
end do
 
print*,' '
print*, 'Make dir:'
print*, 'res = mkdir(''testing_folder'')'
         res = mkdir('testing_folder')
print*, '(RESULT) res = ',res

print*,' '
print*, 'Del file:'
print*, 'res = rm(''testing_file'')'
         res = rm('testing_file')
print*, '(RESULT) res = ',res         

end program   

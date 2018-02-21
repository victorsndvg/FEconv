program test_files
use module_compiler_dependant
use module_os_dependant
use module_files
implicit none

character(maxpath) :: errmsg
integer :: v1(8), v2(8), values(8)

print*, file_exists('a1')
print*, file_exists('a2')
v1 = modification_time('a1', errmsg) 
if (v1(1) <= 0) then 
  print*, 'a1: '//trim(errmsg)
else
  print*, modification_time('a1')
end if
v2 = modification_time('a2', errmsg) 
if (v2(1) <= 0) then 
  print*, 'a2: '//trim(errmsg)
else
  print*, modification_time('a2')
end if
print*, 'a1' .IsNewerThan. 'a2'


CALL DATE_AND_TIME(values=values)
print*, values
end program   

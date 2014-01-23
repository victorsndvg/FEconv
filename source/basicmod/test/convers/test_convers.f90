program test_convers
use module_convers
implicit none
character(200) :: res

print*,' '
print*, 'Trim:'
print*, 'res = trim(''testing.file.txt'',''.'')'
         res = trim('testing.file.txt','.')
print*, '(RESULT) res = ',res         

print*, 'res = trim(''testing.file.txt'',''.'', back=.true.)'
         res = trim('testing.file.txt','.', back=.true.)
print*, '(RESULT) res = ',res         

end program   

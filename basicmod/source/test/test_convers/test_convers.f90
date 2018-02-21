program test_convers
use module_convers
implicit none
character(200) :: res, cad

print*,' '
print*, 'Trim:'
print*, 'res = trim(''testing.file.txt'',''.'')'
         res = trim('testing.file.txt','.')
print*, '(RESULT) res = ', trim(res)         

print*, 'res = trim(''testing.file.txt'',''.'', back=.true.)'
         res = trim('testing.file.txt','.', back=.true.)
print*, '(RESULT) res = ', trim(res)         

print*,' '
!cad = '"u '' no" dos, tres'
cad = '"prue/ba espacio/test/mecano_cubo_n2293.msh" , dos,tres/cuatro'
print*, 'Word_count: ', "'"//trim(cad)//"'"
print*, '(RESULT)', word_count(cad)
print*, 'Word 1: ', word(cad, 1)
print*, 'Word 2: ', word(cad, 2)
print*, 'Word 3: ', word(cad, 3)
print*, 'Word 4: ', word(cad, 4)
end program   

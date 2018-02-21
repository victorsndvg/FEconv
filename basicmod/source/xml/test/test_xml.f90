program test
use module_report, only: report_option
use module_xml_parser, only: real64, fopen, fread, fread_alloc, flist, fclose
implicit none
integer :: ide
real(real64) :: cval, inter(5)
real(real64), allocatable :: ra(:) 
character(20), allocatable :: list(:) 
character(20) :: cc, ca(5) = ['. ', '. ', '. ', '. ', '. ']

!call report_option('info', 'std')
ide = fopen()
!call flist(ide, '/Boundary conditions/', list)
!print*, list
call fread(ide, '/Boundary conditions/Neumann Condition/Constant value', cval)
print*, cval
call fread(ide, '/Boundary conditions/Neumann Condition/Interval', inter)
print*, inter
call fread_alloc(ide, '/Boundary conditions/Neumann Condition/Interval', ra, realloc=.true.)
print*, allocated(ra), size(ra,1), ra
call fread(ide, '/Boundary conditions/Neumann Condition/Label', cc)
print*, cc
call fread(ide, '/Boundary conditions/Neumann Condition/Surfaces', ca)
print*, ca
call fclose(ide)
end program
